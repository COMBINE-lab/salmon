/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Salmon.

    Salmon is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Salmon is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Salmon.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/


#include <iostream>
#include <fstream>
#include <istream>
#include <vector>
#include <array>
#include <unordered_map>
#include <limits>
#include <cmath>
#include <cstdint>
#include <numeric>

#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "Eigen/Dense"
#include "PCA.hpp"

#include "tensemble/TypeDef.h"
#include "tensemble/RandomForestRegressor.h"
#include "tensemble/RandomForestClassifier.h"
#include "tensemble/GBMRegressor.h"
#include "tensemble/GBMClassifier.h"
#include "tensemble/ReadData.h"

#include "CommonTypes.hpp"

#define DEFAULT_N_TREES 100
#define DEFAULT_N_JOBS 1
#define DEFAULT_MAX_FEATURES_RATIO 1.0
#define DEFAULT_MIN_SAMPLE_LEAF 5
#define DEFAULT_MAX_DEPTH   4
#define DEFAULT_SUBSAMPLE   1.0
#define DEFAULT_SPLIT_CRITERION CRITERION_MSE
#define DEFAULT_LOSS SQUARE_LOSS
#define DEFAULT_LEARN_RATE 0.1
#define DEFAULT_OOB 1
#define DEFAULT_VERBOSE 0
#define DEFAULT_BOOTSTRAP 1
#define DEFAULT_COMPUTE_IMPORTANCE 0

namespace bfs = boost::filesystem;
using Kmer = ::uint64_t;
using Sailfish::TranscriptFeatures;
using mpdec = boost::multiprecision::cpp_dec_float_100;

TranscriptFeatures parseFeature(std::ifstream& ifs) {
        TranscriptFeatures tf{};
        ifs >> tf.name;
        ifs >> tf.length;
        ifs >> tf.gcContent;
        for (auto i : boost::irange(size_t{0}, tf.diNucleotides.size())) {
                ifs >> tf.diNucleotides[i];
        }
        // eat the newline
        char junk;
        ifs.get(junk);
        return tf;
}

std::vector<TranscriptFeatures> parseFeatureFile(const bfs::path& featureFile) {
        std::ifstream ifile(featureFile.string());
        std::vector<TranscriptFeatures> feats;
        while (!ifile.eof()) {
                feats.emplace_back( parseFeature(ifile) );
                if (ifile.peek() == EOF) { break; }
        }
        ifile.close();
        return feats;
}

struct TranscriptResult{
        size_t length;
        double tpm;
        double rpkm;
        double kpkm;
        double approxKmerCount;
        double approxCount;
};

struct ExpressionResults {
        std::vector<std::string> comments;
        std::unordered_map<std::string, TranscriptResult> expressions;
};

ExpressionResults parseSalmonFile(const bfs::path& expFile, double& numMappedReads) {
        numMappedReads = 0.0;

        std::ifstream ifile(expFile.string());
        ExpressionResults res;
        while(!ifile.eof()) {

                if (ifile.peek() == '#') {
                        std::string comment;
                        std::getline(ifile, comment);
                        res.comments.emplace_back(comment);
                } else {
                        std::string tname;
                        TranscriptResult tr;
                        ifile >> tname;
                        ifile >> tr.length;
                        ifile >> tr.tpm;
                        ifile >> tr.rpkm;
                        ifile >> tr.approxCount;
                        numMappedReads += tr.approxCount;
                        res.expressions[tname] = tr;
                        // eat the newline
                        char nline; ifile.get(nline);
                }

                if (ifile.peek() == EOF) { break; }
        }

        return res;
}



ExpressionResults parseSailfishFile(const bfs::path& expFile) {
        std::ifstream ifile(expFile.string());
        ExpressionResults res;
        while(!ifile.eof()) {

                if (ifile.peek() == '#') {
                        std::string comment;
                        std::getline(ifile, comment);
                        res.comments.emplace_back(comment);
                } else {
                        std::string tname;
                        TranscriptResult tr;
                        ifile >> tname;
                        ifile >> tr.length;
                        ifile >> tr.tpm;
                        ifile >> tr.rpkm;
                        ifile >> tr.kpkm;
                        ifile >> tr.approxKmerCount;
                        ifile >> tr.approxCount;
                        res.expressions[tname] = tr;
                        // eat the newline
                        char nline; ifile.get(nline);
                }

                if (ifile.peek() == EOF) { break; }
        }

        return res;
}


void populateFromTPMs(vector<mpdec>& tpms,
                      vector<TranscriptFeatures>& features,
                      vector<size_t>& retainedRows,
                      ExpressionResults& sfres,
                      double estimatedReadLength,
                      double kmersPerRead,
                      uint64_t mappedKmers,
                      uint64_t merLen,
                      vector<mpdec>& rpkms,
                      vector<mpdec>& kmerCounts,
                      vector<mpdec>& readCounts) {

  // compute the TPM normalization factor
  mpdec mpzero = 0;
  mpdec sumTPM = std::accumulate(tpms.begin(), tpms.end(), mpzero);
  mpdec norm = 1.0 / sumTPM;

  // compute the relative transcript fractions
  vector<mpdec> tfracs(tpms.size());
  for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
    tfracs[i] = tpms[i] * norm;
  }

  // using the relative transcript fractions, compute the relative
  // nucleotide fractions (transcript fractions * length)
  vector<mpdec> tflens(tpms.size());
  for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
    auto& name = features[i].name;
    auto& r = sfres.expressions[name];
    double l = (r.length - merLen + 1);
    tflens[i] = tfracs[i] * (l);
  }

  // normalize the nucleotide fractions and fill in the estimated k-mer counts
  kmerCounts.clear(); kmerCounts.resize(features.size());
  mpdec tfnorm = 1.0 / std::accumulate(tflens.begin(), tflens.end(), mpzero);
  for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
    tflens[i] *= tfnorm;
    kmerCounts[i] = tflens[i] * mappedKmers;
  }

  uint64_t numReads = mappedKmers / kmersPerRead;

  // use the nucleotide fractions to compute the RPKMs
  double billion = pow(10,9);
  rpkms.clear(); rpkms.resize(features.size());
  readCounts.clear(); readCounts.resize(features.size());

  for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
    auto& name = features[i].name;
    auto& r = sfres.expressions[name];
    double l = (r.length - merLen + 1);
    rpkms[i] = billion * (tflens[i] / l);
    readCounts[i] = (tflens[i] * numReads);
  }

}


int performBiasCorrection(
        bfs::path featureFile,
        bfs::path expressionFile,
        double estimatedReadLength,
        double kmersPerRead,
        uint64_t mappedKmers,
        uint32_t merLen,
        bfs::path outputFile,
        size_t numThreads) {

        auto features = parseFeatureFile(featureFile);
        std::cerr << "parsed " << features.size() << " features\n";

        auto sfres = parseSailfishFile(expressionFile);
        std::cerr << "parsed " << sfres.expressions.size() << " expression values\n";

        std::vector<size_t> retainedRows;
        std::vector<double> retainedRPKMs;
        std::vector<std::string> retainedNames;

        double minLRPKM, maxLRPKM;
        minLRPKM = std::numeric_limits<double>::max();
        maxLRPKM = -minLRPKM;

        for (auto i : boost::irange(size_t{0}, features.size())) {
                auto& tname = features[i].name;
                auto rpkm = sfres.expressions[tname].kpkm; // ALTERATION
                double v;

                if ( rpkm >= 1e-3 ) {
                        retainedRows.emplace_back(i);
                        retainedNames.push_back(tname);
                        v = std::log(rpkm);
                        retainedRPKMs.push_back(v);
                        minLRPKM = std::min(minLRPKM, v);
                        maxLRPKM = std::max(maxLRPKM, v);
                }
        }

        Eigen::MatrixXd featMat(retainedRows.size(), 17);
        std::vector<float> pcavec;
        size_t fnum = 0;
        for (auto r : retainedRows) {
                auto& f = features[r];
                featMat(fnum, 0) = f.gcContent;
                pcavec.push_back(f.gcContent);
                for (auto i : boost::irange(size_t{0}, f.diNucleotides.size())) {
                        featMat(fnum, i+1) =  f.diNucleotides[i];
                        pcavec.push_back(f.diNucleotides[i]);
                }
                ++fnum;
        }

        PCA pca(featMat);

        std::cerr << "Performing PCA decomposition\n";
        pca.performDecomposition();

        auto encodedXSub = pca.projectedData(0.95, true);

        Data train;
        size_t numCols = encodedXSub.cols();
        train.set_size(retainedRows.size(), numCols+1);

        size_t c = 0;
        for (auto r : retainedRows) {
                train.X[c][0] = std::log(static_cast<double>(features[r].length));

                for (auto j : boost::irange(size_t{1}, numCols+1)) {
                        train.X[c][j] = encodedXSub(c, j-1);
                        //train.X[c][j] = scores[le];
                        //train2.X[c][j] = scores[le];
                        //++le;
                        //std::cerr << "Train [" << c <<"][" << j        << "] = " << train.X[c][j] << "\n";
                }
                //le += (17 - dimCutoff);
                train.y[c] = retainedRPKMs[c];
                ++c;
        }

        /** Random Forest Regression **/
        size_t minDepth = 5;
        auto reg = std::unique_ptr<RandomForestRegressor>(new RandomForestRegressor(
                500,
                train.n_features,
                5, // max tree depth
                1, // min_samples_leaf
                1.0, // features ratio
                true, // bootstrap
                true, //out-of-bag
                true, // compute importance
                0, // random seed
                numThreads, // num jobs
                true // verbose
        ));

        std::cerr << "there are " << train.n_samples << " samples\n";
        std::cerr << "there are " << train.n_features << " features\n";
        reg->build(train.X, train.y, train.n_samples);

        std::vector<REAL> pred(train.n_samples, 0.0);
        reg->predict(train.X, &pred[0], train.n_samples, train.n_features);

        REAL trn_rmse=rmse(&pred[0], train.y, train.n_samples);
        REAL trn_r2=R2(&pred[0], train.y, train.n_samples);
        std::cerr << "Train RMSE=" << trn_rmse << ", Correlation Coefficient=" << trn_r2 << "\n";

        double grandMean = 0.0;
        size_t ntrain = train.n_samples;
        for (auto i : boost::irange(size_t{0}, ntrain)) {
                grandMean += retainedRPKMs[i];
        }
        grandMean /= train.n_samples;

        for (auto i : boost::irange(size_t{0}, ntrain)) {
                pred[i] = grandMean + (retainedRPKMs[i] - pred[i]);
        }

        trn_rmse=rmse(&pred[0], train.y, train.n_samples);
        trn_r2=R2(&pred[0], train.y, train.n_samples);
        std::cerr << "Train RMSE=" << trn_rmse << ", Correlation Coefficient=" << trn_r2 << "\n";

        std::ofstream ofile(outputFile.string());
        for (auto& c : sfres.comments) {
                ofile << c << "\n";
        }


        size_t retainedCnt = 0;
        vector<mpdec> kpkms(features.size());
        for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
          auto& name = features[i].name;
          auto& r = sfres.expressions[name];
          if (i == retainedRows[retainedCnt]) {
            kpkms[i] = std::exp(pred[retainedCnt]);
            ++retainedCnt;
          } else {
              kpkms[i] = r.kpkm;
          }
        }


        // compute estimated TPM from the KPKMS
        mpdec mpzero = 0;
        // normalize the KPKMS --- these will estimate the tau_i
        mpdec sumKPKM = std::accumulate(kpkms.begin(), kpkms.end(), mpzero);
        mpdec norm = 1.0 / sumKPKM;

        // then multiply by 10^6 to get TPM_i
        double million = pow(10, 6);
        vector<mpdec> tpms(kpkms.size());
        for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
          tpms[i] = kpkms[i] * norm * million;
        }

        vector<mpdec> rpkms;
        vector<mpdec> kmerCounts;
        vector<mpdec> readCounts;
        populateFromTPMs(tpms, features, retainedRows,
                         sfres, estimatedReadLength, kmersPerRead,
                         mappedKmers, merLen, rpkms, kmerCounts, readCounts);

        for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
          auto& name = features[i].name;
          auto& r = sfres.expressions[name];
          auto length = r.length;
          double effectiveLength = length - merLen + 1;
          ofile << name << '\t' << r.length << '\t' << tpms[i] << '\t'
                << ((length - estimatedReadLength + 1) > 0 ? kpkms[i] : 0.0) << '\t'
                << ((length - merLen + 1) > 0 ? kpkms[i] : 0.0) << '\t'
                << kmerCounts[i] << '\t'
                << readCounts[i] << '\n';
        }
        std::cerr << "retainedCnt = " << retainedCnt << ", nsamps = " << train.n_samples << "\n";

    ofile.close();
    return 0;
}


int performBiasCorrectionSalmon(
        bfs::path featureFile,
        bfs::path expressionFile,
        bfs::path outputFile,
        size_t numThreads) {

        auto features = parseFeatureFile(featureFile);
        std::cerr << "parsed " << features.size() << " features\n";

        double numMappedReads = 0.0;
        auto salmonRes = parseSalmonFile(expressionFile, numMappedReads);
        std::cerr << "parsed " << salmonRes.expressions.size() << " expression values\n";

        std::vector<size_t> retainedRows;
        std::vector<double> retainedTPMs;
        std::vector<std::string> retainedNames;

        double minLTPM, maxLTPM;
        minLTPM = std::numeric_limits<double>::max();
        maxLTPM = -minLTPM;

        for (auto i : boost::irange(size_t{0}, features.size())) {
                auto& tname = features[i].name;
                auto tpm = salmonRes.expressions[tname].tpm; //
                double v;

                if ( tpm >= 1.0 ) {
                        retainedRows.emplace_back(i);
                        retainedNames.push_back(tname);
                        v = std::log(tpm);
                        retainedTPMs.push_back(v);
                        minLTPM = std::min(minLTPM, v);
                        maxLTPM = std::max(maxLTPM, v);
                }
        }


        std::vector<float> pcavec;
        std::vector<double> featData;
        Eigen::MatrixXd featMat(retainedRows.size(), 17);
        size_t fnum = 0;
        size_t linearIndex{0};
        for (auto r : retainedRows) {
                auto& f = features[r];
                pcavec.push_back(f.gcContent);
                featMat(fnum, 0) = f.gcContent;
                for (auto i : boost::irange(size_t{0}, f.diNucleotides.size())) {
                        featMat(fnum, i) = f.diNucleotides[i];
                        pcavec.push_back(f.diNucleotides[i]);
                }
                ++fnum;
        }

        PCA pca(featMat);

        std::cerr << "Performing PCA decomposition\n";
        pca.performDecomposition();

        auto encodedXSub = pca.projectedData(0.95, true);

        Data train;
        size_t numCols = encodedXSub.cols();
        train.set_size(retainedRows.size(), numCols+1);

        size_t c = 0;
        for (auto r : retainedRows) {
                train.X[c][0] = std::log(static_cast<double>(features[r].length));

                for (auto j : boost::irange(size_t{1}, numCols+1)) {
                        train.X[c][j] = encodedXSub(c, j-1);
                }
                train.y[c] = retainedTPMs[c];
                ++c;
        }

        /** Random Forest Regression **/
        size_t minDepth = 5;
        auto reg = std::unique_ptr<RandomForestRegressor>(new RandomForestRegressor(
                500,
                train.n_features,
                5, // max tree depth
                1, // min_samples_leaf
                1.0, // features ratio
                true, // bootstrap
                true, //out-of-bag
                true, // compute importance
                0, // random seed
                numThreads, // num jobs
                true // verbose
        ));

        std::cerr << "there are " << train.n_samples << " samples\n";
        std::cerr << "there are " << train.n_features << " features\n";
        reg->build(train.X, train.y, train.n_samples);

        std::vector<REAL> pred(train.n_samples, 0.0);
        reg->predict(train.X, &pred[0], train.n_samples, train.n_features);

        REAL trn_rmse=rmse(&pred[0], train.y, train.n_samples);
        REAL trn_r2=R2(&pred[0], train.y, train.n_samples);
        std::cerr << "Train RMSE=" << trn_rmse << ", Correlation Coefficient=" << trn_r2 << "\n";

        double grandMean = 0.0;
        size_t ntrain = train.n_samples;
        for (auto i : boost::irange(size_t{0}, ntrain)) {
                grandMean += retainedTPMs[i];
        }
        grandMean /= train.n_samples;

        for (auto i : boost::irange(size_t{0}, ntrain)) {
                pred[i] = grandMean + (retainedTPMs[i] - pred[i]);
        }

        trn_rmse=rmse(&pred[0], train.y, train.n_samples);
        trn_r2=R2(&pred[0], train.y, train.n_samples);
        std::cerr << "Train RMSE=" << trn_rmse << ", Correlation Coefficient=" << trn_r2 << "\n";

        std::ofstream ofile(outputFile.string());
        for (auto& c : salmonRes.comments) {
                ofile << c << "\n";
        }

        size_t retainedCnt = 0;
        vector<mpdec> tpms(features.size());
        double tpmSum{0.0};
        for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
          auto& name = features[i].name;
          auto& r = salmonRes.expressions[name];
          if (i == retainedRows[retainedCnt]) {
            double v = std::exp(pred[retainedCnt]);
            tpms[i] = v;
            tpmSum += v;
            ++retainedCnt;
          } else {
              tpms[i] = r.tpm;
              tpmSum += r.tpm;
          }
        }

        double tpmNorm = 1000000.0 / tpmSum;
        for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
            tpms[i] *= tpmNorm;
        }

        vector<mpdec> fpkms(features.size());
        vector<mpdec> estNumReads(features.size());

        // use TPM estimates to computed estimated read counts
        mpdec mpzero = 0;
        mpdec totalNucDenom = 0;
        for (auto i : boost::irange(size_t{0},  size_t{features.size()})) {
            double len = features[i].length;
            totalNucDenom += tpms[i] * len;
        }

        for (auto i : boost::irange(size_t{0},  size_t{features.size()})) {
            mpdec len = features[i].length;
            estNumReads[i] += ((tpms[i] * len) / totalNucDenom) * numMappedReads;
        }

        double oneBillion = 1000000000.0;
        for (auto i : boost::irange(size_t{0},  size_t{features.size()})) {
            double len = features[i].length;
            fpkms[i] = (estNumReads[i] / (len * numMappedReads)) * oneBillion;
        }

        for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
          auto& name = features[i].name;
          auto& r = salmonRes.expressions[name];
          auto length = r.length;
          ofile << name << '\t'
                << r.length << '\t'
                << tpms[i] << '\t'
                << fpkms[i] << '\t'
                << estNumReads[i] << '\n';
        }
        std::cerr << "retainedCnt = " << retainedCnt << ", nsamps = " << train.n_samples << "\n";

        ofile.close();
    return 0;
}



// GRAVEYARD --- Shark Version

/*int performBiasCorrection(*/
        //bfs::path featureFile,
        //bfs::path expressionFile,
        //double estimatedReadLength,
        //double kmersPerRead,
        //uint64_t mappedKmers,
        //uint32_t merLen,
        //bfs::path outputFile,
        //size_t numThreads) {
        ////int argc, char* argv[]) {

        //using shark::PCA;

        //auto features = parseFeatureFile(featureFile);
        //std::cerr << "parsed " << features.size() << " features\n";

        //auto sfres = parseSailfishFile(expressionFile);
        //std::cerr << "parsed " << sfres.expressions.size() << " expression values\n";

        //std::vector<size_t> retainedRows;
        //std::vector<double> retainedRPKMs;
        //std::vector<std::string> retainedNames;

        //double minLRPKM, maxLRPKM;
        //minLRPKM = std::numeric_limits<double>::max();
        //maxLRPKM = -minLRPKM;

        //for (auto i : boost::irange(size_t{0}, features.size())) {
                //auto& tname = features[i].name;
                //auto rpkm = sfres.expressions[tname].kpkm; // ALTERATION
                //shark::RealVector v(1);

                //if ( rpkm >= 1e-3 ) {
                        //retainedRows.emplace_back(i);
                        //retainedNames.push_back(tname);
                        //v(0) = std::log(rpkm);
                        //retainedRPKMs.push_back(v(0));
                        //minLRPKM = std::min(minLRPKM, v(0));
                        //maxLRPKM = std::max(maxLRPKM, v(0));
                //}
        //}


        ////Pca pca;

        //std::vector<shark::RealVector> featMat;
        //std::vector<float> pcavec;
        //size_t fnum = 0;
        //for (auto r : retainedRows) {
                //auto& f = features[r];
                //shark::RealVector v(17);
                //std::vector<double> fv(17);
                //fv[0] = f.gcContent;
                //v(0) = f.gcContent;
                //pcavec.push_back(f.gcContent);
                //for (auto i : boost::irange(size_t{0}, f.diNucleotides.size())) {
                        //v(i+1) = f.diNucleotides[i];
                        //fv[i+1] = f.diNucleotides[i];
                        //pcavec.push_back(f.diNucleotides[i]);
                //}
                //featMat.emplace_back(v);
                //++fnum;
        //}

        //shark::UnlabeledData<shark::RealVector> Xsub = shark::createDataFromRange(featMat);

        //PCA pcao(Xsub, true);
        //auto evals = pcao.eigenvalues();
        //double totalVariance = 0.0;
        //for ( auto e : evals ) { std::cerr << e << "\n"; totalVariance += e; }
        //std::cerr << "totalVariance: " << totalVariance << "\n";
        //double varCutoff = 0.95;
        //double varSum = 0.0;
        //size_t dimCutoff = 0;
        //size_t currentDim = 0;
        //for ( auto e : evals ) {
                //++currentDim;
                //varSum += e;
                //std::cerr << "ev: " << e <<  "\n";
                //if (varSum / totalVariance >= varCutoff) {
                        //dimCutoff = currentDim;
                        //break;
                //}

        //}
        //std::cerr << varCutoff * 100.0 << "% of the variance is explained by " << dimCutoff << " dimensions\n";
        //shark::LinearModel<> enc;
        //pcao.encoder(enc, dimCutoff);
        //auto encodedXsub = enc(Xsub);

        ////shark::UnlabeledData<shark::RealVector> X(18, features.size());
        //Data train;
        //train.set_size(retainedRows.size(), dimCutoff+1);

        ////size_t le = 0;
        //size_t c = 0;
        //for (auto r : retainedRows) {
                //train.X[c][0] = std::log(static_cast<double>(features[r].length));

                //for (auto j : boost::irange(size_t{1}, dimCutoff+1)) {
                        //train.X[c][j] = encodedXsub.element(c)(j-1);
                        ////train.X[c][j] = scores[le];
                        ////train2.X[c][j] = scores[le];
                        ////++le;
                        ////std::cerr << "Train [" << c <<"][" << j        << "] = " << train.X[c][j] << "\n";
                //}
                ////le += (17 - dimCutoff);
                //train.y[c] = retainedRPKMs[c];// std::log(sfres[features[r].name].expressions.rpkm);
                //++c;
        //}

        //size_t minDepth = 5;
        //auto reg = std::unique_ptr<RandomForestRegressor>(new RandomForestRegressor(
                //500,
                //train.n_features,
                //5, // max tree depth
                //1, // min_samples_leaf
                //1.0, // features ratio
                //true, // bootstrap
                //true, //out-of-bag
                //true, // compute importance
                //0, // random seed
                //numThreads, // num jobs
                //true // verbose
        //));


        //std::cerr << "there are " << train.n_samples << " samples\n";
        //std::cerr << "there are " << train.n_features << " features\n";
        //reg->build(train.X, train.y, train.n_samples);

        //std::vector<REAL> pred(train.n_samples, 0.0);
        //reg->predict(train.X, &pred[0], train.n_samples, train.n_features);

        //REAL trn_rmse=rmse(&pred[0], train.y, train.n_samples);
        //REAL trn_r2=R2(&pred[0], train.y, train.n_samples);
        //std::cerr << "Train RMSE=" << trn_rmse << ", Correlation Coefficient=" << trn_r2 << "\n";

        //double grandMean = 0.0;
        //size_t ntrain = train.n_samples;
        //for (auto i : boost::irange(size_t{0}, ntrain)) {
                //grandMean += retainedRPKMs[i];
        //}
        //grandMean /= train.n_samples;

        //for (auto i : boost::irange(size_t{0}, ntrain)) {
                //pred[i] = grandMean + (retainedRPKMs[i] - pred[i]);
        //}

   //trn_rmse=rmse(&pred[0], train.y, train.n_samples);
   //trn_r2=R2(&pred[0], train.y, train.n_samples);
  //std::cerr << "Train RMSE=" << trn_rmse << ", Correlation Coefficient=" << trn_r2 << "\n";


        ////shark::UnlabeledData<shark::RealVector> X()

        //std::ofstream ofile(outputFile.string());
        //for (auto& c : sfres.comments) {
                //ofile << c << "\n";
        //}


        //size_t retainedCnt = 0;
        //vector<mpdec> kpkms(features.size());
        //for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
          //auto& name = features[i].name;
          //auto& r = sfres.expressions[name];
          //if (i == retainedRows[retainedCnt]) {
            //kpkms[i] = std::exp(pred[retainedCnt]);
            //++retainedCnt;
          //} else {
              //kpkms[i] = r.kpkm;
          //}
        //}


        //// compute estimated TPM from the KPKMS
        //mpdec mpzero = 0;
        //// normalize the KPKMS --- these will estimate the tau_i
        //mpdec sumKPKM = std::accumulate(kpkms.begin(), kpkms.end(), mpzero);
        //mpdec norm = 1.0 / sumKPKM;

        //// then multiply by 10^6 to get TPM_i
        //double million = pow(10, 6);
        //vector<mpdec> tpms(kpkms.size());
        //for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
          //tpms[i] = kpkms[i] * norm * million;
        //}
        //vector<mpdec> rpkms;
        //vector<mpdec> kmerCounts;
        //vector<mpdec> readCounts;
        //populateFromTPMs(tpms, features, retainedRows,
                         //sfres, estimatedReadLength, kmersPerRead,
                         //mappedKmers, merLen, rpkms, kmerCounts, readCounts);

        //for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
          //auto& name = features[i].name;
          //auto& r = sfres.expressions[name];
          //auto length = r.length;
          //double effectiveLength = length - merLen + 1;
          //ofile << name << '\t' << r.length << '\t' << tpms[i] << '\t'
                //<< ((length - estimatedReadLength + 1) > 0 ? kpkms[i] : 0.0) << '\t'
                //<< ((length - merLen + 1) > 0 ? kpkms[i] : 0.0) << '\t'
                //<< kmerCounts[i] << '\t'
                //<< readCounts[i] << '\n';
        //}
        //std::cerr << "retainedCnt = " << retainedCnt << ", nsamps = " << train.n_samples << "\n";

        //ofile.close();
    //return 0;
/*}*/


//int performBiasCorrectionSalmon(
        //bfs::path featureFile,
        //bfs::path expressionFile,
        //bfs::path outputFile,
        //size_t numThreads) {

        //using shark::PCA;

        //auto features = parseFeatureFile(featureFile);
        //std::cerr << "parsed " << features.size() << " features\n";

        //double numMappedReads = 0.0;
        //auto salmonRes = parseSalmonFile(expressionFile, numMappedReads);
        //std::cerr << "parsed " << salmonRes.expressions.size() << " expression values\n";

        //std::vector<size_t> retainedRows;
        //std::vector<double> retainedTPMs;
        //std::vector<std::string> retainedNames;

        //double minLTPM, maxLTPM;
        //minLTPM = std::numeric_limits<double>::max();
        //maxLTPM = -minLTPM;

        //for (auto i : boost::irange(size_t{0}, features.size())) {
                //auto& tname = features[i].name;
                //auto tpm = salmonRes.expressions[tname].tpm; //
                //shark::RealVector v(1);

                //if ( tpm >= 1.0 ) {
                        //retainedRows.emplace_back(i);
                        //retainedNames.push_back(tname);
                        //v(0) = std::log(tpm);
                        //retainedTPMs.push_back(v(0));
                        //minLTPM = std::min(minLTPM, v(0));
                        //maxLTPM = std::max(maxLTPM, v(0));
                //}
        //}


        ////Pca pca;

        //std::vector<shark::RealVector> featMat;
        //std::vector<float> pcavec;
        //size_t fnum = 0;
        //for (auto r : retainedRows) {
                //auto& f = features[r];
                //shark::RealVector v(17);
                //std::vector<double> fv(17);
                //fv[0] = f.gcContent;
                //v(0) = f.gcContent;
                //pcavec.push_back(f.gcContent);
                //for (auto i : boost::irange(size_t{0}, f.diNucleotides.size())) {
                        //v(i+1) = f.diNucleotides[i];
                        //fv[i+1] = f.diNucleotides[i];
                        //pcavec.push_back(f.diNucleotides[i]);
                //}
                //featMat.emplace_back(v);
                //++fnum;
        //}

        //shark::UnlabeledData<shark::RealVector> Xsub = shark::createDataFromRange(featMat);

        //PCA pcao(Xsub, true);
        //auto evals = pcao.eigenvalues();
        //double totalVariance = 0.0;
        //for ( auto e : evals ) { std::cerr << e << "\n"; totalVariance += e; }
        //std::cerr << "totalVariance: " << totalVariance << "\n";
        //double varCutoff = 0.95;
        //double varSum = 0.0;
        //size_t dimCutoff = 0;
        //size_t currentDim = 0;
        //for ( auto e : evals ) {
                //++currentDim;
                //varSum += e;
                //std::cerr << "ev: " << e <<  "\n";
                //if (varSum / totalVariance >= varCutoff) {
                        //dimCutoff = currentDim;
                        //break;
                //}

        //}
        //std::cerr << varCutoff * 100.0 << "% of the variance is explained by " << dimCutoff << " dimensions\n";
        //shark::LinearModel<> enc;
        //pcao.encoder(enc, dimCutoff);
        //auto encodedXsub = enc(Xsub);

        ////shark::UnlabeledData<shark::RealVector> X(18, features.size());
        //Data train;
        //train.set_size(retainedRows.size(), dimCutoff+1);

        ////size_t le = 0;
        //size_t c = 0;
        //for (auto r : retainedRows) {
                //train.X[c][0] = std::log(static_cast<double>(features[r].length));

                //for (auto j : boost::irange(size_t{1}, dimCutoff+1)) {
                        //train.X[c][j] = encodedXsub.element(c)(j-1);
                        ////train.X[c][j] = scores[le];
                        ////train2.X[c][j] = scores[le];
                        ////++le;
                        ////std::cerr << "Train [" << c <<"][" << j        << "] = " << train.X[c][j] << "\n";
                //}
                ////le += (17 - dimCutoff);
                //train.y[c] = retainedTPMs[c];
                //++c;
        //}

        //size_t minDepth = 5;
        //auto reg = std::unique_ptr<RandomForestRegressor>(new RandomForestRegressor(
                //500,
                //train.n_features,
                //5, // max tree depth
                //1, // min_samples_leaf
                //1.0, // features ratio
                //true, // bootstrap
                //true, //out-of-bag
                //true, // compute importance
                //0, // random seed
                //numThreads, // num jobs
                //true // verbose
        //));
        //std::cerr << "there are " << train.n_samples << " samples\n";
        //std::cerr << "there are " << train.n_features << " features\n";
        //reg->build(train.X, train.y, train.n_samples);

        //std::vector<REAL> pred(train.n_samples, 0.0);
        //reg->predict(train.X, &pred[0], train.n_samples, train.n_features);

        //REAL trn_rmse=rmse(&pred[0], train.y, train.n_samples);
        //REAL trn_r2=R2(&pred[0], train.y, train.n_samples);
        //std::cerr << "Train RMSE=" << trn_rmse << ", Correlation Coefficient=" << trn_r2 << "\n";

        //double grandMean = 0.0;
        //size_t ntrain = train.n_samples;
        //for (auto i : boost::irange(size_t{0}, ntrain)) {
                //grandMean += retainedTPMs[i];
        //}
        //grandMean /= train.n_samples;

        //for (auto i : boost::irange(size_t{0}, ntrain)) {
                //pred[i] = grandMean + (retainedTPMs[i] - pred[i]);
        //}

        //trn_rmse=rmse(&pred[0], train.y, train.n_samples);
        //trn_r2=R2(&pred[0], train.y, train.n_samples);
        //std::cerr << "Train RMSE=" << trn_rmse << ", Correlation Coefficient=" << trn_r2 << "\n";

        ////shark::UnlabeledData<shark::RealVector> X()

        //std::ofstream ofile(outputFile.string());
        //for (auto& c : salmonRes.comments) {
                //ofile << c << "\n";
        //}

        //size_t retainedCnt = 0;
        //vector<mpdec> tpms(features.size());
        //for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
          //auto& name = features[i].name;
          //auto& r = salmonRes.expressions[name];
          //if (i == retainedRows[retainedCnt]) {
            //tpms[i] = std::exp(pred[retainedCnt]);
            //++retainedCnt;
          //} else {
              //tpms[i] = r.tpm;
          //}
        //}

        //vector<mpdec> fpkms(features.size());
        //vector<mpdec> estNumReads(features.size());

        //// use TPM estimates to computed estimated read counts
        //mpdec mpzero = 0;
        //mpdec totalNucDenom = 0;
        //for (auto i : boost::irange(size_t{0},  size_t{features.size()})) {
            //double len = features[i].length;
            //totalNucDenom += tpms[i] * len;
        //}

        //for (auto i : boost::irange(size_t{0},  size_t{features.size()})) {
            //mpdec len = features[i].length;
            //estNumReads[i] += ((tpms[i] * len) / totalNucDenom) * numMappedReads;
        //}

        //double oneBillion = 1000000000.0;
        //for (auto i : boost::irange(size_t{0},  size_t{features.size()})) {
            //double len = features[i].length;
            //fpkms[i] = (estNumReads[i] / (len * numMappedReads)) * oneBillion;
        //}

        //for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
          //auto& name = features[i].name;
          //auto& r = salmonRes.expressions[name];
          //auto length = r.length;
          //ofile << name << '\t'
                //<< r.length << '\t'
                //<< tpms[i] << '\t'
                //<< fpkms[i] << '\t'
                //<< estNumReads[i] << '\n';
        //}
        //std::cerr << "retainedCnt = " << retainedCnt << ", nsamps = " << train.n_samples << "\n";

        //ofile.close();
    //return 0;
//}

