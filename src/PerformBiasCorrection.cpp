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

        auto numFeatureVectors = features.size();
        auto numSamples = salmonRes.expressions.size();

        bool skipBiasCorrection{false};

        if (numFeatureVectors != numSamples) {
            std::cerr << "The size of the feature map didn't match the "
                      << "number of transcripts.  Bias correction will not "
                      << "be performed\n";
            skipBiasCorrection = true;
        }

        uint32_t minSamples{1000};
        if (numSamples <= minSamples) {
            std::cerr << "There are an insufficient number of transcripts "
                      << "for post-hoc bias correction.  It will not be performed\n";
            skipBiasCorrection = true;
        }

        if (skipBiasCorrection) {
            std::ofstream ofile(outputFile.string());
            for (auto& c : salmonRes.comments) {
                ofile << c << "\n";
            }
            for (auto& kv : salmonRes.expressions) {
                auto& name = kv.first;
                auto& expRecord = kv.second;
                ofile << name << '\t'
                    << expRecord.length << '\t'
                    << expRecord.tpm << '\t'
                    << expRecord.approxCount << '\n';
            }
            ofile.close();
            return 0;
        }


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

        for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
          auto& name = features[i].name;
          auto& r = salmonRes.expressions[name];
          auto length = r.length;
          ofile << name << '\t'
                << r.length << '\t'
                << tpms[i] << '\t'
                << estNumReads[i] << '\n';
        }
        std::cerr << "retainedCnt = " << retainedCnt << ", nsamps = " << train.n_samples << "\n";

        ofile.close();
    return 0;
}


