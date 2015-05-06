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

#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

#include "shark/Data/Dataset.h"
#include "shark/Algorithms/Trainers/PCA.h"
#include "shark/Algorithms/Trainers/RFTrainer.h"
#include "shark/Algorithms/Trainers/CARTTrainer.h"

#include "tensemble/TypeDef.h"
#include "tensemble/RandomForestRegressor.h"
#include "tensemble/RandomForestClassifier.h"
#include "tensemble/GBMRegressor.h"
#include "tensemble/GBMClassifier.h"
#include "tensemble/ReadData.h"

#define DEFAULT_N_TREES 100
#define DEFAULT_N_JOBS 1
#define DEFAULT_MAX_FEATURES_RATIO 1.0
#define DEFAULT_MIN_SAMPLE_LEAF 10
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
using Kmer = uint64_t;

struct TranscriptFeatures{
	std::string name;
	size_t length;
	double gcContent;
	std::array<Kmer, 16> diNucleotides;
};

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
};

struct ExpressionResults {
	std::vector<std::string> comments;
	std::unordered_map<std::string, TranscriptResult> expressions;
};

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
			res.expressions[tname] = tr;
			// eat the newline
			char nline; ifile.get(nline);
		}

		if (ifile.peek() == EOF) { break; }
	}

	return res;
}

int main(int argc, char* argv[]) {

	using shark::PCA;
	//PCA pca(data)
	bfs::path featureFile = bfs::path(argv[1]);
	auto features = parseFeatureFile(featureFile);
	std::cerr << "parsed " << features.size() << " features\n";

	bfs::path expressionFile = bfs::path(argv[2]);
	auto sfres = parseSailfishFile(expressionFile);
	std::cerr << "parsed " << sfres.expressions.size() << " expression values\n";

	std::vector<shark::RealVector> RPKMLabels;//(features.size(), 1);
	std::vector<size_t> retainedRows;
	std::vector<double> retainedRPKMs;
	std::vector<std::string> retainedNames;

	double minLRPKM, maxLRPKM;
	minLRPKM = std::numeric_limits<double>::max();
	maxLRPKM = -minLRPKM;
	double avgLRPKM = 0.0;

	for (auto i : boost::irange(size_t{0}, features.size())) {
		auto& tname = features[i].name;
		auto rpkm = sfres.expressions[tname].rpkm;
		shark::RealVector v(1);

		if ( rpkm >= 0.001 ) {
			retainedRows.emplace_back(i);
			retainedNames.push_back(tname);
			v(0) = std::log(rpkm);
			retainedRPKMs.push_back(v(0));

			if (std::fabs(rpkm - std::exp(v(0))) > 1e-2) {
				std::cerr << "EXP DOES NOT INVERT LOG!!\n";
			}
			minLRPKM = std::min(minLRPKM, v(0));
			maxLRPKM = std::max(maxLRPKM, v(0));
			avgLRPKM += v(0);
			RPKMLabels.emplace_back(v);
		}
	}
	avgLRPKM /= retainedRows.size();


	std::vector<shark::RealVector> featMat;
	size_t fnum = 0;
	for (auto r : retainedRows) {
		auto& f = features[r];
		shark::RealVector v(17);
		v(0) = f.gcContent;
		for (auto i : boost::irange(size_t{0}, f.diNucleotides.size())) { v(i+1) = f.diNucleotides[i]; }
		featMat.emplace_back(v);
		++fnum;
	}

	shark::UnlabeledData<shark::RealVector> Xsub = shark::createDataFromRange(featMat);
	PCA pca(Xsub, true);
	auto evals = pca.eigenvalues();
	double totalVariance = 0.0;
	for ( auto e : evals ) { totalVariance += e; }
	double varCutoff = 0.95;
	double varSum = 0.0;
	size_t dimCutoff = 0;
	size_t currentDim = 0;
	for ( auto e : evals ) {
		++currentDim;
		varSum += e;
		if (varSum / totalVariance >= varCutoff) {
			dimCutoff = currentDim;
			break;
		}
		//std::cerr << "ev: " << e <<  "\n";
	}
	std::cerr << varCutoff * 100.0 << "% of the variance is explained by " << dimCutoff << " dimensions\n";

	shark::LinearModel<> enc;
	pca.encoder(enc, dimCutoff);
	auto encodedXsub = enc(Xsub);

	std::vector<shark::RealVector> X;
	std::vector<shark::RealVector> labels;

	Data train;
	Data train2;
	train.set_size(retainedRows.size(), dimCutoff+1);
	train2.set_size(retainedRows.size(), dimCutoff+1);

	size_t c = 0;
	for (auto r : retainedRows) {
		shark::RealVector v(dimCutoff + 1);
		shark::RealVector l(1);

		v(0) = std::log(static_cast<double>(features[r].length));
		train.X[c][0] = std::log(static_cast<double>(features[r].length));
		train2.X[c][0] = std::log(static_cast<double>(features[r].length));
		for (auto j : boost::irange(size_t{1}, dimCutoff+1)) {
			train.X[c][j] = encodedXsub.element(c)(j-1);
			train2.X[c][j] = encodedXsub.element(c)(j-1);
			v(j) = encodedXsub.element(c)(j-1);
			//std::cerr << "Train [" << c <<"][" << j	 << "] = " << train.X[c][j] << "\n";
		}
		train.y[c] = retainedRPKMs[c];// std::log(sfres[features[r].name].expressions.rpkm);
		l(0) = retainedRPKMs[c];

		X.emplace_back(v);
		labels.emplace_back(l);

		++c;
	}

	size_t numRetainedSamples = retainedRows.size();

// #define DEFAULT_N_TREES 100
// #define DEFAULT_N_JOBS 1
// #define DEFAULT_MAX_FEATURES_RATIO 1.0
// #define DEFAULT_MIN_SAMPLE_LEAF 2
// #define DEFAULT_MAX_DEPTH   4
// #define DEFAULT_SUBSAMPLE   1.0
// #define DEFAULT_SPLIT_CRITERION CRITERION_MSE
// #define DEFAULT_LOSS SQUARE_LOSS
// #define DEFAULT_LEARN_RATE 0.1
// #define DEFAULT_OOB 1
// #define DEFAULT_VERBOSE 0
// #define DEFAULT_BOOTSTRAP 1
// #define DEFAULT_COMPUTE_IMPORTANCE 0

	/** Random Forest Regression **/
	/*
	auto reg = std::unique_ptr<RandomForestRegressor>(new RandomForestRegressor(
		1000,
		train.n_features,
		DEFAULT_MAX_DEPTH,
		DEFAULT_MIN_SAMPLE_LEAF,
		DEFAULT_MAX_FEATURES_RATIO,
		true, // bootstrap
		true, //out-of-bag
		false,
		0,
		16,
		true
	));
	*/

	/*
	auto reg = std::unique_ptr<GBMRegressor>(new GBMRegressor(
		SQUARE_LOSS,
		500,
		train.n_features,
		DEFAULT_MAX_DEPTH,
		DEFAULT_MIN_SAMPLE_LEAF,
		DEFAULT_MAX_FEATURES_RATIO,
		DEFAULT_SUBSAMPLE, // subsample
		0.25, //learn rate
		true, //out-of-bag
		false, // compute imporance
		34239, // random seed
		30, // num jobs
		true // verbose
	));*/

	/*
	std::cerr << "there are " << numRetainedSamples << " samples\n";
	std::cerr << "there are " << train.n_features << " features\n";
	reg->build(train.X, train.y, numRetainedSamples);

	REAL* pred = new REAL[numRetainedSamples];
	for (auto i : boost::irange(size_t{0}, size_t{numRetainedSamples})) { pred[i] = 0.0; }
	reg->predict(train2.X, pred, numRetainedSamples, train.n_features);

 REAL trn_rmse=rmse(pred, train.y, numRetainedSamples);
  REAL trn_r2=R2(pred, train.y, numRetainedSamples);
  std::cerr << "Train RMSE=" << trn_rmse << ", Correlation Coefficient=" << trn_r2 << "\n";
	*/


	shark::Data<shark::RealVector> input = shark::createDataFromRange(X);
	shark::Data<shark::RealVector> regLabels = shark::createDataFromRange(labels);

	shark::LabeledData<shark::RealVector, shark::RealVector> regressionData(input, regLabels);

	shark::RFTrainer trainer;
	shark::RFClassifier model;

	trainer.setNodeSize(20);
	trainer.setNTrees(500);
	trainer.setOOBratio(0.001);

	//std::cerr << "number of classes = " << shark::numberOfClasses(regressionSubset) << "\n";
	//std::cerr << "number of classes = " << shark::numberOfClasses(regressionSubset) << "\n";
	trainer.train(model, regressionData);

	auto predictionData = model(regressionData.inputs());
	std::cerr << "inputs\n";
	std::vector<double> pred(retainedRows.size(), 0.0);
	//std::cerr << regressionSubset << "\n";
	std::cerr << predictionData << "\n";
	size_t ctr = 0;
	std::cerr << "NUM ELEMENTS: " << predictionData.numberOfElements() << "\n";
	for (auto i : boost::irange(size_t{0}, predictionData.numberOfElements())) {
		pred[i] = predictionData.element(i)(0);
	}
	std::cerr << "cha!\n";

	//auto& inputLabels = regressionSubset.labels();
	/*
	auto mm = std::minmax_element(train.y, train.y + numRetainedSamples);
	double minLRPKM = std::get<0>(mm);
	double maxLRPKM = std::get<1>(mm);
	*/

	for (auto i : boost::irange(size_t{0}, size_t{numRetainedSamples})) {
		pred[i] = RPKMLabels[i](0) - pred[i];
	}

	auto mmpred = std::minmax_element(pred.begin(), pred.begin() + numRetainedSamples);
	double minPred = *(mmpred.first);
	double maxPred = *(mmpred.second);

	double scale = std::fabs(maxLRPKM - minLRPKM) / std::fabs(maxPred - minPred);

	std::cerr << "min,max LRPKM : " << minLRPKM << ", " << maxLRPKM << "\n";
	std::cerr << "min, max pred : " << minPred << ", " << maxPred  << "\n";
	std::cerr << "SCALE: " << scale << "\n";


	minPred = std::numeric_limits<double>::max();
	for (auto i : boost::irange(size_t{0}, size_t{numRetainedSamples})) {
		pred[i] *= scale;
		minPred = std::min(minPred, pred[i]);
	}

	double shift{minLRPKM - minPred};
	minPred = std::numeric_limits<double>::max();
	maxPred = -minPred;
	for (auto i : boost::irange(size_t{0}, size_t{numRetainedSamples})) {
		pred[i] += shift;
		minPred = std::min(minPred, pred[i]);
		maxPred = std::max(maxPred, pred[i]);
	}

	/*
   trn_rmse=rmse(pred, train.y, numRetainedSamples);
   trn_r2=R2(pred, train.y, numRetainedSamples);
  std::cerr << "Train RMSE=" << trn_rmse << ", Correlation Coefficient=" << trn_r2 << "\n";
	*/

	std::ofstream ofile(argv[3]);
	for (auto& c : sfres.comments) {
		ofile << c << "\n";
	}

	size_t retainedCnt = 0;
	for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {

		double rpkm = 0.0;
		auto& name = features[i].name;
		auto& r = sfres.expressions[name];

		if (i == retainedRows[retainedCnt]) {
			//rpkm = std::exp(pred[retainedCnt]);
			rpkm = std::exp(pred[retainedCnt]);
			++retainedCnt;
		} else {
			rpkm = r.rpkm;
		}
		/*
		if (i == retainedRows[retainedCnt]) {
			if (retainedNames[retainedCnt] != features[i].name) {
				std::cerr << "AHHH!!!\n";
			}

			rpkm = std::exp(RPKMLabels[retainedCnt](0));// pred[retainedCnt]);

			++retainedCnt;
		} else {
			rpkm = sfres.expressions[features[i].name].rpkm;
		}
		*/
		//std::cerr << "Feature = " << features[i].name << ", rpkm = " << rpkm << "\n";

		ofile << name << '\t' << r.length << '\t' << r.tpm << '\t' << rpkm << '\n';
	}
	std::cerr << "retainedCnt = " << retainedCnt << ", nsamps = " << numRetainedSamples << "\n";
	//for ( auto i : retainedRows ) { std::cerr << i << " "; }
	ofile.close();

}
