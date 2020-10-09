#include "WhiteList.hpp"
#include "tbb/task_scheduler_init.h"
#include <cassert>
#include <sstream>

namespace alevin {
  namespace whitelist {

    double get_log_likelihood(double prior,
                              DoubleVectorT& sigma,
                              DoubleVectorT& theta,
                              DoubleVectorT& query){
      double logPrior = prior ? std::log(prior) : 0.0;
      double mean, variance, likelihood {0.5};

      // Gaussian Probability
      for (size_t i=0; i<sigma.size(); i++){
        variance = sigma[i];
        mean = theta[i];

        likelihood += variance ? std::log(2.0 * M_PI * variance) : 0.0;
        likelihood += variance ? std::pow((query[i] - mean), 2) / variance : 0.0;
      }
      likelihood *= -0.5;
      return likelihood+logPrior;
    }

    bool naive_bayes_predict(DoubleMatrixT& sigma,
                             DoubleMatrixT& theta,
                             DoubleVectorT& query,
                             DoubleVectorT& classPrior,
                             double& trueProbability,
                             double& falseProbability){
      trueProbability = get_log_likelihood(classPrior[0],
                                           sigma[0],
                                           theta[0],
                                           query);
      falseProbability = get_log_likelihood(classPrior[1],
                                            sigma[1],
                                            theta[1],
                                            query);

      if ( trueProbability < falseProbability ){
        return false;
      }
      else{
        return true;
      }
    }

    void update_mean_variance(DoubleMatrixT features,
                              DoubleMatrixT& theta,
                              DoubleMatrixT& sigma,
                              std::vector<uint32_t>& classCount,
                              uint32_t is_true){
      size_t numFeatures = features[0].size();

      for ( auto& feature: features ){
        for (size_t i=0; i<numFeatures; i++){
          theta[is_true][i] += feature[i];
        }
      }

      // get mean of the dataset
      for ( size_t i=0; i<theta[is_true].size(); i++ ){
        theta[is_true][i] /= classCount[is_true];
      }

      // get variance
      for ( auto& feature: features ){
        for (size_t i=0; i<numFeatures; i++){
          sigma[is_true][i] += std::pow((feature[i] - theta[is_true][i]), 2);
        }
      }

      // get variance of the dataset
      for ( size_t i=0; i<sigma[is_true].size(); i++ ){
        sigma[is_true][i] /= classCount[is_true];
      }
    }

    void naive_bayes_fit(DoubleMatrixT& features,
                         DoubleMatrixT& theta,
                         DoubleMatrixT& sigma,
                         std::vector<uint32_t>& classCount,
                         DoubleVectorT& classPrior,
                         /*size_t numClasses,*/
                         size_t numTrueCells,
                         size_t numAmbiguousCells,
                         size_t numFalseCells){
      //std::assert( theta.size() == sigma.size() );
      size_t numCells = features.size();
      size_t numFeatures = features[0].size();

      size_t trueStartIdx{0}, trueEndIdx{numTrueCells-1};
      size_t falseStartIdx{numTrueCells+numAmbiguousCells};
      size_t falseEndIdx{numCells-1};

      //std::assert( falseEndIdx-falseStartIdx == numFalseCells);

      classCount[0] = numTrueCells;
      classCount[1] = numFalseCells;

      update_mean_variance(DoubleMatrixT (features.begin()+trueStartIdx, features.begin()+trueEndIdx),
                           theta, sigma, classCount, 0);
      update_mean_variance(DoubleMatrixT (features.begin()+falseStartIdx, features.begin()+falseEndIdx),
                           theta, sigma, classCount, 1);

      // calculate prior
      size_t norm = classCount[0] + classCount[1];
      classPrior[0] = classCount[0] / norm;
      classPrior[1] = classCount[1] / norm;
    }

    // Implementation from https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/naive_bayes.py
    std::vector<uint32_t> classifyBarcodes(DoubleMatrixT& featureCountsMatrix,
                                           size_t numCells, size_t numFeatures,
                                           size_t numLowConfidentBarcode,
                                           size_t numTrueCells,
                                           std::vector<double>& trueProb,
                                           std::vector<double>& falseProb){

      size_t numFalseCells { numLowConfidentBarcode };
      size_t numAmbiguousCells { numCells - numTrueCells - numLowConfidentBarcode };
      size_t i, numClasses { 2 };

      std::vector<uint32_t> selectedBarcodes;
      std::vector<uint32_t> classCount (numClasses, 0);
      std::vector<double> classPrior (numClasses, 0.0);

      DoubleMatrixT theta (numClasses, DoubleVectorT (numFeatures, 0.0));
      DoubleMatrixT sigma (numClasses, DoubleVectorT (numFeatures, 0.0));

      naive_bayes_fit(featureCountsMatrix, theta, sigma,
                      classCount, classPrior, /*numClasses,*/
                      numTrueCells, numAmbiguousCells, numFalseCells);

      //trueProb.resize(numAmbiguousCells, 0.0);
      //falseProb.resize(numAmbiguousCells, 0.0);

      for (auto vec: sigma){
        for (auto cell : vec){
          std::cout<<cell<<"\t";
        }
        std::cout<<"\n";
      }

      for (i=0; i<numTrueCells; i++){
        selectedBarcodes.emplace_back(i);
      }
      auto ambiguousCellOffset = numTrueCells + numAmbiguousCells;
      for (; i<ambiguousCellOffset; i++){
        double trPb{0.0}, flPb{0.0};
        if (naive_bayes_predict(sigma, theta,
                                featureCountsMatrix[i],
                                classPrior, trPb, flPb)){
          selectedBarcodes.emplace_back(i);

        }
        trueProb.emplace_back( trPb );
        falseProb.emplace_back( flPb );
      }

      return selectedBarcodes;
    }

    template <typename ProtocolT>
    bool performWhitelisting(AlevinOpts<ProtocolT>& aopt,
                             std::vector<std::string>& trueBarcodes,
                             bool useRibo, bool useMito,
                             size_t numLowConfidentBarcode){
      size_t numCells = trueBarcodes.size();
      size_t numFeatures {5};
      if (useMito) { ++numFeatures; }
      if (useRibo) { ++numFeatures; }
      aopt.numFeatures = numFeatures;

      aopt.jointLog->info("Starting to make feature Matrix");
      size_t numTrueCells = ( numCells - numLowConfidentBarcode ) / 2;
      DoubleMatrixT featureCountsMatrix( numCells, DoubleVectorT (numFeatures, 0.0));

      {
        // populating features
        boost::filesystem::path featuresFileName = aopt.outputDirectory / "featureDump.txt";
        std::ifstream featuresfile(featuresFileName.string());

        std::string line, token;
        size_t cellId {0};
        while( std::getline( featuresfile, line) ) {
          cellId += 1;
          if (cellId == 1) { continue; }

          size_t featureId {0};
          std::vector<double> features;
          std::stringstream buffer(line);
          while( getline( buffer, token, '\t') ) {
            featureId += 1;
            if (aopt.dumpUmiGraph and featureId > 4+numFeatures) { break; }
            if (featureId > 4) { features.emplace_back( std::strtod(token.c_str(),
                                                                    nullptr)); }
          }

          if (features.size() != numFeatures) {
            aopt.jointLog->error("Size of the feature file doesn't match.\n"
                                 "Please report the error on github.");
            aopt.jointLog->flush();
            exit(84);
          }

          featureCountsMatrix[cellId-2] = features;
        }

        if (cellId-1 != numCells) {
          aopt.jointLog->error("The features file has less number of cells.\n"
                               "Please report the error on github.");
          aopt.jointLog->flush();
          exit(84);
        }
        aopt.jointLog->info("Done making feature Matrix");
        aopt.jointLog->flush();
      } // done populating features

      std::vector<double> trueProbs, falseProbs;
      std::vector<uint32_t> whitelistBarcodes =
        classifyBarcodes(featureCountsMatrix, numCells,
                         numFeatures, numLowConfidentBarcode,
                         numTrueCells, trueProbs, falseProbs);

      auto whitelistFileName = aopt.outputDirectory / "whitelist.txt";
      std::ofstream whitelistStream(whitelistFileName.string());
      for (auto i: whitelistBarcodes){
        whitelistStream << trueBarcodes[i] << std::endl;
      }
      whitelistStream.close();
      aopt.intelligentCutoff = whitelistBarcodes.size();

      if (aopt.dumpfeatures) {
        std::ofstream predictionStream;
        auto predictionFileName = aopt.outputDirectory / "predictions.txt";
        predictionStream.open(predictionFileName.string());
        predictionStream << "cb\ttrue_prob\tFalse_prob\n";

        for (auto i: whitelistBarcodes){
          if (i >= numTrueCells) {
            predictionStream << trueBarcodes[i] << "\t"
                             << trueProbs[i-numTrueCells] << "\t"
                             << falseProbs[i-numTrueCells] << "\t\n";
          }
        }
        predictionStream.close();
      }//end-if dumpfeatures

      return true;
    }
    template bool performWhitelisting(AlevinOpts<alevin::protocols::DropSeq>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      bool useRibo, bool useMito,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::CITESeq>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      bool useRibo, bool useMito,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::InDrop>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      bool useRibo, bool useMito,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::ChromiumV3>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      bool useRibo, bool useMito,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Chromium>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      bool useRibo, bool useMito,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Gemcode>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      bool useRibo, bool useMito,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::CELSeq>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      bool useRibo, bool useMito,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::CELSeq2>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      bool useRibo, bool useMito,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::QuartzSeq2>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      bool useRibo, bool useMito,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Custom>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      bool useRibo, bool useMito,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::CustomGeometry>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      bool useRibo, bool useMito,
                                      size_t numLowConfidentBarcode);
  }
}
