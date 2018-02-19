#include "WhiteList.hpp"
#include <cmath>
#include <cassert>
#include <fstream>

namespace alevin {
  namespace whitelist {

    void getCorrelation(std::vector<double>& A,
                        std::vector<double>& B){
      //calculate pearsonr
    }

    using DoubleMatrixT = std::vector<std::vector<double>> ;
    using DoubleVectorT = std::vector<double> ;

    double get_log_likelihood(double prior,
                              DoubleVectorT& sigma,
                              DoubleVectorT& theta,
                              DoubleVectorT& query){
      double logPrior = std::log(prior);
      double mean, variance, likelihood {0.5};

      // Gaussian Probability
      for (size_t i=0; i<sigma.size(); i++){
        variance = sigma[i];
        mean = theta[i];

        likelihood += std::log(2.0 * M_PI * variance);
        likelihood += std::pow((query[i] - mean), 2) / variance;
      }
      likelihood *= -0.5;

      return likelihood+logPrior;
    }

    bool naive_bayes_predict(DoubleMatrixT& sigma,
                             DoubleMatrixT& theta,
                             DoubleVectorT& query,
                             DoubleVectorT& classPrior){
      //std::assert(sigma.size() == theta.size());
      //std::assert(query.size() == theta[0].size());

      double trueProbability = get_log_likelihood(classPrior[0],
                                                  sigma[0],
                                                  theta[0],
                                                  query);
      double falseProbability = get_log_likelihood(classPrior[1],
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
                         size_t numClasses,
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
                                           size_t numLowConfidentBarcode){

      size_t numFalseCells { numLowConfidentBarcode };
      size_t numTrueCells = ( numCells - numFalseCells ) / 2;
      size_t numAmbiguousCells { numTrueCells };
      size_t i, numClasses { 2 };

      std::vector<uint32_t> selectedBarcodes;
      std::vector<uint32_t> classCount (numClasses, 0);
      std::vector<double> classPrior (numClasses, 0.0);

      DoubleMatrixT theta (numClasses, DoubleVectorT (numFeatures, 0.0));
      DoubleMatrixT sigma (numClasses, DoubleVectorT (numFeatures, 0.0));

      naive_bayes_fit(featureCountsMatrix, theta, sigma,
                      classCount, classPrior, numClasses,
                      numTrueCells, numAmbiguousCells, numFalseCells);

      for (i=0; i<numTrueCells; i++){
        selectedBarcodes.emplace_back(i);
      }
      for (; i<numTrueCells*2; i++){
        if (naive_bayes_predict(sigma, theta,
                                featureCountsMatrix[i], classPrior)){
          selectedBarcodes.emplace_back(i);
        }
      }

      return selectedBarcodes;
    }

    template <typename ProtocolT>
    bool performWhitelisting(AlevinOpts<ProtocolT>& aopt,
                             std::vector<uint32_t>& umiCount,
                             std::vector<std::string>& trueBarcodes,
                             CFreqMapT& freqCounter, size_t numGenes,
                             std::vector<std::vector<double>>& countMatrix,
                             spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                             size_t numLowConfidentBarcode){
      // freqCounter has frequency of reads for all detected Barcodes
      // umiCount has frequency of UMI per-cell after knee selection
      // Count matrix file after the deduplicated counts
      // TODO::
      // 1. No mitrochondrial genes
      // 2. No correlation
      // 3. No rRNA information
      // 4. Using all txps i.e. not ignoring txp with 0 values in all the cells

      size_t numCells = trueBarcodes.size();
      size_t numFeatures{4};

      DoubleMatrixT geneCountsMatrix( numCells, DoubleVectorT (numGenes, 0.0));
      DoubleMatrixT featureCountsMatrix( numCells, DoubleVectorT (numFeatures, 0.0));

      // loop over each barcode
      for (size_t i=0; i<trueBarcodes.size(); i++){
        std::vector<double> featureVector(numFeatures);
        std::string currBarcodeName = trueBarcodes[i];
        uint32_t rawBarcodeFrequency;

        // Alignment Rate
        bool indexOk = freqCounter.find(currBarcodeName, rawBarcodeFrequency);
        if ( not indexOk ){
          aopt.jointLog->error("Error: index not find in freq Counter\n"
                               "Please Report the issue on github");
          exit(1);
        }
        featureVector[0] = umiCount[i] / static_cast<double>(rawBarcodeFrequency);

        size_t numNonZeroGeneCount{0};
        double totalCellCount{0.0}, maxCount{0};
        std::vector<double> geneCounts(numGenes, 0.0);
        auto& countVec = countMatrix[i];

        for (size_t j=0; j<countVec.size(); j++){
          auto count = countVec[j];
          numNonZeroGeneCount += 1;
          geneCounts[ txpToGeneMap[j] ] += count;
          totalCellCount += count;
          if (count>maxCount){
            maxCount = count;
          }
        }

        double meanCount = totalCellCount / numNonZeroGeneCount;
        // meanMaxCount
        featureVector[1] = meanCount / maxCount ;
        // dedup Rate
        featureVector[2] = 1.0 - (totalCellCount / umiCount[i]);

        //count of genes over mean
        size_t overMeanCount{0};
        for (auto count: geneCounts){
          if (count > meanCount){
            overMeanCount += 1;
          }
        }
        featureVector[3] = overMeanCount;
        featureCountsMatrix[i] = featureVector;
        geneCountsMatrix[i] = geneCounts;
      }

      std::vector<uint32_t> whitelistBarcodes =
        classifyBarcodes(featureCountsMatrix, numCells,
                         numFeatures, numLowConfidentBarcode);

      std::ofstream whitelistStream;
      auto whitelistFileName = aopt.outputDirectory / "whitelist.txt";
      whitelistStream.open(whitelistFileName.string());
      for (auto i: whitelistBarcodes){
        whitelistStream << trueBarcodes[i] << "\n";
      }
      whitelistStream.close();

      return true;
    }
    template bool performWhitelisting(AlevinOpts<alevin::protocols::DropSeq>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter, size_t numGenes,
                                      std::vector<std::vector<double>>& countMatrix,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::InDrop>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter, size_t numGenes,
                                      std::vector<std::vector<double>>& countMatrix,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Chromium>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter, size_t numGenes,
                                      std::vector<std::vector<double>>& countMatrix,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Custom>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter, size_t numGenes,
                                      std::vector<std::vector<double>>& countMatrix,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      size_t numLowConfidentBarcode);
  }
}
