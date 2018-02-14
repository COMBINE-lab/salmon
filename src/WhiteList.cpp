#include "WhiteList.hpp"
#include <cmath>

namespace alevin {
  namespace whitelist {

    void getCorrelation(std::vector<double>& A,
                        std::vector<double>& B){
      //calculate pearsonr
    }

    use DoubleMatrixT std::vector<std::vector<double>> ;
    use DoubleVectorT std::vector<double> ;

    double get_log_likelihood(double prior,
                              DoubleVectorT& sigma,
                              DoubleVectorT& theta,
                              DoubleVectorT& query){
      double logPrior = std::log(prior);
      double mean, variance, likelihood {0.5};

      // Gaussian Probability
      for (size_t i=0; i<sigma.len(); i++){
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
      std::assert(sigma.size() == theta.size());
      std::assert(query.size() == theta[0].size());

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

    void update_mean_variance(DoubleMatrixT& features,
                              DoubleMatrixT& theta,
                              DoubleMatrixT& sigma,
                              DoubleVectorT& classCount,
                              uint8_t is_true){
      std::assert( theta.size() == sigma.size() );
      size_t numFeatures = features[0].size();

      for ( auto& feature: features ){
        for (size_t j=0; j<numFeatures; j++){
          theta[is_true][i] += feature[i];
        }
      }
    }

    // Implementation from https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/naive_bayes.py
    void naive_bayes_fit(DoubleMatrixT& positiveFeatures,
                         DoubleMatrixT& negativeFeatures){
      size_t numFeatures { positiveFeatures[0].size() };
      size_t numClasses {2};

      DoubleMatrixT theta (numClasses, DoubleVectorT (numFeatures, 0.0));
      DoubleMatrixT sigma (numClasses, DoubleVectorT (numFeatures, 0.0));

      std::vector<uint32_t> classCount (numClasses, 0);
      std::vector<double> classPrior (numClasses, 0.0);

      update_mean_variance(positiveFeatures, theta, sigma, classCount, 1);
      update_mean_variance(negativeFeatures, theta, sigma, classCount, 0);

      double var_smoothing {1e-9};
      // get mean across full dataset (i.e. both +ve and -ve)
      double maxIntraFeatureVariance = std::max(mean);
      double epsilon = var_smoothing * maxIntraFeatureVariance;
    }

    std::vector<uint32_t> classifyBarcodes(spp::sparse_hash_map<int32_t, std::vector<double>>& featureCountsMatrix,
                                           size_t numCells){
      size_t numFalseCells {1000};
      size_t numTrueCells = ( numCells - numFalseCells ) / 2;
      size_t numAmbiguousCells {numTrueCells};
      std::vector<uint32_t> selectedBarcodes;

      return selectedBarcodes;
    }

    template <typename ProtocolT>
    bool performWhitelisting(AlevinOpts<ProtocolT>& aopt,
                             std::vector<uint32_t>& umiCount,
                             std::vector<std::string>& trueBarcodes,
                             CFreqMapT& freqCounter, size_t numGenes,
                             std::vector<std::vector<double>>& countMatrix,
                             spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap){
      // freqCounter has frequency of reads for all detected Barcodes
      // umiCount has frequency of UMI per-cell after knee selection
      // Count matrix file after the deduplicated counts
      // TODO::
      // 1. No mitrochondrial genes
      // 2. No correlation
      // 3. No rRNA information
      // 4. Using all txps i.e. not ignoring txp with 0 values in all the cells

      size_t numCells = trueBarcodes.size();

      spp::sparse_hash_map<int32_t, std::vector<double>> geneCountsMatrix;
      spp::sparse_hash_map<int32_t, std::vector<double>> featureCountsMatrix;
      size_t numFeatures(4);

      // loop over each barcode
      for (size_t i=0; i<trueBarcodes.size(); i++){
        std::vector<double> featureVector(numFeatures);
        std::string currBarcodeName = trueBarcodes[i];
        uint32_t rawBarcodeFrequency;

        // Alignment Rate
        bool indexOk = freqCounter.find(currBarcodeName, rawBarcodeFrequency);
        featureVector[0] = umiCount[i] / static_cast<double>(rawBarcodeFrequency);

        size_t numNonZeroGeneCount{0};
        double totalCellCount{0.0}, maxCount{0};
        std::vector<double> geneCounts(numGenes);
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
        classifyBarcodes(featureCountsMatrix, numCells);

      return true;
    }
    template bool performWhitelisting(AlevinOpts<alevin::protocols::DropSeq>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter, size_t numGenes,
                                      std::vector<std::vector<double>>& countMatrix,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::InDrop>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter, size_t numGenes,
                                      std::vector<std::vector<double>>& countMatrix,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Chromium>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter, size_t numGenes,
                                      std::vector<std::vector<double>>& countMatrix,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Custom>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter, size_t numGenes,
                                      std::vector<std::vector<double>>& countMatrix,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap);
  }
}
