#include "WhiteList.hpp"

namespace alevin {
  namespace whitelist {
    template <typename ProtocolT>
    bool performWhitelisting(AlevinOpts<ProtocolT>& aopt,
                             std::vector<uint32_t>& umiCount,
                             std::vector<std::string>& trueBarcodes,
                             CFreqMapT& freqCounter,
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
      size_t numFalseCells {1000};
      size_t numTrueCells = ( numCells - numFalseCells ) / 2;
      size_t numAmbiguousCells = numTrueCells;

      spp::sparse_hash_map<int32_t, std::vector<double>> geneCountsMatrix;
      spp::sparse_hash_map<int32_t, std::vector<double>> featureCountsMatrix;
      size_t numFeatures(4);

      // loop over each barcode
      for (size_t i=0; i<trueBarcodes.size(); i++){
        std::vector<double> featureVector(numFeatures);
        std::string currBarcodeName = trueBarcodes[i];
        // Alignment Rate
        featureVector[0] = umiCount[i] / static_cast<double>(freqCounter[currBarcodeName]);

        size_t numNonZeroGeneCount{0};
        double totalCellCount{0.0}, maxCount{0};
        std::vector<double> geneCounts(numGenes);
        auto& countVec = countMatrix[i];

        for (auto j=0; j<countVec.size(); j++){
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
        featureVector[2] = 1.0 - (cellCount / umiCount[i]);

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

      return true;
    }
    template bool performWhitelisting(AlevinOpts<alevin::protocols::DropSeq>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::InDrop>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Chromium>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Custom>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter);
  }
}
