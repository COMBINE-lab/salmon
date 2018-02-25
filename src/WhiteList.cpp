#include "WhiteList.hpp"
#include <cmath>
#include <cassert>
#include <fstream>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

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
                             DoubleVectorT& classPrior){
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
                                featureCountsMatrix[i],
                                classPrior)){
          selectedBarcodes.emplace_back(i);
        }
      }

      return selectedBarcodes;
    }

    void populate_count_matrix(boost::filesystem::path& outDir,
                               DoubleMatrixT& countMatrix){
      boost::iostreams::filtering_istream countMatrixStream;
      countMatrixStream.push(boost::iostreams::gzip_decompressor());

      auto countMatFilename = outDir / "quants_mat.gz";
      countMatrixStream.push(boost::iostreams::file_source(countMatFilename.string(),
                                                           std::ios_base::in | std::ios_base::binary));
      size_t numTxps = countMatrix[0].size();
      size_t elSize = sizeof(typename std::vector<double>::value_type);
      for (auto& cell: countMatrix){
        countMatrixStream.read(reinterpret_cast<char*>(cell.data()), elSize * numTxps);
      }
    }

    template <typename ProtocolT>
    bool performWhitelisting(AlevinOpts<ProtocolT>& aopt,
                             std::vector<uint32_t>& umiCount,
                             std::vector<std::string>& trueBarcodes,
                             CFreqMapT& freqCounter,
                             spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                             spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                             size_t numLowConfidentBarcode){
      // freqCounter has frequency of reads for all detected Barcodes
      // umiCount has frequency of UMI per-cell after knee selection
      // Count matrix file after the deduplicated counts
      // TODO::
      // 2. No correlation
      // 4. Using all txps i.e. not ignoring txp with 0 values in all the cells

      size_t numCells = trueBarcodes.size();
      size_t numGenes = geneIdxMap.size();
      size_t numTxps = txpToGeneMap.size();
      size_t numFeatures{4};

      DoubleMatrixT countMatrix(numCells, DoubleVectorT (numTxps, 0.0));
      //DoubleMatrixT geneCountsMatrix( numCells, DoubleVectorT (numGenes, 0.0));
      populate_count_matrix(aopt.outputDirectory, countMatrix);

      std::ofstream qFile;
      boost::filesystem::path qFilePath = aopt.outputDirectory / "txt_count_mat_new.txt";
      qFile.open(qFilePath.string());
      for (auto& row : countMatrix) {
        for (auto cell : row) {
          qFile << cell << ',';
        }
        qFile << "\n";
      }

      spp::sparse_hash_set<uint32_t> mRnaGenes, rRnaGenes;
      bool useMitoRna {false}, useRiboRna{false};

      if(boost::filesystem::exists(aopt.mRnaFile)){
        numFeatures += 1;
        useMitoRna = true;
        std::ifstream mRnaFile(aopt.mRnaFile.string());
        std::string gene;
        if(mRnaFile.is_open()) {
          while(getline(mRnaFile, gene)) {
            if (geneIdxMap.contains(gene)){
              mRnaGenes.insert(geneIdxMap[ gene ]);
            }
            else{
              aopt.jointLog->error("{} mrna gene not found in txp tp gene map", gene);
              aopt.jointLog->flush();
              exit(1);
            }
          }
          mRnaFile.close();
        }
        aopt.jointLog->info("Done importing mRnaFile");
        aopt.jointLog->info("Total {} mRna genes", mRnaGenes.size());
      }
      else{
        aopt.jointLog->warn("mrna file not provided; using is 1 less feature for whitelisting");
      }

      if(boost::filesystem::exists(aopt.rRnaFile)){
        numFeatures += 1;
        useRiboRna = true;
        std::ifstream rRnaFile(aopt.rRnaFile.string());
        std::string gene;
        if(rRnaFile.is_open()) {
          while(getline(rRnaFile, gene)) {
            if (geneIdxMap.contains(gene)){
              rRnaGenes.insert(geneIdxMap[ gene ]);
            }
            else{
              aopt.jointLog->error("{} rrna gene not found in txp tp gene map", gene);
              aopt.jointLog->flush();
              exit(1);
            }
          }
          rRnaFile.close();
        }
        aopt.jointLog->info("Done importing rRnaFile");
        aopt.jointLog->info("Total {} rRna genes", rRnaGenes.size());
      }
      else{
        aopt.jointLog->warn("rrna file not provided; using is 1 less feature for whitelisting");
      }

      aopt.jointLog->info("Starting to make feature Matrix");
      DoubleMatrixT featureCountsMatrix( numCells, DoubleVectorT (numFeatures, 0.0));

      // loop over each barcode
      for (size_t i=0; i<trueBarcodes.size(); i++){
        std::vector<double> featureVector(numFeatures);
        std::string currBarcodeName = trueBarcodes[i];
        uint32_t rawBarcodeFrequency{0};

        // Alignment Rate
        bool indexOk = freqCounter.find(currBarcodeName, rawBarcodeFrequency);
        if ( not indexOk ){
          aopt.jointLog->error("Error: index {} not found in freq Counter\n"
                               "Please Report the issue on github", currBarcodeName,
                               freqCounter.size());
          aopt.jointLog->flush();
          exit(1);
        }
        featureVector[0] = umiCount[i] / static_cast<double>(rawBarcodeFrequency);

        size_t numNonZeroGeneCount{0};
        double totalCellCount{0.0}, maxCount{0};
        std::vector<double> geneCounts(numGenes, 0.0);
        auto& countVec = countMatrix[i];
        double mitoCount {0.0}, riboCount {0.0};

        for (size_t j=0; j<countVec.size(); j++){
          auto count = countVec[j];
          numNonZeroGeneCount += 1;
          geneCounts[ txpToGeneMap[j] ] += count;
          totalCellCount += count;
          if (count>maxCount){
            maxCount = count;
          }
          if (useMitoRna and mRnaGenes.contains(j)){
            mitoCount += count;
          }
          if (useRiboRna and rRnaGenes.contains(j)){
            riboCount += count;
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

        if (useMitoRna) {
          featureVector[4] = mitoCount ? mitoCount/totalCellCount : 0.0;
        }
        if (useRiboRna) {
          featureVector[5] = riboCount ? riboCount/totalCellCount : 0.0;
        }

        //geneCountsMatrix[i] = geneCounts;
        featureCountsMatrix[i] = featureVector;
      }
      aopt.jointLog->info("Done making feature Matrix");

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
                                      CFreqMapT& freqCounter,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::InDrop>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Chromium>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Custom>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      size_t numLowConfidentBarcode);
  }
}
