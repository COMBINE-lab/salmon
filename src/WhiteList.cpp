#include "WhiteList.hpp"
#include "tbb/task_scheduler_init.h"
#include <cassert>

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
                      classCount, classPrior, numClasses,
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

    double getPearsonCorrelation(std::vector<double>& first,
                                 std::vector<double>& second){
      if (first.size() != second.size()){
        std::cout << "correlation vectors size does not match" << std::flush;
        exit(1);
      }

      size_t num_elem = first.size();
      double f_sum = std::accumulate(first.begin(), first.end(), 0.0);
      double f_sq_sum = std::accumulate( first.begin(), first.end(), 0.0,
                                         [](const double& l, const double& r){
                                           return (l + r*r);} );
      double s_sum = std::accumulate(second.begin(), second.end(), 0.0);
      double s_sq_sum = std::accumulate( second.begin(), second.end(), 0.0,
                                         [](const double& l, const double& r){
                                           return (l + r*r);} );

      double p_sum = std::inner_product(first.begin(), first.end(),
                                        second.begin(), 0.0);

      double numrtr = p_sum-( (f_sum*s_sum) / num_elem);

      double f_std = f_sq_sum - pow(f_sum, 2)/num_elem;
      double s_std = s_sq_sum - pow(s_sum, 2)/num_elem;

      double denom = pow(f_std * s_std , 0.5);

      if (denom == 0.0 or numrtr == 0.0){
        return 0.0;
      }

      return numrtr / denom;
    }

    uint32_t populate_count_matrix(boost::filesystem::path& outDir,
                               bool inDebugMode,
                               size_t numElem,
                               DoubleMatrixT& countMatrix) {
      boost::iostreams::filtering_istream countMatrixStream;
      countMatrixStream.push(boost::iostreams::gzip_decompressor());
      uint32_t zerod_cells {0};

      auto countMatFilename = outDir / "quants_mat.gz";
      if(not boost::filesystem::exists(countMatFilename)){
        std::cout<<"ERROR: Can't import Binary file quants.mat.gz, it doesn't exist" <<std::flush;
        exit(1);
      }
      countMatrixStream.push(boost::iostreams::file_source(countMatFilename.string(),
                                                           std::ios_base::in | std::ios_base::binary));
      size_t elSize = sizeof(typename std::vector<double>::value_type);
      size_t cellCount {0};
      for (auto& cell: countMatrix){
        cellCount += 1;
        countMatrixStream.read(reinterpret_cast<char*>(cell.data()), elSize * numElem);
        double readCount = std::accumulate(cell.begin(), cell.end(), 0.0);

        if (readCount == 0){
          if (not inDebugMode) {
            std::cout<<"ERROR: Importing counts from binary\n"
                     <<"Cell has 0 reads, should not happen.\n"
                     <<"Saw total "<< cellCount << " Cells before Error"
                     <<std::flush;
            exit(1);
          }
          else {
            zerod_cells += 1;
          }
        }
      }

      return zerod_cells;
    }

    template <typename ProtocolT>
    bool performWhitelisting(AlevinOpts<ProtocolT>& aopt,
                             std::vector<uint32_t>& umiCount,
                             DoubleMatrixT& geneCountsMatrix,
                             std::vector<std::string>& trueBarcodes,
                             CFreqMapT& freqCounter,
                             spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                             size_t numLowConfidentBarcode){
      // freqCounter has frequency of reads for all detected Barcodes
      // umiCount has frequency of UMI per-cell after knee selection
      // Count matrix file after the deduplicated counts
      // TODO::
      // 4. Using all txps i.e. not ignoring txp with 0 values in all the cells
      size_t numCells = trueBarcodes.size();
      size_t numGenes = geneIdxMap.size();
      size_t numFeatures{4};
      if (aopt.useCorrelation){
        numFeatures += 1;
      }

      spp::sparse_hash_set<uint32_t> mRnaGenes, rRnaGenes;
      bool useMitoRna {false}, useRiboRna{false};
      boost::filesystem::path barcodeOrderFile = aopt.outputDirectory / "quants_mat_rows.txt";
      spp::sparse_hash_map<std::string, uint32_t> bcOrderMap;

      if(boost::filesystem::exists(barcodeOrderFile)){
        std::ifstream bcFile(barcodeOrderFile.string());
        std::string bc;
        if(bcFile.is_open()) {
          uint32_t count = 0;
          while(getline(bcFile, bc)) {
            bcOrderMap[bc] = count;
            count += 1;
          }
          bcFile.close();
        }
        aopt.jointLog->info("Done importing order of barcodes \"quants_mat_rows.txt\" file.");
        aopt.jointLog->info("Total {} barcodes found", bcOrderMap.size());
        assert(bcOrderMap.size() == numCells);
      }
      else{
        aopt.jointLog->error("quants_mat_rows.txt file not found. It contains the order of "
                             "the cells counts was dumped.");
        aopt.jointLog->flush();
        exit(1);
      }

      if(boost::filesystem::exists(aopt.mRnaFile)){
        numFeatures += 1;
        useMitoRna = true;
        std::ifstream mRnaFile(aopt.mRnaFile.string());
        std::string gene;
        size_t skippedGenes {0};
        if(mRnaFile.is_open()) {
          while(getline(mRnaFile, gene)) {
            if (geneIdxMap.contains(gene)){
              mRnaGenes.insert(geneIdxMap[ gene ]);
            }
            else{
              skippedGenes += 1;
            }
          }
          mRnaFile.close();
        }
        if (skippedGenes > 0){
          aopt.jointLog->warn("{} mitorna gene(s) does not have transcript in the reference",
                              skippedGenes);
        }
        aopt.jointLog->info("Total {} usable mRna genes", mRnaGenes.size());
      }
      else{
        aopt.jointLog->warn("mrna file not provided; using is 1 less feature for whitelisting");
      }

      if(boost::filesystem::exists(aopt.rRnaFile)){
        numFeatures += 1;
        useRiboRna = true;
        std::ifstream rRnaFile(aopt.rRnaFile.string());
        std::string gene;
        size_t skippedGenes {0};
        if(rRnaFile.is_open()) {
          while(getline(rRnaFile, gene)) {
            if (geneIdxMap.contains(gene)){
              rRnaGenes.insert(geneIdxMap[ gene ]);
            }
            else{
              skippedGenes += 1;
            }
          }
          rRnaFile.close();
        }
        if (skippedGenes > 0){
          aopt.jointLog->warn("{} ribosomal rna gene(s) does not have transcript in the reference",
                              skippedGenes);
        }
        aopt.jointLog->info("Total {} usable rRna genes", rRnaGenes.size());
      }
      else{
        aopt.jointLog->warn("rrna file not provided; using is 1 less feature for whitelisting");
      }

      aopt.jointLog->info("Starting to make feature Matrix");
      DoubleMatrixT featureCountsMatrix( numCells, DoubleVectorT (numFeatures, 0.0));

      // loop over each barcode
      tbb::task_scheduler_init tbbScheduler(aopt.numThreads);
      tbb::parallel_for(
                        BlockedIndexRange(size_t(0), size_t(trueBarcodes.size())),
                        [&featureCountsMatrix, &trueBarcodes,
                         &freqCounter, &aopt, &umiCount, &geneCountsMatrix,
                         useMitoRna, useRiboRna, &bcOrderMap,
                         &mRnaGenes, &rRnaGenes, numGenes](const BlockedIndexRange& range) -> void {
                          for (auto i : boost::irange(range.begin(), range.end())) {
                            uint32_t count_matrix_i = bcOrderMap[ trueBarcodes[i] ];
                            std::vector<double>& featureVector = featureCountsMatrix[i];
                            const std::string& currBarcodeName = trueBarcodes[i];
                            uint32_t rawBarcodeFrequency{0};

                            // Alignment Rate
                            auto indexIt = freqCounter.find(currBarcodeName);
                            bool indexOk = indexIt != freqCounter.end();
                            if ( not indexOk ){
                              aopt.jointLog->error("Error: index {} not found in freq Counter\n"
                                                   "Please Report the issue on github", currBarcodeName,
                                                   freqCounter.size());
                              aopt.jointLog->flush();
                              exit(1);
                            }
                            rawBarcodeFrequency = *indexIt;
                            featureVector[0] = rawBarcodeFrequency ?
                              umiCount[i] / static_cast<double>(rawBarcodeFrequency) : 0.0;

                            double totalCellCount{0.0}, maxCount{0.0};
                            double mitoCount {0.0}, riboCount {0.0};
                            DoubleVectorT& geneCounts = geneCountsMatrix[count_matrix_i];

                            for (size_t j=0; j<geneCounts.size(); j++){
                              auto count = geneCounts[j];
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

                            double meanCount = totalCellCount / numGenes;
                            // meanMaxCount
                            featureVector[1] = maxCount ? meanCount / maxCount : 0.0 ;
                            // dedup Rate
                            featureVector[2] = umiCount[i] ? 1.0 - (totalCellCount / umiCount[i]) : 0.0;

                            //count of genes over mean
                            size_t overMeanCount{0};
                            for (auto count: geneCounts){
                              if (count > meanCount){
                                overMeanCount += 1;
                              }
                            }
                            featureVector[3] = overMeanCount;

                            size_t curFeatIdx = 4;
                            if (useMitoRna) {
                              featureVector[curFeatIdx] = mitoCount ? mitoCount/totalCellCount : 0.0;
                              curFeatIdx += 1;
                            }
                            if (useRiboRna) {
                              featureVector[curFeatIdx] = riboCount ? riboCount/totalCellCount : 0.0;
                            }
                          }
                        });

      aopt.jointLog->info("Done making regular featues");
      size_t numTrueCells = ( numCells - numLowConfidentBarcode ) / 2;

      if (aopt.useCorrelation){
        aopt.jointLog->info("Starting to Make correlation Matrix (Heavy!)");
        tbb::parallel_for(
                          BlockedIndexRange(size_t(0), size_t(numCells)),
                          [&featureCountsMatrix, numTrueCells,
                           &geneCountsMatrix, numFeatures](const BlockedIndexRange& range) -> void {
                            for (auto i : boost::irange(range.begin(), range.end())) {
                              double maxCorr = 0.0;
                              for(size_t j=0; j<numTrueCells; j++){
                                if (i == j){
                                  continue;
                                }
                                double currCorr = getPearsonCorrelation(geneCountsMatrix[i],
                                                                        geneCountsMatrix[j]);
                                if (currCorr > maxCorr){
                                  maxCorr = currCorr;
                                }
                              }
                              featureCountsMatrix[i][numFeatures-1] = maxCorr;
                            }
                          });
      }

      aopt.jointLog->info("Done making feature Matrix");

      if (aopt.dumpfeatures){
        std::ofstream featureStream;
        auto featureFileName = aopt.outputDirectory / "featureDump.txt";
        featureStream.open(featureFileName.string());
        /*for (size_t i=0; i<numTrueCells; i++){
          whitelistStream1 << trueBarcodes[i] << "\n";
          }*/
        /*for (size_t i=0; i<geneCountsMatrix.size(); i++){
          whitelistStream1<< trueBarcodes[i] << "\t";
          for ( auto cell: geneCountsMatrix[i] ){
          whitelistStream1 << cell << "\t";
          }
          whitelistStream1 << "\n";
          }*/
        for(size_t i=0; i<featureCountsMatrix.size(); i++){
          featureStream << trueBarcodes[i];
          for(size_t j=0; j<numFeatures; j++) {
            featureStream << "\t" << featureCountsMatrix[i][j];
          }
          featureStream << "\n";
        }
        featureStream.close();
      }

      std::vector<double> trueProbs;
      std::vector<double> falseProbs;

      std::vector<uint32_t> whitelistBarcodes =
        classifyBarcodes(featureCountsMatrix, numCells,
                         numFeatures, numLowConfidentBarcode,
                         numTrueCells, trueProbs, falseProbs);

      std::ofstream whitelistStream;
      auto whitelistFileName = aopt.outputDirectory / "whitelist.txt";
      whitelistStream.open(whitelistFileName.string());

      std::ofstream predictionStream;
      auto predictionFileName = aopt.outputDirectory / "predictions.txt";
      predictionStream.open(predictionFileName.string());
      predictionStream << "cb\ttrue_prob\tFalse_prob\n";

      for (auto i: whitelistBarcodes){
        whitelistStream << trueBarcodes[i] << "\n";
        if (i >= numTrueCells) {
          predictionStream << trueBarcodes[i] << "\t"
                           << trueProbs[i-numTrueCells] << "\t"
                           << falseProbs[i-numTrueCells] << "\t\n";
        }
      }
      whitelistStream.close();
      predictionStream.close();

      return true;
    }
    template bool performWhitelisting(AlevinOpts<alevin::protocols::DropSeq>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      DoubleMatrixT& geneCountsMatrix,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::InDrop>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      DoubleMatrixT& geneCountsMatrix,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::ChromiumV3>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      DoubleMatrixT& geneCountsMatrix,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Chromium>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      DoubleMatrixT& geneCountsMatrix,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Gemcode>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      DoubleMatrixT& geneCountsMatrix,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::CELSeq>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      DoubleMatrixT& geneCountsMatrix,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Custom>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      DoubleMatrixT& geneCountsMatrix,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      size_t numLowConfidentBarcode);
  }
}
