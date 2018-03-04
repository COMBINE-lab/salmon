#include "WhiteList.hpp"
#include "tbb/task_scheduler_init.h"

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
                                           size_t numLowConfidentBarcode,
                                           size_t numTrueCells){

      size_t numFalseCells { numLowConfidentBarcode };
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

      //for (auto vec: sigma){
      //  for (auto cell : vec){
      //    std::cout<<cell<<"\t";
      //  }
      //  std::cout<<"\n";
      //}

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

      double f_mean = f_sum / num_elem;
      double s_mean = s_sum / num_elem;

      double f_std = pow( (f_sq_sum / num_elem) - pow(f_mean, 2), 0.5);
      double s_std = pow( (s_sq_sum / num_elem) - pow(s_mean, 2), 0.5);

      double denom = f_std * s_std * num_elem;
      if (denom == 0){
        return 0;
      }
      double numer = 0.0;
      std::for_each(first.begin(), first.end(), [&](double& d) { numer += (d-f_mean) ;});
      std::for_each(second.begin(), second.end(), [&](double& d) { numer += (d-f_mean) ;});

      return numer / denom ;
    }

    void populate_count_matrix(boost::filesystem::path& outDir,
                               DoubleMatrixT& countMatrix){
      boost::iostreams::filtering_istream countMatrixStream;
      countMatrixStream.push(boost::iostreams::gzip_decompressor());

      auto countMatFilename = outDir / "quants_mat.gz";
      if(not boost::filesystem::exists(countMatFilename)){
        std::cout<<"ERROR: Can't import Binary file quants.mat.gz, it doesn't exist" <<std::flush;
        exit(1);
      }
      countMatrixStream.push(boost::iostreams::file_source(countMatFilename.string(),
                                                           std::ios_base::in | std::ios_base::binary));
      size_t numTxps = countMatrix[0].size();
      size_t elSize = sizeof(typename std::vector<double>::value_type);
      size_t cellCount {0};
      for (auto& cell: countMatrix){
        cellCount += 1;
        countMatrixStream.read(reinterpret_cast<char*>(cell.data()), elSize * numTxps);
        if (std::accumulate(cell.begin(), cell.end(), 0.0) == 0){
          std::cout<<"ERROR: Importing counts from binary\n"
                   <<"Cell has 0 reads, should not happen.\n"
                   <<"Saw total "<< cellCount << " Cells before Error"
                   <<std::flush;
          exit(1);
        }
      }
    }

    template <typename ProtocolT>
    bool performWhitelisting(AlevinOpts<ProtocolT>& aopt,
                             std::vector<uint32_t>& umiCount,
                             DoubleMatrixT& countMatrix,
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

      DoubleMatrixT geneCountsMatrix( numCells, DoubleVectorT (numGenes, 0.0));

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
              aopt.jointLog->error("{} mitorna gene not found in txp tp gene map", gene);
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
      // TODO:: This can be parallelized
      tbb::task_scheduler_init tbbScheduler(aopt.numThreads);
      tbb::parallel_for(
                        BlockedIndexRange(size_t(0), size_t(trueBarcodes.size())),
                        [&featureCountsMatrix, &trueBarcodes,
                         &freqCounter, &aopt, &umiCount, &geneCountsMatrix,
                         &countMatrix, &txpToGeneMap, useMitoRna, useRiboRna,
                         &mRnaGenes, &rRnaGenes, numGenes](const BlockedIndexRange& range) -> void {
                          for (auto i : boost::irange(range.begin(), range.end())) {
                            std::vector<double>& featureVector = featureCountsMatrix[i];
                            std::string& currBarcodeName = trueBarcodes[i];
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
                            featureVector[0] = rawBarcodeFrequency ?
                              umiCount[i] / static_cast<double>(rawBarcodeFrequency) : 0.0;

                            double totalCellCount{0.0}, maxCount{0};
                            std::vector<double>& geneCounts = geneCountsMatrix[i];
                            auto& countVec = countMatrix[i];
                            double mitoCount {0.0}, riboCount {0.0};

                            for (size_t j=0; j<countVec.size(); j++){
                              auto count = countVec[j];
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
      aopt.jointLog->info("Starting to Make correlation Matrix (Heavy!)");

      size_t numTrueCells = ( numCells - numLowConfidentBarcode ) / 2;
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

      aopt.jointLog->info("Done making feature Matrix");

      std::vector<uint32_t> whitelistBarcodes =
        classifyBarcodes(featureCountsMatrix, numCells,
                         numFeatures, numLowConfidentBarcode,
                         numTrueCells);

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
                                      DoubleMatrixT& countMatrix,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::InDrop>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      DoubleMatrixT& countMatrix,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Chromium>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      DoubleMatrixT& countMatrix,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      size_t numLowConfidentBarcode);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Custom>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      DoubleMatrixT& countMatrix,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      size_t numLowConfidentBarcode);
  }
}
