#include "CollapsedCellOptimizer.hpp"
#include "EMUtils.hpp"
#include <assert.h>

CollapsedCellOptimizer::CollapsedCellOptimizer() {}
/*
 * Use the "relax" EM algorithm over gene equivalence
 * classes to estimate the latent variables (alphaOut)
 * given the current estimates (alphaIn).
 */
void CellEMUpdate_(std::vector<SalmonEqClass>& eqVec,
                   const CollapsedCellOptimizer::SerialVecType& alphaIn,
                   CollapsedCellOptimizer::SerialVecType& alphaOut) {
  assert(alphaIn.size() == alphaOut.size());

  for (size_t eqID=0; eqID < eqVec.size(); eqID++) {
    auto& kv = eqVec[eqID];

    uint32_t count = kv.count;
    // for each label in this class
    const std::vector<uint32_t>& genes = kv.labels;
    size_t groupSize = genes.size();

    if (BOOST_LIKELY(groupSize > 1)) {
      double denom = 0.0;
      for (size_t i = 0; i < groupSize; ++i) {
        auto gid = genes[i];
        denom += alphaIn[gid];
      }

      if (denom > 0.0) {
        double invDenom = count / denom;
        for (size_t i = 0; i < groupSize; ++i) {
          auto gid = genes[i];
          double v = alphaIn[gid];
          if (!std::isnan(v)) {
            alphaOut[gid] += v * invDenom;
          }
        }//end-for
      }//endif for denom>0
    }//end-if boost gsize>1
    else if (groupSize == 1){
      alphaOut[genes.front()] += count;
    }
    else{
      std::cerr<<"0 Group size for salmonEqclasses in EM\n"
               <<"Please report this on github";
      exit(1);
    }
  }//end-outer for
}


double truncateAlphas(VecT& alphas, double cutoff) {
  // Truncate tiny expression values
  double alphaSum = 0.0;

  for (size_t i = 0; i < alphas.size(); ++i) {
    if (alphas[i] < cutoff) {
      alphas[i] = 0.0;
    }
    alphaSum += alphas[i];
  }
  return alphaSum;
}

bool runPerCellEM(double& totalNumFrags, size_t numGenes,
                  CollapsedCellOptimizer::SerialVecType& alphas,
                  std::vector<SalmonEqClass>& salmonEqclasses,
                  std::shared_ptr<spdlog::logger>& jointlog,
                  bool initUniform){

  // An EM termination criterion, adopted from Bray et al. 2016
  uint32_t minIter {50};
  double relDiffTolerance {0.01};
  uint32_t maxIter {10000};
  size_t numClasses = salmonEqclasses.size();

  if ( initUniform ) {
    double uniformPrior = 1.0/numGenes;
    std::fill(alphas.begin(), alphas.end(), uniformPrior);
  }
  CollapsedCellOptimizer::SerialVecType alphasPrime(numGenes, 0.0);

  assert( numGenes == alphas.size() );
  for (size_t i = 0; i < numGenes; ++i) {
    alphas[i] += 0.5;
    alphas[i] *= 1e-3;
  }

  bool converged{false};
  double maxRelDiff = -std::numeric_limits<double>::max();
  size_t itNum = 0;

  // EM termination criteria, adopted from Bray et al. 2016
  double minAlpha = 1e-8;
  double alphaCheckCutoff = 1e-2;
  constexpr double minWeight = std::numeric_limits<double>::denorm_min();

  while (itNum < minIter or (itNum < maxIter and !converged)) {
    CellEMUpdate_(salmonEqclasses, alphas, alphasPrime);

    converged = true;
    maxRelDiff = -std::numeric_limits<double>::max();
    for (size_t i = 0; i < numGenes; ++i) {
      if (alphasPrime[i] > alphaCheckCutoff) {
        double relDiff =
          std::abs(alphas[i] - alphasPrime[i]) / alphasPrime[i];
        maxRelDiff = (relDiff > maxRelDiff) ? relDiff : maxRelDiff;
        if (relDiff > relDiffTolerance) {
          converged = false;
        }
      }
      alphas[i] = alphasPrime[i];
      alphasPrime[i] = 0.0;
    }

    ++itNum;
  }

  // Truncate tiny expression values
  totalNumFrags = truncateAlphas(alphas, minAlpha);

  if (totalNumFrags < minWeight) {
    jointlog->error("Total alpha weight was too small! "
                    "Make sure you ran salmon correctly.");
    jointlog->flush();
    return false;
  }

  return true;
}

bool runBootstraps(size_t numGenes,
                   CollapsedCellOptimizer::SerialVecType& geneAlphas,
                   std::vector<SalmonEqClass>& salmonEqclasses,
                   std::shared_ptr<spdlog::logger>& jointlog,
                   uint32_t numBootstraps,
                   CollapsedCellOptimizer::SerialVecType& variance,
                   bool useAllBootstraps,
                   std::vector<std::vector<double>>& sampleEstimates,
                   bool initUniform){

  // An EM termination criterion, adopted from Bray et al. 2016
  uint32_t minIter {50};
  double relDiffTolerance {0.01};
  uint32_t maxIter {10000};
  size_t numClasses = salmonEqclasses.size();

  CollapsedCellOptimizer::SerialVecType mean(numGenes, 0.0);
  CollapsedCellOptimizer::SerialVecType squareMean(numGenes, 0.0);
  CollapsedCellOptimizer::SerialVecType alphas(numGenes, 0.0);
  CollapsedCellOptimizer::SerialVecType alphasPrime(numGenes, 0.0);

  //extracting weight of eqclasses for making discrete distribution
  uint32_t totalNumFrags = 0;
  std::vector<uint64_t> eqCounts;
  for (auto& eqclass: salmonEqclasses) {
    totalNumFrags += eqclass.count;
    eqCounts.emplace_back(eqclass.count);
  }

  // Multinomial Sampler
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<uint64_t> csamp(eqCounts.begin(),
                                             eqCounts.end());

  uint32_t bsNum {0};
  while ( bsNum++ < numBootstraps) {
    csamp.reset();

    for (size_t sc = 0; sc < numClasses; ++sc) {
      salmonEqclasses[sc].count = 0;
    }

    for (size_t fn = 0; fn < totalNumFrags; ++fn) {
      salmonEqclasses[csamp(gen)].count += 1;
    }

    for (size_t i = 0; i < numGenes; ++i) {
      if ( initUniform ) {
        alphas[i] = 1.0 / numGenes;
      } else {
        alphas[i] = (geneAlphas[i] + 0.5) * 1e-3;
      }
    }

    bool converged{false};
    double maxRelDiff = -std::numeric_limits<double>::max();
    size_t itNum = 0;

    // EM termination criteria, adopted from Bray et al. 2016
    double minAlpha = 1e-8;
    double alphaCheckCutoff = 1e-2;
    constexpr double minWeight = std::numeric_limits<double>::denorm_min();

    while (itNum < minIter or (itNum < maxIter and !converged)) {
      CellEMUpdate_(salmonEqclasses, alphas, alphasPrime);

      converged = true;
      maxRelDiff = -std::numeric_limits<double>::max();
      for (size_t i = 0; i < numGenes; ++i) {
        if (alphasPrime[i] > alphaCheckCutoff) {
          double relDiff =
            std::abs(alphas[i] - alphasPrime[i]) / alphasPrime[i];
          maxRelDiff = (relDiff > maxRelDiff) ? relDiff : maxRelDiff;
          if (relDiff > relDiffTolerance) {
            converged = false;
          }
        }
        alphas[i] = alphasPrime[i];
        alphasPrime[i] = 0.0;
      }

      ++itNum;
    }//end-EM-while

    // Truncate tiny expression values
    double alphaSum = 0.0;
    // Truncate tiny expression values
    alphaSum = truncateAlphas(alphas, minAlpha);

    if (alphaSum < minWeight) {
      jointlog->error("Total alpha weight was too small! "
                      "Make sure you ran salmon correclty.");
      jointlog->flush();
      return false;
    }

    for(size_t i=0; i<numGenes; i++) {
      double alpha = alphas[i];
      mean[i] += alpha;
      squareMean[i] += alpha * alpha;
    }

    if (useAllBootstraps) {
      sampleEstimates.emplace_back(alphas);
    }
  }//end-boot-while

  // calculate mean and variance of the values
  for(size_t i=0; i<numGenes; i++) {
    double meanAlpha = mean[i] / numBootstraps;
    geneAlphas[i] = meanAlpha;
    variance[i] = (squareMean[i]/numBootstraps) - (meanAlpha*meanAlpha);
  }

  return true;
}

void optimizeCell(std::vector<std::string>& trueBarcodes,
                  std::atomic<uint32_t>& barcode,
                  size_t totalCells, eqMapT& eqMap,
                  std::deque<std::pair<TranscriptGroup, uint32_t>>& orderedTgroup,
                  std::shared_ptr<spdlog::logger>& jointlog,
                  bfs::path& outDir, std::vector<uint32_t>& umiCount,
                  std::vector<CellState>& skippedCB,
                  bool verbose, GZipWriter& gzw, size_t umiLength, bool noEM,
                  bool quiet, tbb::atomic<double>& totalDedupCounts,
                  tbb::atomic<uint32_t>& totalExpGeneCounts,
                  spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                  uint32_t numGenes, uint32_t numBootstraps,
                  bool naiveEqclass, bool dumpUmiGraph, bool useAllBootstraps,
                  bool initUniform, CFreqMapT& freqCounter,
                  spp::sparse_hash_set<uint32_t>& mRnaGenes,
                  spp::sparse_hash_set<uint32_t>& rRnaGenes,
                  std::atomic<uint64_t>& totalUniEdgesCounts,
                  std::atomic<uint64_t>& totalBiEdgesCounts){
  size_t numCells {trueBarcodes.size()};
  size_t trueBarcodeIdx;

  // looping over until all the cells
  while((trueBarcodeIdx = barcode++) < totalCells) {
    // per-cell level optimization
    if ( umiCount[trueBarcodeIdx] == 0 ) {
      //skip the barcode if no mapped UMI
      skippedCB[trueBarcodeIdx].inActive = true;
      continue;
    }

    // extracting the sequence of the barcode
    auto& trueBarcodeStr = trueBarcodes[trueBarcodeIdx];

    //extracting per-cell level eq class information
    double totalCount{0.0};
    double totalExpGenes{0};
    std::vector<uint32_t> eqIDs;
    std::vector<uint32_t> eqCounts;
    std::vector<UGroupT> umiGroups;
    std::vector<tgrouplabelt> txpGroups;
    std::vector<double> geneAlphas(numGenes, 0.0);
    std::vector<uint8_t> tiers (numGenes, 0);

    for (auto& key : orderedTgroup) {
      //traversing each class and copying relevant data.
      bool isKeyPresent = eqMap.find_fn(key.first, [&](const SCTGValue& val){
          auto& bg = val.barcodeGroup;
          auto bcIt = bg.find(trueBarcodeIdx);

          // sub-selecting bgroup of this barcode only
          if (bcIt != bg.end()){
            // extracting txp labels
            const std::vector<uint32_t>& txps = key.first.txps;

            // original counts of the UMI
            uint32_t eqCount {0};
            for(auto& ugroup: bcIt->second){
              eqCount += ugroup.second;
            }

            txpGroups.emplace_back(txps);
            umiGroups.emplace_back(bcIt->second);

            // for dumping per-cell eqclass vector
            if(verbose){
              eqIDs.push_back(static_cast<uint32_t>(key.second));
              eqCounts.push_back(eqCount);
            }
          }
        });

      if(!isKeyPresent){
        jointlog->error("Not able to find key in Cuckoo hash map."
                        "Please Report this issue on github");
        jointlog->flush();
        exit(1);
      }
    }

    if ( !naiveEqclass ) {
      // perform the UMI deduplication step
      std::vector<SalmonEqClass> salmonEqclasses;
      spp::sparse_hash_map<uint16_t, uint32_t> numMolHash;
      bool dedupOk = dedupClasses(geneAlphas, totalCount, txpGroups,
                                  umiGroups, salmonEqclasses,
                                  txpToGeneMap, tiers, gzw,
                                  dumpUmiGraph, trueBarcodeStr, numMolHash,
                                  totalUniEdgesCounts, totalBiEdgesCounts);
      if( !dedupOk ){
        jointlog->error("Deduplication for cell {} failed \n"
                        "Please Report this on github.", trueBarcodeStr);
        jointlog->flush();
        std::exit(74);
      }

      if ( numBootstraps and noEM ) {
        jointlog->error("Cannot perform bootstrapping with noEM");
        jointlog->flush();
        exit(1);
      }

      // perform EM for resolving ambiguity
      if ( !noEM ) {
        bool isEMok = runPerCellEM(totalCount,
                                   numGenes,
                                   geneAlphas,
                                   salmonEqclasses,
                                   jointlog,
                                   initUniform);
        if( !isEMok ){
          jointlog->error("EM iteration for cell {} failed \n"
                          "Please Report this on github.", trueBarcodeStr);
          jointlog->flush();
          std::exit(74);
        }
      }

      std::string features;
      uint8_t featureCode {0};
      {
        std::stringstream featuresStream;
        featuresStream << trueBarcodeStr;

        // Making features
        double totalUmiCount {0.0};
        double maxNumUmi {0.0};
        for (auto count: geneAlphas) {
          if (count>0.0) {
            totalUmiCount += count;
            totalExpGenes += 1;
            if (count > maxNumUmi) { maxNumUmi = count; }
          }
        }

        uint32_t numGenesOverMean {0};
        double mitoCount {0.0}, riboCount {0.0};
        double meanNumUmi {totalUmiCount / totalExpGenes};
        double meanByMax = maxNumUmi ? meanNumUmi / maxNumUmi : 0.0;
        for (size_t j=0; j<geneAlphas.size(); j++){
          auto count = geneAlphas[j];
          if (count > meanNumUmi) { ++numGenesOverMean; }

          if (mRnaGenes.contains(j)){
            mitoCount += count;
          }
          if (rRnaGenes.contains(j)){
            riboCount += count;
          }
        }

        auto indexIt = freqCounter.find(trueBarcodeStr);
        bool indexOk = indexIt != freqCounter.end();
        if ( not indexOk ){
          jointlog->error("Error: index {} not found in freq Counter\n"
                          "Please Report the issue on github", trueBarcodeStr);
          jointlog->flush();
          exit(84);
        }

        uint64_t numRawReads = *indexIt;
        uint64_t numMappedReads { umiCount[trueBarcodeIdx] };
        double mappingRate = numRawReads ?
          numMappedReads / static_cast<double>(numRawReads) : 0.0;
        double deduplicationRate = numMappedReads ?
          1.0 - (totalUmiCount / numMappedReads) : 0.0;

        // Feature created after discussion with Mehrtash
        double averageNumMolPerArbo {0.0};
        size_t totalNumArborescence {0};
        std::stringstream arboString ;
        for (auto& it: numMolHash) {
          totalNumArborescence += it.second;
          averageNumMolPerArbo += (it.first * it.second);
          if (dumpUmiGraph) {
            arboString << "\t" << it.first << ":" << it.second;
          }
        }
        averageNumMolPerArbo /= totalNumArborescence;

        featuresStream << "\t" << numRawReads
                       << "\t" << numMappedReads
                       << "\t" << totalUmiCount
                       << "\t" << mappingRate
                       << "\t" << deduplicationRate
                       << "\t" << meanByMax
                       << "\t" << totalExpGenes
                       << "\t" << numGenesOverMean;

        if (dumpUmiGraph) {
          featuresStream << arboString.rdbuf();
        } else {
          featuresStream << "\t" << averageNumMolPerArbo;
        }

        if (mRnaGenes.size() > 1) {
          featureCode += 1;
          featuresStream << "\t" << mitoCount / totalUmiCount;
        }

        if (rRnaGenes.size() > 1) {
          featureCode += 2;
          featuresStream << "\t" << riboCount / totalUmiCount;
        }

        features = featuresStream.str();
      } // end making features


      // write the abundance for the cell
      gzw.writeSparseAbundances( trueBarcodeStr,
                                 features,
                                 featureCode,
                                 geneAlphas,
                                 tiers,
                                 dumpUmiGraph );


      // maintaining count for total number of predicted UMI
      salmon::utils::incLoop(totalDedupCounts, totalCount);
      totalExpGeneCounts += totalExpGenes;

      if ( numBootstraps > 0 ){
        std::vector<std::vector<double>> sampleEstimates;
        std::vector<double> bootVariance(numGenes, 0.0);

        bool isBootstrappingOk = runBootstraps(numGenes,
                                               geneAlphas,
                                               salmonEqclasses,
                                               jointlog,
                                               numBootstraps,
                                               bootVariance,
                                               useAllBootstraps,
                                               sampleEstimates,
                                               initUniform);
        if( not isBootstrappingOk or
            (useAllBootstraps and sampleEstimates.size()!=numBootstraps)
            ){
          jointlog->error("Bootstrapping failed \n"
                          "Please Report this on github.");
          jointlog->flush();
          std::exit(74);
        }

        // write the abundance for the cell
        gzw.writeBootstraps( trueBarcodeStr,
                             geneAlphas, bootVariance,
                             useAllBootstraps, sampleEstimates);
      }//end-if
    }
    else {
      // doing per eqclass level naive deduplication
      for (size_t eqId=0; eqId<umiGroups.size(); eqId++) {
        spp::sparse_hash_set<uint64_t> umis;

        for(auto& it: umiGroups[eqId]) {
          umis.insert( it.first );
        }
        totalCount += umis.size();

        // filling in the eqclass level deduplicated counts
        if (verbose) {
          eqCounts[eqId] = umis.size();
        }
      }

      // maintaining count for total number of predicted UMI
      salmon::utils::incLoop(totalDedupCounts, totalCount);
    }

    if (verbose) {
      gzw.writeCellEQVec(trueBarcodeIdx, eqIDs, eqCounts, true);
    }

    //printing on screen progress
    const char RESET_COLOR[] = "\x1b[0m";
    char green[] = "\x1b[30m";
    green[3] = '0' + static_cast<char>(fmt::GREEN);
    char red[] = "\x1b[30m";
    red[3] = '0' + static_cast<char>(fmt::RED);

    double cellCount {static_cast<double>(barcode)};//numCells-jqueue.size_approx()};
    if (cellCount > totalCells) { cellCount = totalCells; }
    double percentCompletion {cellCount*100/numCells};
    if (not quiet){
      fmt::print(stderr, "\033[A\r\r{}Analyzed {} cells ({}{}%{} of all).{}\n",
                 green, cellCount, red, round(percentCompletion), green, RESET_COLOR);
    }
  }
}

template <typename ProtocolT>
bool CollapsedCellOptimizer::optimize(EqMapT& fullEqMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      AlevinOpts<ProtocolT>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode){
  size_t numCells = trueBarcodes.size();
  size_t numGenes = geneIdxMap.size();
  size_t numWorkerThreads{1};
  bool hasWhitelist = boost::filesystem::exists(aopt.whitelistFile);

  if (aopt.numThreads > 1) {
    numWorkerThreads = aopt.numThreads - 1;
  }

  //get the keys of the map
  std::deque<std::pair<TranscriptGroup, uint32_t>> orderedTgroup;
  //spp::sparse_hash_set<uint64_t> uniqueUmisCounter;
  uint32_t eqId{0};
  for(const auto& kv : fullEqMap.lock_table()){
    // assuming the iteration through lock table is always same
    if(kv.first.txps.size() == 1){
      orderedTgroup.push_front(std::make_pair(kv.first, eqId));
    }
    else{
      orderedTgroup.push_back(std::make_pair(kv.first, eqId));
    }
    eqId++;
  }

  if (aopt.noEM) {
    aopt.jointLog->warn("Not performing EM; this may result in discarding ambiguous reads\n");
    aopt.jointLog->flush();
  }

  if (aopt.initUniform) {
    aopt.jointLog->warn("Using uniform initialization for EM");
    aopt.jointLog->flush();
  }

  spp::sparse_hash_set<uint32_t> mRnaGenes, rRnaGenes;
  bool useMito {false}, useRibo {false};
  if(boost::filesystem::exists(aopt.mRnaFile)) {
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
    if (mRnaGenes.size() > 0 ) { useMito = true; }
  }
  else if (hasWhitelist) {
    aopt.jointLog->warn("mrna file not provided; using is 1 less feature for whitelisting");
  }

  if(boost::filesystem::exists(aopt.rRnaFile)){
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
    if (rRnaGenes.size() > 0 ) { useRibo = true; }
  }
  else if (hasWhitelist) {
    aopt.jointLog->warn("rrna file not provided; using is 1 less feature for whitelisting");
  }


  std::vector<CellState> skippedCB (numCells);
  std::atomic<uint32_t> bcount{0};
  tbb::atomic<double> totalDedupCounts{0.0};
  tbb::atomic<uint32_t> totalExpGeneCounts{0};
  std::atomic<uint64_t> totalBiEdgesCounts{0};
  std::atomic<uint64_t> totalUniEdgesCounts{0};

  std::vector<std::thread> workerThreads;
  for (size_t tn = 0; tn < numWorkerThreads; ++tn) {
    workerThreads.emplace_back(optimizeCell,
                               std::ref(trueBarcodes),
                               std::ref(bcount),
                               numCells,
                               std::ref(fullEqMap),
                               std::ref(orderedTgroup),
                               std::ref(aopt.jointLog),
                               std::ref(aopt.outputDirectory),
                               std::ref(umiCount),
                               std::ref(skippedCB),
                               aopt.dumpBarcodeEq,
                               std::ref(gzw),
                               aopt.protocol.umiLength,
                               aopt.noEM,
                               aopt.quiet,
                               std::ref(totalDedupCounts),
                               std::ref(totalExpGeneCounts),
                               std::ref(txpToGeneMap),
                               numGenes,
                               aopt.numBootstraps,
                               aopt.naiveEqclass,
                               aopt.dumpUmiGraph,
                               aopt.dumpfeatures,
                               aopt.initUniform,
                               std::ref(freqCounter),
                               std::ref(rRnaGenes),
                               std::ref(mRnaGenes),
                               std::ref(totalUniEdgesCounts),
                               std::ref(totalBiEdgesCounts));
  }

  for (auto& t : workerThreads) {
    t.join();
  }
  aopt.jointLog->info("Total {0:.2f} UMI after deduplicating.",
                      totalDedupCounts);
  aopt.jointLog->info("Total {} BiDirected Edges.",
                      totalBiEdgesCounts);
  aopt.jointLog->info("Total {} UniDirected Edges.",
                      totalUniEdgesCounts);

  //adjusting for float
  aopt.totalDedupUMIs = totalDedupCounts+1;
  aopt.totalExpGenes = totalExpGeneCounts;

  uint32_t skippedCBcount {0};
  for(auto cb: skippedCB){
    if (cb.inActive) {
      skippedCBcount += 1;
    }
  }

  if( skippedCBcount > 0 ) {
    aopt.jointLog->warn("Skipped {} barcodes due to No mapped read",
                        skippedCBcount);
    auto lowRegionCutoffIdx = numCells - numLowConfidentBarcode;
    for (size_t idx=0; idx < numCells; idx++){
      // not very efficient way but assuming the size is small enough
      if (skippedCB[idx].inActive) {
        trueBarcodes.erase(trueBarcodes.begin() + idx);
        if (idx > lowRegionCutoffIdx){
          numLowConfidentBarcode--;
        }
      }
    }
    numCells = trueBarcodes.size();
  }

  gzw.close_all_streams();

  std::vector<std::string> geneNames(numGenes);
  for (auto geneIdx : geneIdxMap) {
    geneNames[geneIdx.second] = geneIdx.first;
  }
  boost::filesystem::path gFilePath = aopt.outputDirectory / "quants_mat_cols.txt";
  std::ofstream gFile(gFilePath.string());
  std::ostream_iterator<std::string> giterator(gFile, "\n");
  std::copy(geneNames.begin(), geneNames.end(), giterator);
  gFile.close();

  if( not hasWhitelist ){
    aopt.jointLog->info("Clearing EqMap; Might take some time.");
    fullEqMap.clear();

    aopt.jointLog->info("Starting white listing");
    bool whitelistingSuccess = alevin::whitelist::performWhitelisting(aopt,
                                                                      umiCount,
                                                                      trueBarcodes,
                                                                      freqCounter,
                                                                      useRibo,
                                                                      useMito,
                                                                      numLowConfidentBarcode);
    if (!whitelistingSuccess) {
      aopt.jointLog->error(
                           "The white listing algorithm failed. This is likely the result of "
                           "bad input (or a bug). If you cannot track down the cause, please "
                           "report this issue on GitHub.");
      aopt.jointLog->flush();
      return false;
    }

    aopt.jointLog->info("Finished white listing");
    aopt.jointLog->flush();
  } //end-if whitelisting

  if (aopt.dumpMtx){
    aopt.jointLog->info("Starting dumping cell v gene counts in mtx format");
    boost::filesystem::path qFilePath = aopt.outputDirectory / "quants_mat.mtx.gz";

    boost::iostreams::filtering_ostream qFile;
    qFile.push(boost::iostreams::gzip_compressor(6));
    qFile.push(boost::iostreams::file_sink(qFilePath.string(),
                                           std::ios_base::out | std::ios_base::binary));

    // mtx header
    qFile << "%%MatrixMarket\tmatrix\tcoordinate\treal\tgeneral" << std::endl
          << numCells << "\t" << numGenes << "\t" << totalExpGeneCounts << std::endl;

    {

      auto popcount = [](uint8_t n) {
        size_t count {0};
        while (n) {
          n &= n-1;
          ++count;
        }
        return count;
      };

      uint32_t zerod_cells {0};
      size_t numFlags = (numGenes/8)+1;
      std::vector<uint8_t> alphasFlag (numFlags, 0);
      size_t flagSize = sizeof(decltype(alphasFlag)::value_type);

      std::vector<float> alphasSparse;
      alphasSparse.reserve(numFlags/2);
      size_t elSize = sizeof(decltype(alphasSparse)::value_type);

      auto countMatFilename = aopt.outputDirectory / "quants_mat.gz";
      if(not boost::filesystem::exists(countMatFilename)){
        std::cout<<"ERROR: Can't import Binary file quants.mat.gz, it doesn't exist" << std::flush;
        exit(84);
      }

      boost::iostreams::filtering_istream countMatrixStream;
      countMatrixStream.push(boost::iostreams::gzip_decompressor());
      countMatrixStream.push(boost::iostreams::file_source(countMatFilename.string(),
                                                           std::ios_base::in | std::ios_base::binary));

      for (size_t cellCount=0; cellCount<numCells; cellCount++){
        countMatrixStream.read(reinterpret_cast<char*>(alphasFlag.data()), flagSize * numFlags);

        size_t numExpGenes {0};
        std::vector<size_t> indices;
        for (size_t j=0; j<alphasFlag.size(); j++) {
          uint8_t flag = alphasFlag[j];
          size_t numNonZeros = popcount(flag);
          numExpGenes += numNonZeros;

          for (size_t i=0; i<8; i++){
            if (flag & (128 >> i)) {
              indices.emplace_back( (i*8)+j );
            }
          }
        }

        if (indices.size() != numExpGenes) {
          aopt.jointLog->error("binary format reading error {}: {}: {}",
                               indices.size(), numExpGenes);
          aopt.jointLog->flush();
          exit(84);
        }


        alphasSparse.clear();
        alphasSparse.resize(numExpGenes);
        countMatrixStream.read(reinterpret_cast<char*>(alphasSparse.data()), elSize * numExpGenes);

        float readCount {0.0};
        readCount += std::accumulate(alphasSparse.begin(), alphasSparse.end(), 0.0);

        for(size_t i=0; i<numExpGenes; i++) {
          qFile << cellCount+1 << "\t"
                << indices[i] << "\t"
                << alphasSparse[i] <<  std::endl;
        }

        //size_t alphasSparseCounter {0};
        //for (size_t i=0; i<numGenes; i+=8) {
        //  uint8_t flag = alphasFlag[i];
        //  for (size_t j=0; j<8; j++) {
        //    size_t vectorIndex = i+j;
        //    if (vectorIndex >= numGenes) { break; }

        //    if ( flag & (1<<(7-j)) ) {
        //      if (alphasSparseCounter >= numExpGenes) {
        //        aopt.jointLog->error("binary format reading error {}: {}: {}",
        //                             alphasSparseCounter, numExpGenes, readCount);
        //        aopt.jointLog->flush();
        //        exit(84);
        //      }

        //      float count = alphasSparse[alphasSparseCounter];
        //      readCount += count;
        //      qFile << cellCount+1 << "\t"
        //            << vectorIndex+1 << "\t"
        //            << count << std::endl;

        //      ++alphasSparseCounter;
        //    }
        //  }
        //}

        if (readCount == 0.0){
          zerod_cells += 1;
        }
      } // end-for each cell

      if (zerod_cells > 0) {
        aopt.jointLog->warn("Found {} cells with 0 counts", zerod_cells);
      }
    }

    boost::iostreams::close(qFile);
    aopt.jointLog->info("Finished dumping counts into mtx");
  }

  return true;
} //end-optimize


namespace apt = alevin::protocols;
template
bool CollapsedCellOptimizer::optimize(EqMapT& fullEqMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                      AlevinOpts<apt::DropSeq>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode);
template
bool CollapsedCellOptimizer::optimize(EqMapT& fullEqMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,

                                      AlevinOpts<apt::InDrop>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode);
template
bool CollapsedCellOptimizer::optimize(EqMapT& fullEqMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,

                                      AlevinOpts<apt::ChromiumV3>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode);
template
bool CollapsedCellOptimizer::optimize(EqMapT& fullEqMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,

                                      AlevinOpts<apt::Chromium>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode);
template
bool CollapsedCellOptimizer::optimize(EqMapT& fullEqMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,

                                      AlevinOpts<apt::Gemcode>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode);
template
bool CollapsedCellOptimizer::optimize(EqMapT& fullEqMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,

                                      AlevinOpts<apt::CELSeq>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode);
template
bool CollapsedCellOptimizer::optimize(EqMapT& fullEqMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,

                                      AlevinOpts<apt::CELSeq2>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode);
template
bool CollapsedCellOptimizer::optimize(EqMapT& fullEqMap,
                                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,

                                      AlevinOpts<apt::Custom>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode);
