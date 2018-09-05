#include "CollapsedCellOptimizer.hpp"
#include "EMUtils.hpp"

CollapsedCellOptimizer::CollapsedCellOptimizer() {}

double truncateAlphas(VecT& alphas, double cutoff) {
  // Truncate tiny expression values
  double alphaSum = 0.0;

  for (size_t i = 0; i < alphas.size(); ++i) {
    if (alphas[i] <= cutoff) {
      alphas[i] = 0.0;
    }
    alphaSum += alphas[i];
  }
  return alphaSum;
}

bool runPerCellEM(
                  std::vector<std::vector<uint32_t>>& txpGroups,
                  std::vector<std::vector<double>>& txpGroupCombinedWeights,
                  std::vector<uint64_t>& txpGroupCounts,
                  const std::vector<Transcript>& transcripts,
                  uint64_t totalNumFrags,
                  CollapsedCellOptimizer::SerialVecType& alphas,
                  std::shared_ptr<spdlog::logger>& jointlog,
                  std::unordered_set<uint32_t>& activeTranscriptIds){

  // An EM termination criterion, adopted from Bray et al. 2016
  uint32_t minIter {50};
  double relDiffTolerance {0.01};
  uint32_t maxIter {10000};
  size_t numClasses = txpGroups.size();
  double uniformTxpWeight = 1.0 / activeTranscriptIds.size();

  std::vector<uint64_t> eqCounts(numClasses, 0);
  CollapsedCellOptimizer::SerialVecType alphasPrime(transcripts.size(), 0.0);

  for (size_t i = 0; i < transcripts.size(); ++i) {
    auto got = activeTranscriptIds.find( i );
    if ( got == activeTranscriptIds.end() ){
      alphas[i] = 0.0;
    }
    else{
      alphas[i] = uniformTxpWeight * totalNumFrags;
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

    EMUpdate_(txpGroups, txpGroupCombinedWeights,
              txpGroupCounts, const_cast<std::vector<Transcript>&>(transcripts),
              alphas, alphasPrime);

    converged = true;
    maxRelDiff = -std::numeric_limits<double>::max();
    for (size_t i = 0; i < transcripts.size(); ++i) {
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
  double alphaSum = 0.0;
  // Truncate tiny expression values
  alphaSum = truncateAlphas(alphas, minAlpha);

  if (alphaSum < minWeight) {
    jointlog->error("Total alpha weight was too small! "
                    "Make sure you ran salmon correclty.");
    jointlog->flush();
    return false;
  }

  return true;
}

void getMinSetTxps(std::vector<tgrouplabelt>& txpgroups,
                   std::vector<UGroupT>& umigroups,
                   size_t numTxps,
                   spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap){
  // get left and right set from Rob's figure and frequency of each txp
  std::vector<uint32_t> txpCount(numTxps);
  spp::sparse_hash_set<uint32_t> tgroupSet;
  for(size_t i=0; i<txpgroups.size(); i++){
    for (auto txp: txpgroups[i]){
      txpCount[txp]++;
    }
    tgroupSet.insert(i);
  }

  // initialize original index locations
  std::vector<size_t> idx(numTxps);
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in txpCount
  sort(idx.begin(), idx.end(),
       [&txpCount](size_t i1, size_t i2) {return txpCount[i1] > txpCount[i2];});

  size_t head{0};
  spp::sparse_hash_set<uint32_t> minTxps;
  while(tgroupSet.size() != 0){
    // extract transcript id for most frequent txp
    uint32_t tid = idx[head];
    std::vector<uint32_t> removableTgroups;

    // select which tgroup to remove
    for (auto tgroupIdx: tgroupSet){
      auto& tgroup = txpgroups[tgroupIdx];
      if (std::find(tgroup.begin(), tgroup.end(), tid) != tgroup.end() ){
        removableTgroups.emplace_back(tgroupIdx);
      }
    }

    // remove tgroups relevant for this tid
    for(auto tgroupIdx: removableTgroups){
      tgroupSet.erase(tgroupIdx);
    }

    // if removed at least one tgroup then add to mintxps
    if(removableTgroups.size() > 0){
      minTxps.insert(tid);
    }

    // get next highest txp
    head++;
  }

  // modify txp labels of eqclasses
  // Very inefficient way O(n^2) but doing for checking
  std::deque<spp::sparse_hash_set<uint32_t>> newTgroups;
  std::deque<UGroupT> newUgroups;
  for (size_t i=0; i<txpgroups.size(); i++){
    auto& tgroup = txpgroups[i];
    // extract new label by removing txps not in mintxps
    spp::sparse_hash_set<uint32_t> newlabel;
    spp::sparse_hash_set<uint32_t> geneSet;
    for(auto txp: tgroup){
      if (minTxps.contains(txp)){
        newlabel.insert(txp);
      }
      if (txpToGeneMap.contains(txp)){
        geneSet.insert(txpToGeneMap[txp]);
      }
      else{
        std::cerr<< "Wrong Txp To gene Map"<<std::flush;
        exit(1);
      }
    }

    //TODO: Major Update required here when doing ambiguity resolution
    if (geneSet.size() > 1){
      continue;
    }

    // check if this label is already present in the vector
    auto it = std::find(newTgroups.begin(), newTgroups.end(), newlabel);
    if ( it != newTgroups.end() ){
      // if found then just append UMI
      auto index = std::distance(newTgroups.begin(), it);
      auto& ugroup = newUgroups[index];
      for (auto umi: umigroups[i]){
        ugroup[umi.first] += umi.second;
      }
    }
    else{
      // if 1 length txp insert in front
      if (newlabel.size() == 1){
        newUgroups.push_front(umigroups[i]);
        newTgroups.push_front(newlabel);
      }
      // else insert t last
      else{
        newUgroups.push_back(umigroups[i]);
        newTgroups.push_back(newlabel);
      }
    }
  }

  // replace the present group with a new one
  txpgroups.clear();
  for(auto& tgroupSet: newTgroups){
    std::vector<uint32_t> tgroup (tgroupSet.begin(), tgroupSet.end());
    txpgroups.emplace_back(tgroup);
  }

  umigroups.clear();
  for(auto& ugroup: newUgroups){
    umigroups.emplace_back(ugroup);
  }
}

void optimizeCell(SCExpT& experiment,
                  std::vector<std::string>& trueBarcodes,
                  std::atomic<uint32_t>& barcode,
                  size_t totalCells, eqMapT& eqMap,
                  std::deque<std::pair<TranscriptGroup, uint32_t>>& orderedTgroup,
                  std::shared_ptr<spdlog::logger>& jointlog,
                  bfs::path& outDir, std::vector<uint32_t>& umiCount,
                  std::vector<CellState>& skippedCB,
                  bool verbose, GZipWriter& gzw, size_t umiLength, bool doEM,
                  bool quiet, std::atomic<uint64_t>& totalDedupCounts,
                  spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                  uint32_t numGenes, bool txpLevel, bool naive,
                  bool eqClassLevel, bool inDebugMode){
  // HACK: todo: major refactoring needed
  if (eqClassLevel){
    txpLevel = eqClassLevel;
  }
  size_t numCells {trueBarcodes.size()};
  size_t trueBarcodeIdx;

  // looping over until the end of the file
  while((trueBarcodeIdx = barcode++) < totalCells) {
    // per-cell level optimization
    if ( (not inDebugMode && umiCount[trueBarcodeIdx] < 10) or
         (inDebugMode && umiCount[trueBarcodeIdx] == 0) ) {
      //skip the barcode if no mapped UMI
      skippedCB[trueBarcodeIdx].inActive = true;
      continue;
    }

    auto& trueBarcodeStr = trueBarcodes[trueBarcodeIdx];

    //extracting per-cell level eq class information
    const std::vector<Transcript>& transcripts = experiment.transcripts();
    std::unordered_set<uint32_t> activetranscriptids;
    std::vector<tgrouplabelt> txpgroups;
    std::vector<UGroupT> umigroups;
    std::vector<tgroupweightvec> txpgroupcombinedweights;
    std::vector<uint64_t> origcounts;
    uint64_t totalcount{0};

    spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint64_t>> umiBiasList;

    std::vector<double> geneAlphas(numGenes, 0.0);
    spp::sparse_hash_map<uint32_t, spp::sparse_hash_map<uint64_t, uint32_t>> geneUmiGroup;

    // equivalence class vector encoding for this cell (i.e. row)
    std::vector<uint32_t> eqIDs;
    std::vector<uint32_t> counts;

    for (auto& key : orderedTgroup) {
      //traversing each class and copying relevant data.
      bool isKeyPresent = eqMap.find_fn(key.first, [&](const SCTGValue& val){
          auto& bg = val.barcodeGroup;
          auto bcIt = bg.find(trueBarcodeIdx);

          if (bcIt != bg.end()){
            // sub-selecting bgroup of this barcode only
            const std::vector<uint32_t>& txps = key.first.txps;
            uint32_t eqCount {0};
            for(auto& ugroup: bcIt->second){
                eqCount += ugroup.second;
            }
            // Avi -> Major Hack
            // Basically probStartpos is 1/effec_length.
            // effec_length for all txp is 100
            // Eqclass weights are set to 1 since sopt.noRichEqClasses is true
            // Making aux as 1/#oftxps_in_class
            // Ideally should be taken from eqclass.combinedWeights but
            // ignoring for now
            txpgroups.emplace_back(txps);
            umigroups.emplace_back(bcIt->second);

            if(verbose){
              eqIDs.push_back(static_cast<uint32_t>(key.second));
              counts.push_back(static_cast<uint32_t>(eqCount));
            }
            if (not txpLevel){
              uint32_t geneId;
              if ( alevin::utils::hasOneGene(txps, geneId,
                                             txpToGeneMap, numGenes) ){
                for (auto& umiIt:bcIt->second){
                  geneUmiGroup[geneId][umiIt.first] += umiIt.second;
                }
              }
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

    if (verbose) {
      gzw.writeCellEQVec(trueBarcodeIdx, eqIDs, counts, true);
    }

    if (not txpLevel){
      std::vector<uint32_t> dummyTxpGroup {0,1};
      for (auto& geneUmis: geneUmiGroup){
        auto eqCount = dedupReads(umiLength,
                                  jointlog,
                                  geneUmis.second,
                                  umiBiasList,
                                  dummyTxpGroup);
        totalDedupCounts += eqCount;
        geneAlphas[geneUmis.first] = eqCount;
      }
    }

    if (not eqClassLevel){
      getMinSetTxps(txpgroups, umigroups,
                    transcripts.size(),
                    txpToGeneMap);
    }

    if (txpLevel or doEM){
      // this scope contains code for deduplicating and per-Eq class level
      spp::sparse_hash_set<uint32_t> removableTgroups;
      for (size_t i=0; i<txpgroups.size(); i++){
        // sub-selecting bgroup of this barcode only
        auto eqCount = dedupReads(umiLength,
                                  jointlog,
                                  umigroups[i],
                                  umiBiasList,
                                  txpgroups[i]);

        if ( eqCount != 0 ) {
          const std::vector<uint32_t>& txps = txpgroups[i];
          // Avi -> Major Hack
          // Basically probStartpos is 1/effec_length.
          // effec_length for all txp is 100
          // Eqclass weights are set to 1 since sopt.noRichEqClasses is true
          // Making aux as 1/#oftxps_in_class
          // Ideally should be taken from eqclass.combinedWeights but
          // ignoring for now
          if (doEM){
            const tgroupweightvec auxs(txps.size(), 1.0/txps.size());
            // convert to non-atomic
            txpgroupcombinedweights.emplace_back(auxs.begin(), auxs.end());
          }
          origcounts.emplace_back(eqCount);

          totalcount += eqCount;

          // currently add only 1 length eqclass to active txps
          for (auto& t : txps) {
            activetranscriptids.insert(t);
          }
        }
        else{
          removableTgroups.insert(i);
        }
      }

      // remove 0 count eqclass
      std::vector<tgrouplabelt> newTgroups;
      for(size_t i=0; i<txpgroups.size();i++){
        if (not removableTgroups.contains(i)){
          newTgroups.emplace_back(txpgroups[i]);
        }
      }
      txpgroups = newTgroups;
    }// end-if txpLevel

    if( not doEM ){
      // no em i.e. only eqclass mode
      if ( txpLevel ){
        for (size_t k=0; k<txpgroups.size(); k++){
          uint32_t geneId;
          if ( alevin::utils::hasOneGene(txpgroups[k], geneId,
                                         txpToGeneMap, numGenes) ){
            geneAlphas[geneId] += origcounts[k];
            totalDedupCounts += origcounts[k];
          }
        }
      }
      else if( not naive ) {
        // collision resolution based on the UMI in 1-length
        // eqclasses after minset-ting
        geneUmiGroup.clear();
        for (size_t k=0; k<txpgroups.size(); k++){
          uint32_t geneId;
          // redundant if should be removed later
          if ( alevin::utils::hasOneGene(txpgroups[k], geneId,
                                         txpToGeneMap, numGenes)
               and txpgroups[k].size() == 1 ){
            if ( geneUmiGroup.contains(geneId) ){
              for (auto umiIt: umigroups[k]){
                if (geneUmiGroup[geneId].contains(umiIt.first)){
                  geneAlphas[geneId] += 1;
                  totalDedupCounts += 1;
                }
                else{
                  geneUmiGroup[geneId][umiIt.first] += 1;
                }
              }
            }
            else{
              for (auto umiIt: umigroups[k]){
                geneUmiGroup[geneId][umiIt.first] += 1;
              }
            }
          }
        }
      }
      gzw.writeAbundances(inDebugMode, trueBarcodeStr, geneAlphas);
    }
    else {
      CollapsedCellOptimizer::SerialVecType alphas(transcripts.size(), 0.0);
      bool isEMok = runPerCellEM(txpgroups,
                                 txpgroupcombinedweights,
                                 origcounts,
                                 transcripts,
                                 totalcount,
                                 alphas,
                                 jointlog,
                                 activetranscriptids);
      if( !isEMok ){
        jointlog->error("EM iteration for cell {} failed \n"
                        "Please Report this on github.", trueBarcodeStr);
        jointlog->flush();
        std::exit(1);
      }

      if (txpLevel){
        // dump txp level counts in matrix
        gzw.writeAbundances( inDebugMode, trueBarcodeStr, alphas);
      }
      else{
        //dump gene level counts in matrix
        std::vector<double> geneAlphas(numGenes, 0.0);
        for (size_t tid=0; tid<alphas.size(); tid++){
          uint32_t gid;
          if(txpToGeneMap.contains(tid)){
            gid = txpToGeneMap.at(tid);
          }
          else{
            std::cerr << "Out of Range error for txp to gene Map: " << '\n' << std::flush;
            std::cerr << tid << "\t not found" << std::flush;
            exit(1);
          }
          if (gid > numGenes){
            std::cerr<< "Gene id out of range"
                     << "Please check txp2gene has the write entries"
                     << std::flush;
            exit(1);
          }
          geneAlphas[gid] += alphas[tid];
        }
        gzw.writeAbundances( inDebugMode, trueBarcodeStr, geneAlphas);
      }
      totalDedupCounts += totalcount;
    }

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

uint32_t getTxpToGeneMap(spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                         const std::vector<Transcript>& transcripts,
                         const std::string& geneMapFile,
                         spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap){
  std::string fname = geneMapFile;
  std::ifstream t2gFile(fname);

  spp::sparse_hash_map<std::string, uint32_t> txpIdxMap(transcripts.size());

  for (size_t i=0; i<transcripts.size(); i++){
    txpIdxMap[ transcripts[i].RefName ] = i;
  }

  uint32_t tid, gid, geneCount{0};
  std::string tStr, gStr;
  if(t2gFile.is_open()) {
    while( not t2gFile.eof() ) {
      t2gFile >> tStr >> gStr;

      if(not txpIdxMap.contains(tStr)){
        continue;
      }
      tid = txpIdxMap[tStr];

      if (geneIdxMap.contains(gStr)){
        gid = geneIdxMap[gStr];
      }
      else{
        gid = geneCount;
        geneIdxMap[gStr] = gid;
        geneCount++;
      }

      txpToGeneMap[tid] = gid;
    }
    t2gFile.close();
  }
  if(txpToGeneMap.size() < transcripts.size()){
    std::cerr << "ERROR: "
              << "Txp to Gene Map not found for "
              << transcripts.size() - txpToGeneMap.size()
              <<" transcripts. Exiting" << std::flush;
    exit(1);
  }

  return geneCount;
}


template <typename ProtocolT>
bool CollapsedCellOptimizer::optimize(SCExpT& experiment,
                                      AlevinOpts<ProtocolT>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode){
  double relDiffTolerance{0.01};
  uint32_t maxIter {10000};

  auto& fullEqMap = experiment.equivalenceClassBuilder().eqMap();
  size_t numCells = trueBarcodes.size();
  size_t numWorkerThreads{1};

  if (aopt.numThreads > 1) {
    numWorkerThreads = aopt.numThreads - 1;
  }

  //get the keys of the map
  std::deque<std::pair<TranscriptGroup, uint32_t>> orderedTgroup;
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

  spp::sparse_hash_map<uint32_t, uint32_t> txpToGeneMap;
  spp::sparse_hash_map<std::string, uint32_t> geneIdxMap;

  uint32_t numGenes = getTxpToGeneMap(txpToGeneMap,
                                      experiment.transcripts(),
                                      aopt.geneMapFile.string(),
                                      geneIdxMap);

  if (aopt.dumpBarcodeEq){
    std::ofstream oFile;
    boost::filesystem::path oFilePath = aopt.outputDirectory / "cell_eq_order.txt";
    oFile.open(oFilePath.string());
    for (auto& bc : trueBarcodes) {
      oFile << bc << "\n";
    }
    oFile.close();

    {//dump transcripts names
      boost::filesystem::path tFilePath = aopt.outputDirectory / "transcripts.txt";
      std::ofstream tFile(tFilePath.string());
      for (auto& txp: experiment.transcripts()) {
        tFile << txp.RefName << "\n";
      }
      tFile.close();
    }
  }

  std::vector<CellState> skippedCB (numCells);
  std::atomic<uint32_t> bcount{0};
  std::atomic<uint64_t> totalDedupCounts{0};

  std::vector<std::thread> workerThreads;
  for (size_t tn = 0; tn < numWorkerThreads; ++tn) {
    workerThreads.emplace_back(optimizeCell,
                               std::ref(experiment),
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
                               aopt.doEM,
                               aopt.quiet,
                               std::ref(totalDedupCounts),
                               std::ref(txpToGeneMap),
                               numGenes,
                               aopt.txpLevel,
                               aopt.naive,
                               aopt.eqClassLevel,
                               aopt.debug);
  }

  for (auto& t : workerThreads) {
    t.join();
  }
  aopt.jointLog->info("Total {} UMI after deduplicating.",
                      totalDedupCounts);

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
        else if ( not aopt.debug ){
          std::cout<< "Skipped Barcodes are from High Confidence Region\n"
                   << " Should not happen"<<std::flush;
          exit(1);
        }
      }
    }
    numCells = trueBarcodes.size();
  }

  gzw.close_all_streams();

  bool txpLevel = aopt.txpLevel;
  if (not txpLevel){//dump gene names
    std::vector<std::string> geneNames(numGenes);
    for (auto geneIdx : geneIdxMap) {
      geneNames[geneIdx.second] = geneIdx.first;
    }
    boost::filesystem::path gFilePath = aopt.outputDirectory / "quants_mat_cols.txt";
    std::ofstream gFile(gFilePath.string());
    std::ostream_iterator<std::string> giterator(gFile, "\n");
    std::copy(geneNames.begin(), geneNames.end(), giterator);
    gFile.close();
  }


  std::vector<std::vector<double>> countMatrix;

  bool hasWhitelist = boost::filesystem::exists(aopt.whitelistFile);
  if(not aopt.nobarcode){
    if(not hasWhitelist  or aopt.dumpCsvCounts){
      aopt.jointLog->info("Clearing EqMap; Might take some time.");
      fullEqMap.clear();
      size_t numElem = txpLevel ? experiment.transcripts().size() : numGenes;
      std::string mode = txpLevel ? "transcript" : "gene";

      aopt.jointLog->info("Starting Import of the {} count matrix.", mode);
      countMatrix.resize(trueBarcodes.size(),
                         std::vector<double> (numElem, 0.0));
      auto zerod_cells = alevin::whitelist::populate_count_matrix(aopt.outputDirectory,
                                                                  aopt.debug,
                                                                  numElem,
                                                                  countMatrix);
      if (zerod_cells > 0) {
        aopt.jointLog->warn("Found {} cells with no reads,"
                            " ignoring due to debug mode.", zerod_cells);
      }

      aopt.jointLog->info("Done Importing {} count matrix for dimension {}x{}",
                          mode, numCells, numGenes);

      if (aopt.dumpCsvCounts){
        aopt.jointLog->info("Starting dumping cell v {} counts in csv format",
                            mode);
        std::ofstream qFile;
        boost::filesystem::path qFilePath = aopt.outputDirectory / "quants_mat.csv";
        qFile.open(qFilePath.string());
        for (auto& row : countMatrix) {
          for (auto cell : row) {
            qFile << cell << ',';
          }
          qFile << "\n";
        }
        qFile.close();

        aopt.jointLog->info("Finished dumping csv counts");
      }

      if(not hasWhitelist and not txpLevel){
        aopt.jointLog->info("Starting white listing");
        bool whitelistingSuccess = alevin::whitelist::performWhitelisting(aopt,
                                                                          umiCount,
                                                                          countMatrix,
                                                                          trueBarcodes,
                                                                          freqCounter,
                                                                          geneIdxMap,
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
      }
      else if(txpLevel){
        aopt.jointLog->warn("can't perform whitelisting on txp level cell count matrix");
      }
    }
  } // end-if no barcode

  return true;
} //end-optimize


namespace apt = alevin::protocols;
template
bool CollapsedCellOptimizer::optimize(SCExpT& experiment,
                                      AlevinOpts<apt::DropSeq>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode);
template
bool CollapsedCellOptimizer::optimize(SCExpT& experiment,
                                      AlevinOpts<apt::InDrop>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode);
template
bool CollapsedCellOptimizer::optimize(SCExpT& experiment,
                                      AlevinOpts<apt::Chromium>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode);
template
bool CollapsedCellOptimizer::optimize(SCExpT& experiment,
                                      AlevinOpts<apt::Gemcode>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode);
template
bool CollapsedCellOptimizer::optimize(SCExpT& experiment,
                                      AlevinOpts<apt::Celseq>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode);
template
bool CollapsedCellOptimizer::optimize(SCExpT& experiment,
                                      AlevinOpts<apt::Custom>& aopt,
                                      GZipWriter& gzw,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount,
                                      CFreqMapT& freqCounter,
                                      size_t numLowConfidentBarcode);
