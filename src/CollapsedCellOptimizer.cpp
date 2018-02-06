#include "CollapsedCellOptimizer.hpp"

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
                  std::shared_ptr<spdlog::logger>& jointlog,
                  bfs::path& outDirPath,
                  std::unordered_set<uint32_t>& activeTranscriptIds) {

  // An EM termination criterion, adopted from Bray et al. 2016
  uint32_t minIter {50};
  double relDiffTolerance {0.01};
  uint32_t maxIter {10000};
  size_t numClasses = txpGroups.size();
  double uniformTxpWeight = 1.0 / activeTranscriptIds.size();

  std::vector<uint64_t> eqCounts(numClasses, 0);
  CollapsedCellOptimizer::SerialVecType alphas(transcripts.size(), 0.0);
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
    return false;
  }

  GZipWriter gzw(outDirPath, jointlog);
  gzw.writeAbundances(alphas, const_cast<std::vector<Transcript>&>(transcripts));

  return true;
}


void optimizeCell(SCExpT& experiment,
                  std::vector<std::string>& trueBarcodes,
                  std::atomic<uint32_t>& barcode,
                  size_t totalCells,
                  eqMapT& eqMap,
                  std::deque<TranscriptGroup>& orderedTgroup,
                  std::shared_ptr<spdlog::logger>& jointlog,
                  bfs::path& outDir, std::vector<uint32_t>& umiCount,
                  tbb::atomic<uint32_t>& skippedCBcount,
                  bool verbose, size_t umiLength, bool noEM,
                  spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap){
  size_t numCells {trueBarcodes.size()};
  size_t trueBarcodeIdx;

  // looping over until the end of the file
  while((trueBarcodeIdx = barcode++) < totalCells) {
    // per-cell level optimization
    if(umiCount[trueBarcodeIdx] == 0){
      //skip the barcode if no mapped UMI
      skippedCBcount += 1;
      continue;
    }
    auto& trueBarcodeStr = trueBarcodes[trueBarcodeIdx];
    //creating quant directory for each cell and file

    bfs::path qDirPath = outDir / "cell" / trueBarcodeStr ;
    if (!bfs::exists(qDirPath)) {
      bool dirSuccess = boost::filesystem::create_directories(qDirPath);
      if (!dirSuccess) {
        fmt::print(stderr,"\nCould not create output directory {}\nExiting Now.",
                   qDirPath.string());
        exit(1);
      }
    }


    //extracting per-cell level eq class information
    const std::vector<Transcript>& transcripts = experiment.transcripts();
    std::unordered_set<uint32_t> activetranscriptids;
    std::vector<tgrouplabelt> txpgroups;
    std::vector<tgroupweightvec> txpgroupcombinedweights;
    std::vector<uint64_t> origcounts;
    uint64_t totalcount{0};

    std::unordered_map<uint32_t, std::unordered_set<uint64_t>> umiBiasList;
    std::ofstream qFile;

    if(verbose){
      bfs::path qFilePath = qDirPath / "cell_eq_classes.txt";
      qFile.open(qFilePath.string());

      // Number of transcripts
      qFile << transcripts.size() << '\n';

      // Number of equivalence classes
      //qFile << "XX" << '\n';

      // dump transcript names
      for (auto& t : transcripts) {
        qFile << t.RefName << '\n';
      }
    }


    for (auto& key : orderedTgroup) {
      //traversing each class and copying relevant data.
      bool isKeyPresent = eqMap.find_fn(key, [&](const SCTGValue& val){
          auto& bg = val.barcodeGroup;
          auto bcIt = bg.find(trueBarcodeIdx);

          if (bcIt != bg.end()){
            // sub-selecting bgroup of this barcode only
            auto eqCount = dedupReads(umiLength,
                                      jointlog,
                                      bcIt->second,
                                      umiBiasList,
                                      key);

            if ( eqCount != 0 ) {
              const std::vector<uint32_t>& txps = key.txps;
              // Avi -> Major Hack
              // Basically probStartpos is 1/effec_length.
              // effec_length for all txp is 100
              // Eqclass weights are set to 1 since sopt.noRichEqClasses is true
              // Making aux as 1/#oftxps_in_class
              // Ideally should be taken from eqclass.combinedWeights but
              // ignoring for now
              const tgroupweightvec auxs(txps.size(), 1.0/txps.size());
              // convert to non-atomic
              txpgroupcombinedweights.emplace_back(auxs.begin(), auxs.end());
              txpgroups.push_back(txps);
              origcounts.emplace_back(eqCount);

              totalcount += eqCount;

              if(verbose){
                // group size
                qFile << txps.size();
                // dump txp-group members
                for (auto tid : txps) { qFile << "\t" << tid; }
                qFile << "\t" << eqCount << "\n";
              }

              // currently add only 1 length eqclass to active txps
              if (txps.size() == 1) {
                for (auto& t : txps) {
                  activetranscriptids.insert(t);
                }
              }
            }
          }
        });

      if(!isKeyPresent){
        jointlog->error("Not able to find key in Cuckoo hash map."
                        "Please Report this issue on github");
        exit(1);
      }
    }

    // get the list of active gene ids
    std::unordered_set<uint32_t> activeGeneIds;
    for(auto& tid: activetranscriptids){
      if(txpToGeneMap.contains(tid)){
        uint32_t gid = txpToGeneMap[tid];
        activeGeneIds.insert(gid);
      }
      else{
        std::cerr << "Out of Range error for txp to gene Map in 1st: " << '\n';
        std::cerr << tid << "\t" << transcripts[tid].RefName << " not found";
        exit(1);
      }
    }

    // parse through eqclass once more to filter out non-relevant txps
    for (auto& txps: txpgroups){
      if (txps.size() != 1){
        for (auto& tid: txps){
          uint32_t gid;
          if(txpToGeneMap.contains(tid)){
            gid = txpToGeneMap.at(tid);
          }
          else{
            std::cerr << "Out of Range error for txp to gene Map: " << '\n';
            std::cerr << tid << "\t not found";
            exit(1);
          }
          if (activeGeneIds.find(gid) != activeGeneIds.end()){
            activetranscriptids.insert(tid);
          }
        }
      }
    }

    if(verbose){
      qFile.close();
      //jointlog->info("optimizing over {} equivalence classes", txpgroups.size());
    }

    double totalnumfrags{static_cast<double>(totalcount)};

    if (activetranscriptids.size() == 0) {
      jointlog->error("it seems that no transcripts are expressed; something is "
                      "likely wrong!");
      std::exit(1);
    }

    if(not noEM){
      bool isEMok = runPerCellEM(txpgroups,
                                 txpgroupcombinedweights,
                                 origcounts,
                                 transcripts,
                                 totalnumfrags,
                                 jointlog,
                                 qDirPath,
                                 activetranscriptids);
      if( !isEMok ){
        jointlog->error("EM iteration for cell {} failed \n"
                        "Please Report this on github.", trueBarcodeStr);
        std::exit(1);
      }
    }

    const char RESET_COLOR[] = "\x1b[0m";
    char green[] = "\x1b[30m";
    green[3] = '0' + static_cast<char>(fmt::GREEN);
    char red[] = "\x1b[30m";
    red[3] = '0' + static_cast<char>(fmt::RED);

    double cellCount {barcode};//numCells-jqueue.size_approx()};
    if (cellCount > totalCells) { cellCount = totalCells; }
    double percentCompletion {cellCount*100/numCells};
    fmt::print(stderr, "\033[A\r\r{}Analyzed {} cells ({}{}%{} of all).{}\n",
               green, cellCount, red, round(percentCompletion), green, RESET_COLOR);

    //found = jqueue.try_dequeue(trueBarcodeIdx);
  }
}

void getTxpToGeneMap(spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                     const std::vector<Transcript>& transcripts,
                     const std::string& geneMapFile){
  std::string fname = geneMapFile;
  std::ifstream t2gFile(fname);

  spp::sparse_hash_map<std::string, uint32_t> txpIdxMap(transcripts.size());
  spp::sparse_hash_map<std::string, uint32_t> geneIdxMap;

  for (size_t i=0; i<transcripts.size(); i++){
    txpIdxMap[ transcripts[i].RefName ] = i;
  }

  uint32_t tid, gid, geneCount{0};
  std::string tStr, gStr;
  if(t2gFile.is_open()) {
    while( not t2gFile.eof() ) {
      t2gFile >> tStr >> gStr;
      if (geneIdxMap.contains(gStr)){
        gid = geneIdxMap[gStr];
      }
      else{
        gid = geneCount;
        geneIdxMap[gStr] = gid;
        geneCount++;
      }

      if(not txpIdxMap.contains(tStr)){
        continue;
      }
      tid = txpIdxMap[tStr];

      txpToGeneMap[tid] = gid;
    }
    t2gFile.close();
  }
}


template <typename ProtocolT>
bool CollapsedCellOptimizer::optimize(SCExpT& experiment,
                                      AlevinOpts<ProtocolT>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount){
  double relDiffTolerance{0.01};
  uint32_t maxIter {10000};

  auto& fullEqMap = experiment.equivalenceClassBuilder().eqMap();
  size_t numCells = trueBarcodes.size();
  size_t numWorkerThreads{1};

  if (aopt.numThreads > 1) {
    numWorkerThreads = aopt.numThreads - 1;
  }

  //get the keys of the map
  std::deque<TranscriptGroup> orderedTgroup;
  for(const auto& kv : fullEqMap.lock_table()){
    if(kv.first.txps.size() == 1){
      orderedTgroup.push_front(kv.first);
    }
    else{
      orderedTgroup.push_back(kv.first);
    }
  }

  spp::sparse_hash_map<uint32_t, uint32_t> txpToGeneMap;
  getTxpToGeneMap(txpToGeneMap,
                  experiment.transcripts(),
                  aopt.geneMapFile.string());

  tbb::atomic<uint32_t> skippedCBcount{0};
  std::atomic<uint32_t> bcount{0};

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
                               std::ref(skippedCBcount),
                               aopt.dumpBarcodeEq,
                               aopt.protocol.umiLength,
                               aopt.noEM,
                               std::ref(txpToGeneMap));
  }

  for (auto& t : workerThreads) {
    t.join();
  }
  if(skippedCBcount>0){
    aopt.jointLog->warn("Skipped {} barcodes due to No mapped read",
                        skippedCBcount);
  }

  return true;
} //end-optimize


namespace apt = alevin::protocols;
template
bool CollapsedCellOptimizer::optimize(SCExpT& experiment,
                                      AlevinOpts<apt::DropSeq>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount);
template
bool CollapsedCellOptimizer::optimize(SCExpT& experiment,
                                      AlevinOpts<apt::InDrop>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount);
template
bool CollapsedCellOptimizer::optimize(SCExpT& experiment,
                                      AlevinOpts<apt::Chromium>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount);
template
bool CollapsedCellOptimizer::optimize(SCExpT& experiment,
                                      AlevinOpts<apt::Custom>& aopt,
                                      std::vector<std::string>& trueBarcodes,
                                      std::vector<uint32_t>& umiCount);
