#include "AlevinHash.hpp"

void getTxpToGeneMap(spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                     std::vector<std::string>& transcripts,
                     const std::string& geneMapFile,
                     spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap){
  std::ifstream t2gFile(geneMapFile);

  spp::sparse_hash_map<std::string, uint32_t> txpIdxMap(transcripts.size());

  for (size_t i=0; i<transcripts.size(); i++){
    txpIdxMap[ transcripts[i] ] = i;
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
}

size_t readBfh(bfs::path& eqFilePath,
               std::vector<std::string>& txpNames,
               size_t bcLength,
               EqMapT &countMap,
               std::vector<std::string>& bcNames,
               CFreqMapT& freqCounter,
               TrueBcsT& trueBarcodes,
               size_t& totalNormalized,
               bool hasWhitelist
               ) {
  if (hasWhitelist && trueBarcodes.size() == 0) {
    fmt::print(stderr, "whitelist file had 0 CB");
    std::cerr<<std::endl<<std::flush;;
    return 0;
  }

  std::ifstream equivFile(eqFilePath.string());

  size_t numReads{0};
  size_t numTxps, numBcs, numEqclasses;

  // Number of transcripts
  equivFile >> numTxps;

  // Number of barcodes
  equivFile >> numBcs;

  // Number of equivalence classes
  equivFile >> numEqclasses;

  txpNames.resize(numTxps);
  for (size_t i=0; i<numTxps; i++) {
    equivFile >> txpNames[i] ;
  }

  bcNames.resize(numBcs);
  for (size_t i=0; i<numBcs; i++) {
    equivFile >> bcNames[i] ;

    if (bcNames[i].size() != bcLength) {
      fmt::print(stderr, "CB {} has wrong length", bcNames[i]);
      std::cerr<<std::endl<<std::flush;;
      return 0;
    }
  }

  spp::sparse_hash_map<size_t, size_t> barcodeMap;
  { // bcname and bc index rearrangement based on external whitelist
    if (not hasWhitelist) {
      for (size_t i=0; i<bcNames.size(); i++) {
        barcodeMap[i] = i;
      }
    } else {
      // convert set to indexed vector
      size_t idx{0};
      spp::sparse_hash_map<std::string, size_t> trueBarcodeMap;
      for (auto& bc: trueBarcodes) {
        trueBarcodeMap[bc] = idx;
        idx += 1;
      }

      // extracting relevant barcodes
      for (size_t i=0; i<bcNames.size(); i++) {
        if (trueBarcodeMap.contains(bcNames[i])) {
          barcodeMap[i] = trueBarcodeMap[ bcNames[i] ];
        }
      }

      bcNames.clear();
      bcNames.resize(trueBarcodeMap.size());
      for(auto& it: trueBarcodeMap) {
        bcNames[it.second] = it.first;
      }

    } // end else case of not hasWhitelist
  } // end name/index rearrangement

  countMap.max_num_worker_threads(1);
  countMap.reserve(1000000);

  alevin::types::AlevinUMIKmer umiObj;
  //printing on screen progress
  const char RESET_COLOR[] = "\x1b[0m";
  char green[] = "\x1b[30m";
  green[3] = '0' + static_cast<char>(fmt::GREEN);
  char red[] = "\x1b[30m";
  red[3] = '0' + static_cast<char>(fmt::RED);
  std::cerr<<std::endl;

  for (size_t i=0; i<numEqclasses; i++) {
    uint32_t count;
    size_t labelSize ;
    equivFile >> labelSize;

    std::vector<uint32_t> txps(labelSize);
    for (auto& tid : txps) { equivFile >> tid; }
    auto txGroup = TranscriptGroup (txps);

    size_t bgroupSize;
    equivFile >> count >> bgroupSize;

    size_t normalizer{0};
    uint32_t countValidator {0};
    for (size_t j=0; j<bgroupSize; j++){
      size_t ugroupSize;
      std::string bcName;
      bool skipCB {false};
      uint32_t old_bc, new_bc;

      equivFile >> old_bc >> ugroupSize;
      if (not barcodeMap.contains(old_bc)) {
        skipCB = true;
      } else {
        new_bc = barcodeMap[old_bc];
        bcName = bcNames[new_bc];
      }

      for (size_t k=0; k<ugroupSize; k++){
        std::string umiSeq;
        uint32_t umiCount;

        equivFile >> umiSeq >> umiCount;
        countValidator += umiCount;
        if (skipCB) {normalizer += umiCount; continue;}

        uint64_t umiIndex;
        bool isUmiIdxOk = umiObj.fromChars(umiSeq);
        if(isUmiIdxOk){
          umiIndex = umiObj.word(0);
          auto upfn = [new_bc, umiIndex, umiCount](SCTGValue& x) -> void {
            // update the count
            x.count += umiCount;
            x.updateBarcodeGroup(new_bc, umiIndex, umiCount);
          };

          SCTGValue value(umiCount, new_bc, umiIndex, true);
          countMap.upsert(txGroup, upfn, value);
          freqCounter[bcName] += umiCount;
        }
      }// end-ugroup for
    }//end-bgroup for

    if (count != countValidator){
      fmt::print(stderr, "BFH eqclass count mismatch"
                 "{} Orignial, validator {} "
                 "Eqclass number {}",
                 count, countValidator, i);
      std::cerr<<std::endl<<std::flush;;
      return 0;
    }

    totalNormalized += normalizer;
    numReads += (countValidator - normalizer);
    double completionFrac = i*100.0/numEqclasses;
    uint32_t percentCompletion {static_cast<uint32_t>(completionFrac)};
    if ( percentCompletion % 10 == 0 || percentCompletion > 95) {
      fmt::print(stderr, "\r{}Done Reading : {}{}%{} and skipped reads: {}{}{}",
                 green, red, percentCompletion, green,
                 red, normalizer, RESET_COLOR);
    }
  }//end-eqclass for
  std::cerr<<std::endl;
  equivFile.close();

  return numReads;
}

template <typename ProtocolT>
int salmonHashQuantify(AlevinOpts<ProtocolT>& aopt,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter) {
  TrueBcsT trueBarcodes;
  bool hasWhitelist = boost::filesystem::exists(aopt.whitelistFile);
  if( hasWhitelist ){
    alevin::utils::readWhitelist(aopt.whitelistFile,
                                 trueBarcodes);

    aopt.jointLog->info("Done importing white-list Barcodes");
    aopt.jointLog->info("Total {} white-listed Barcodes", trueBarcodes.size());
  }

  EqMapT countMap;
  size_t numReads {0};
  size_t totalNormalized {0};
  std::vector<std::string> txpNames, bcNames;
  { // Populating Bfh
    aopt.jointLog->info("Reading BFH");
    aopt.jointLog->flush();

    numReads = readBfh( aopt.bfhFile,
                        txpNames,
                        aopt.protocol.barcodeLength,
                        countMap, bcNames,
                        freqCounter,
                        trueBarcodes,
                        totalNormalized,
                        hasWhitelist );
    if( numReads == 0 ){
      aopt.jointLog->error("error reading bfh");
      aopt.jointLog->flush();
      std::exit(74);
    }

    if (totalNormalized > 0) {
      aopt.jointLog->warn("Skipped {} reads due to external whitelist",
                          totalNormalized);
    }
    aopt.jointLog->info("Fount total {} reads in bfh",
                        numReads);
    aopt.jointLog->flush();
  } // Done populating Bfh

  // extracting meta data for calling alevinOptimize
  aopt.jointLog->info("Reading transcript to gene Map");
  spp::sparse_hash_map<uint32_t, uint32_t> txpToGeneMap;
  spp::sparse_hash_map<std::string, uint32_t> geneIdxMap;
  getTxpToGeneMap(txpToGeneMap,
                  txpNames,
                  aopt.geneMapFile.string(),
                  geneIdxMap);

  GZipWriter gzw(outputDirectory, aopt.jointLog);
  alevinOptimize(bcNames, txpToGeneMap, geneIdxMap,
                 countMap, aopt, gzw, freqCounter, 0);
  return 0;
}


template
int salmonHashQuantify(AlevinOpts<apt::Chromium>& aopt,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::ChromiumV3>& aopt,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::Gemcode>& aopt,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::DropSeq>& aopt,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::CITESeq>& aopt,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::InDrop>& aopt,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::CELSeq>& aopt,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::CELSeq2>& aopt,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::QuartzSeq2>& aopt,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::Custom>& aopt,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);

template
int salmonHashQuantify(AlevinOpts<apt::CustomGeometry>& aopt,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);

