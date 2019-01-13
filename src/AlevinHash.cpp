#include "AlevinHash.hpp"

void loadIndexGetTranscripts(bfs::path& indexDirectory,
                             std::vector<Transcript>& transcripts,
                             std::shared_ptr<spdlog::logger>& jointLog) {
  // ==== Figure out the index type
  boost::filesystem::path versionPath = indexDirectory / "versionInfo.json";
  SalmonIndexVersionInfo versionInfo;
  versionInfo.load(versionPath);
  if (versionInfo.indexVersion() == 0) {
    jointLog->error("Error: The index version file  {}"
                    " doesn't seem to exist.  "
                    "Please try re-building the salmon index.",
                    versionPath.string());
    jointLog->flush();
    exit(1);
  }

  // Check index version compatibility here
  auto indexType = versionInfo.indexType();
  // ==== Figure out the index type

  SalmonIndex salmonIndex(jointLog, indexType);
  salmonIndex.load(indexDirectory);

  // Now we'll have either an FMD-based index or a QUASI index
  // dispatch on the correct type.
  switch (salmonIndex.indexType()) {
  case SalmonIndexType::QUASI:
    if (salmonIndex.is64BitQuasi()) {
      jointLog->error("Error: This version of salmon does not support "
                      "the index mode.");
      jointLog->flush();
      exit(1);

      //if (salmonIndex.isPerfectHashQuasi()) {
      //  loadTranscriptsFromQuasi(salmonIndex_->quasiIndexPerfectHash64(),
      //                           sopt);
      //} else {
      //  loadTranscriptsFromQuasi(salmonIndex_->quasiIndex64(), sopt);
      //}
    } else {
      if (salmonIndex.isPerfectHashQuasi()) {
        auto* index = salmonIndex.quasiIndexPerfectHash32();
        size_t numRecords = index->txpNames.size();
        jointLog->info("Index contained {:n} targets", numRecords);

        double alpha = 0.005;
        for (auto i : boost::irange(size_t(0), numRecords)) {
          uint32_t id = i;
          const char* name = index->txpNames[i].c_str();
          uint32_t len = index->txpLens[i];
          // copy over the length, then we're done.
          transcripts.emplace_back(id, name, len, alpha);
        }
        //loadTranscriptsFromQuasi(salmonIndex_->quasiIndexPerfectHash32(),
        //                         sopt);
      } else {
        auto* index = salmonIndex.quasiIndex32();
        size_t numRecords = index->txpNames.size();
        jointLog->info("Index contained {:n} targets", numRecords);

        double alpha = 0.005;
        for (auto i : boost::irange(size_t(0), numRecords)) {
          uint32_t id = i;
          const char* name = index->txpNames[i].c_str();
          uint32_t len = index->txpLens[i];
          // copy over the length, then we're done.
          transcripts.emplace_back(id, name, len, alpha);
        }
      //  loadTranscriptsFromQuasi(salmonIndex_->quasiIndex32(), sopt);
      }
    }
    break;
  case SalmonIndexType::FMD:
    jointLog->error("Error: This version of salmon does not support "
                   "the FMD index mode.");
    jointLog->flush();
    exit(1);
  }
}

template <typename ProtocolT>
int salmonHashQuantify(AlevinOpts<ProtocolT>& aopt,
                       bfs::path& indexDirectory,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter) {
  aopt.jointLog->info("Reading BFH");
  aopt.jointLog->flush();

  bfs::path eqFilePath = aopt.bfhFile;
  std::ifstream equivFile(eqFilePath.string());

  size_t numTxps, numBcs, numEqclasses;

  // Number of transcripts
  equivFile >> numTxps;

  // Number of barcodes
  equivFile >> numBcs;

  // Number of equivalence classes
  equivFile >> numEqclasses;

  std::vector<std::string> txpNames (numTxps);
  for (size_t i=0; i<numTxps; i++) {
    equivFile >> txpNames[i] ;
  }

  size_t bcLength {aopt.protocol.barcodeLength};
  std::vector<std::string> bcNames (numBcs);
  for (size_t i=0; i<numBcs; i++) {
    equivFile >> bcNames[i] ;
    if (bcNames[i].size() != bcLength) {
      aopt.jointLog->error("CB {} has wrong length", bcNames[i]);
      aopt.jointLog->flush();
      exit(1);
    }
  }

  EqMapT countMap;
  countMap.set_max_resize_threads(1);
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

    uint32_t countValidator {0};
    for (size_t j=0; j<bgroupSize; j++){
      uint32_t bc;
      size_t ugroupSize;

      equivFile >> bc >> ugroupSize;
      std::string bcName = bcNames[bc];

      for (size_t k=0; k<ugroupSize; k++){
        std::string umiSeq;
        uint64_t umiIndex;
        uint32_t umiCount;
        equivFile >> umiSeq >> umiCount;

        bool isUmiIdxOk = umiObj.fromChars(umiSeq);
        if(isUmiIdxOk){
          umiIndex = umiObj.word(0);
        } else {
          aopt.jointLog->error("Umi Alevin Object conversion error");
          aopt.jointLog->flush();
          exit(1);
        }

        auto upfn = [bc, umiIndex, umiCount](SCTGValue& x) -> void {
          // update the count
          x.count += umiCount;
          x.updateBarcodeGroup(bc, umiIndex, umiCount);
        };

        countValidator += umiCount;
        SCTGValue value(bc, umiIndex, umiCount);
        countMap.upsert(txGroup, upfn, value);

        freqCounter[bcName] += umiCount;
      }// end-ugroup for
    }//end-bgroup for

    if (count != countValidator){
      aopt.jointLog->error("BFH eqclass count mismatch"
                           "{} Orignial, validator {} "
                           "Eqclass number {}",
                           count, countValidator, i);
      aopt.jointLog->flush();
      exit(1);
    }

    double completionFrac = i*100.0/numEqclasses;
    uint32_t percentCompletion {static_cast<uint32_t>(completionFrac)};
    if ( percentCompletion % 10 == 0 || percentCompletion > 95) {
      fmt::print(stderr, "\r{}Done Reading : {}{}%{}",
                 green, red, percentCompletion, RESET_COLOR);
    }
  }//end-eqclass for
  std::cerr<<std::endl;
  equivFile.close();


  std::vector<Transcript> transcripts;
  loadIndexGetTranscripts(indexDirectory,
                          transcripts,
                          aopt.jointLog);
  GZipWriter gzw(outputDirectory, aopt.jointLog);

  if(boost::filesystem::exists(aopt.whitelistFile)){
    aopt.jointLog->warn("can't use whitelist file in bfh Mode"
                        ";Ignroing the file");
    aopt.jointLog->flush();
  }

  alevinOptimize(bcNames, transcripts, countMap,
                 aopt, gzw, freqCounter, 0);
  return 0;
}


template
int salmonHashQuantify(AlevinOpts<apt::Chromium>& aopt,
                       bfs::path& indexDirectory,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::ChromiumV3>& aopt,
                       bfs::path& indexDirectory,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::Gemcode>& aopt,
                       bfs::path& indexDirectory,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::DropSeq>& aopt,
                       bfs::path& indexDirectory,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::InDrop>& aopt,
                       bfs::path& indexDirectory,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::CELSeq>& aopt,
                       bfs::path& indexDirectory,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::CELSeq2>& aopt,
                       bfs::path& indexDirectory,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::Custom>& aopt,
                       bfs::path& indexDirectory,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
