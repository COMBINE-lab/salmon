template <typename EQBuilderT>
ReadExperiment<EQBuilderT>::ReadExperiment(
    std::vector<ReadLibrary>& readLibraries, SalmonIndex* salmonIndex,
    SalmonOpts& sopt)
    : state_(readLibraries, sopt.numConditionalGCBins, sopt.numFragGCBins),
      salmonIndex_(salmonIndex),
      transcripts_(std::vector<Transcript>()),
      totalAssignedFragments_(0),
      eqBuilder_(sopt.jointLog, sopt.maxHashResizeThreads) {
  for (auto& rl : state_.library) {
    rl.checkValid();
  }

  size_t maxFragLen = sopt.fragLenDistMax;
  double meanFragLen = sopt.fragLenDistPriorMean;
  double fragLenStd = sopt.fragLenDistPriorSD;
  size_t fragLenKernelN = 4;
  double fragLenKernelP = 0.5;
  fragLengthDist_.reset(
      new FragmentLengthDistribution(1.0, maxFragLen, meanFragLen, fragLenStd,
                                     fragLenKernelN, fragLenKernelP, 1));

  if (state_.library.front().getFormat().type == ReadType::SINGLE_END) {
    std::vector<double> logPMF;
    size_t minVal;
    size_t maxVal;
    fragLengthDist_->dumpPMF(logPMF, minVal, maxVal);
    double sum = salmon::math::LOG_0;
    for (auto v : logPMF) {
      sum = salmon::math::logAdd(sum, v);
    }
    for (auto& v : logPMF) {
      v -= sum;
    }

    std::vector<double> pmf(maxVal + 1, 0.0);
    for (size_t i = minVal; i < maxVal; ++i) {
      pmf[i] = 100.0 * std::exp(logPMF[i - minVal]);
    }

    using distribution_utils::DistributionSpace;
    state_.conditionalMeans = distribution_utils::correctionFactorsFromMass(
        pmf, DistributionSpace::LINEAR);
  }

  fmt::MemoryWriter infostr;
  switch (salmonIndex_->indexType()) {
  case SalmonIndexType::QUASI:
    infostr << "Error: This version of salmon does not support a RapMap-based index.";
    throw std::invalid_argument(infostr.str());
  case SalmonIndexType::FMD:
    infostr << "Error: This version of salmon does not support the FMD index mode.";
    throw std::invalid_argument(infostr.str());
  case SalmonIndexType::PUFF:
    if (salmonIndex_->isSparse()) {
      loadTranscriptsFromPuff(salmonIndex_->puffSparseIndex(), sopt);
    } else {
      loadTranscriptsFromPuff(salmonIndex_->puffIndex(), sopt);
    }
  }

  salmon::utils::markAuxiliaryTargets(sopt, transcripts_);
  clusters_.reset(new ClusterForest(transcripts_.size(), transcripts_));
}

template <typename EQBuilderT>
void ReadExperiment<EQBuilderT>::updateTranscriptLengthsAtomic(
    std::atomic<bool>& done) {
  if (sl_.try_lock()) {
    if (!done) {
      auto fld = fragLengthDist_.get();
      std::vector<double> logPMF;
      size_t minVal;
      size_t maxVal;
      fld->dumpPMF(logPMF, minVal, maxVal);
      double sum = salmon::math::LOG_0;
      for (auto v : logPMF) {
        sum = salmon::math::logAdd(sum, v);
      }
      for (auto& v : logPMF) {
        v -= sum;
      }

      std::vector<double> pmf(maxVal + 1, 0.0);
      for (size_t i = minVal; i < maxVal; ++i) {
        pmf[i] = 100.0 * std::exp(logPMF[i - minVal]);
      }

      using distribution_utils::DistributionSpace;
      auto correctionFactors = distribution_utils::correctionFactorsFromMass(
          pmf, DistributionSpace::LINEAR);
      distribution_utils::computeSmoothedEffectiveLengths(
          pmf.size(), transcripts_, correctionFactors, DistributionSpace::LOG);
      done = true;
    }
    sl_.unlock();
  }
}

template <typename EQBuilderT>
template <typename PuffIndexT>
void ReadExperiment<EQBuilderT>::loadTranscriptsFromPuff(
    PuffIndexT* idx_, const SalmonOpts& sopt) {
  const auto& refNames = idx_->getFullRefNames();
  const auto& refLengths = idx_->getFullRefLengths();
  const auto& completeRefLengths = idx_->getFullRefLengthsComplete();

  size_t numRecords = refNames.size();
  auto log = sopt.jointLog.get();
  numDecoys_ = 0;

  log->info("Index contained {:n} targets", numRecords);
  transcripts_.reserve(numRecords);
  std::vector<uint32_t> lengths;
  lengths.reserve(numRecords);
  auto k = idx_->k();
  int64_t numShort{0};
  double alpha = 0.005;
  auto& allRefSeq = idx_->refseq_;
  auto& refAccumLengths = idx_->refAccumLengths_;
  for (auto i : boost::irange(size_t(0), numRecords)) {
    uint32_t id = i;
    bool isShort = refLengths[i] <= k;
    bool isDecoy = idx_->isDecoy(i - numShort);

    if (isDecoy and !sopt.validateMappings) {
      log->warn("The index contains decoy targets, but these should not be used in the "
                "absence of selective-alignment (--validateMappings, --mimicBT2 or --mimicStrictBT2). "
                "Skipping loading of decoys.");
      break;
    }

    const char* name = refNames[i].c_str();
    uint32_t len = refLengths[i];
    transcripts_.emplace_back(id, name, len, alpha);
    auto& txp = transcripts_.back();
    txp.setCompleteLength(completeRefLengths[i]);

    if (!isDecoy and !isShort and (sopt.biasCorrect or sopt.gcBiasCorrect)) {
      auto tid = i - numShort;
      int64_t refAccPos = tid > 0 ? refAccumLengths[tid - 1] : 0;
      int64_t refTotalLength = refAccumLengths[tid] - refAccPos;
      if (len != refTotalLength) {
        log->warn("len : {:n}, but txp.RefLength : {:n} :: refTotalLength : {:n}",
                  len, txp.RefLength, refTotalLength);
      }
      char* tseq =
          pufferfish::util::getRefSeqOwned(allRefSeq, refAccPos, refTotalLength);
      txp.setSequenceOwned(tseq, sopt.gcBiasCorrect, sopt.reduceGCMemory);
    }
    txp.setDecoy(isDecoy);
    numShort += isShort ? 1 : 0;
    if (isDecoy) {
      ++numDecoys_;
    } else {
      lengths.push_back(txp.RefLength);
    }
  }
  sopt.jointLog->info("Number of decoys : {:n}", numDecoys_);
  auto firstDecoyIndex = idx_->firstDecoyIndex();
  if (firstDecoyIndex < numRecords) {
    sopt.jointLog->info("First decoy index : {:n} ", firstDecoyIndex);
  }
  setTranscriptLengthClasses_(lengths, state_.posBiasFW.size());
}

template <typename EQBuilderT>
template <typename QuasiIndexT>
void ReadExperiment<EQBuilderT>::loadTranscriptsFromQuasi(
    QuasiIndexT* idx_, const SalmonOpts& sopt) {
  size_t numRecords = idx_->txpNames.size();
  auto log = sopt.jointLog.get();
  numDecoys_ = 0;

  log->info("Index contained {:n} targets", numRecords);
  transcripts_.reserve(numRecords);
  std::vector<uint32_t> lengths;
  lengths.reserve(numRecords);
  double alpha = 0.005;
  for (auto i : boost::irange(size_t(0), numRecords)) {
    uint32_t id = i;
    bool isDecoy = idx_->isDecoy(i);

    if (isDecoy and !sopt.validateMappings) {
      log->warn("The index contains decoy targets, but these should not be used in the "
                "absence of selective-alignment (--validateMappings, --mimicBT2 or --mimicStrictBT2). "
                "Skipping loading of decoys.");
      break;
    }

    const char* name = idx_->txpNames[i].c_str();
    uint32_t len = idx_->txpLens[i];
    transcripts_.emplace_back(id, name, len, alpha);
    auto& txp = transcripts_.back();
    txp.setCompleteLength(idx_->txpCompleteLens[i]);
    txp.setSequenceBorrowed(idx_->seq.c_str() + idx_->txpOffsets[i],
                            sopt.gcBiasCorrect, sopt.reduceGCMemory);
    txp.setDecoy(isDecoy);
    if (isDecoy) {
      ++numDecoys_;
    }

    lengths.push_back(txp.RefLength);
  }
  setTranscriptLengthClasses_(lengths, state_.posBiasFW.size());
}

template <typename EQBuilderT>
std::string ReadExperiment<EQBuilderT>::readFilesAsString() {
  std::stringstream sstr;
  size_t ln{0};
  size_t numReadLibraries{state_.library.size()};

  for (auto& rl : state_.library) {
    sstr << rl.readFilesAsString();
    if (ln++ < numReadLibraries) {
      sstr << "; ";
    }
  }
  return sstr.str();
}

template <typename EQBuilderT>
void ReadExperiment<EQBuilderT>::summarizeLibraryTypeCounts(
    boost::filesystem::path& opath) {
  LibraryFormat fmt1(ReadType::SINGLE_END, ReadOrientation::NONE,
                     ReadStrandedness::U);
  LibraryFormat fmt2(ReadType::SINGLE_END, ReadOrientation::NONE,
                     ReadStrandedness::U);

  std::ofstream os(opath.string());
  cereal::JSONOutputArchive oa(os);

  fmt::MemoryWriter errstr;
  auto log = spdlog::get("jointLog");

  uint64_t numFmt1{0};
  uint64_t numFmt2{0};
  uint64_t numAgree{0};
  uint64_t numDisagree{0};
  uint64_t numDisagreeUnstranded{0};
  uint64_t numDisagreeStranded{0};

  for (auto& rl : state_.library) {
    auto fmt = rl.format();
    auto& counts = rl.libTypeCounts();

    oa(cereal::make_nvp("read_files", rl.readFilesAsString()));
    std::string expectedFormat = rl.format().toString();
    oa(cereal::make_nvp("expected_format", expectedFormat));

    double compatFragmentRatio =
        rl.numCompat() / static_cast<double>(numAssignedFragments_);
    oa(cereal::make_nvp("compatible_fragment_ratio", compatFragmentRatio));
    oa(cereal::make_nvp("num_compatible_fragments", rl.numCompat()));
    oa(cereal::make_nvp("num_assigned_fragments", numAssignedFragments_.load()));

    std::vector<ReadStrandedness> strands;
    switch (fmt.orientation) {
    case ReadOrientation::SAME:
    case ReadOrientation::NONE:
      strands.push_back(ReadStrandedness::S);
      strands.push_back(ReadStrandedness::A);
      break;
    case ReadOrientation::AWAY:
    case ReadOrientation::TOWARD:
      strands.push_back(ReadStrandedness::SA);
      strands.push_back(ReadStrandedness::AS);
      break;
    }

    fmt1.type = fmt.type;
    fmt1.orientation = fmt.orientation;
    fmt1.strandedness = strands[0];
    fmt2.type = fmt.type;
    fmt2.orientation = fmt.orientation;
    fmt2.strandedness = strands[1];

    numFmt1 = 0;
    numFmt2 = 0;

    for (size_t i = 0; i < counts.size(); ++i) {
      if (i == fmt.formatID()) {
        numAgree = counts[i];
      } else {
        numDisagreeStranded += counts[i];
      }

      if (i == fmt1.formatID()) {
        numFmt1 = counts[i];
      } else if (i == fmt2.formatID()) {
        numFmt2 = counts[i];
      } else {
        numDisagreeUnstranded += counts[i];
      }
    }

    double ratio = 0.0;
    if (fmt.strandedness == ReadStrandedness::U) {
      numAgree = numFmt1 + numFmt2;
      numDisagree = numDisagreeUnstranded;
      ratio = numAgree > 0
                  ? (static_cast<double>(numFmt1) / (numFmt1 + numFmt2))
                  : 0.0;

      if (numAgree == 0) {
        errstr << "NOTE: Read Lib [" << rl.readFilesAsString() << "] :\n";
        errstr << "\nFound no concordant and consistent mappings. "
                  "If this is a paired-end library, are you sure the reads are properly paired? "
                  "check the file: "
               << opath.string() << " for details\n";
        log->warn(errstr.str());
        errstr.clear();
      } else if (std::abs(ratio - 0.5) > 0.01) {
        errstr << "NOTE: Read Lib [" << rl.readFilesAsString() << "] :\n";
        errstr << "\nDetected a *potential* strand bias > 1\% in an "
                  "unstranded protocol "
               << "check the file: " << opath.string() << " for details\n";
        log->warn(errstr.str());
        errstr.clear();
      }
    } else {
      numDisagree = numDisagreeStranded;
      ratio = numAgree > 0
                  ? (static_cast<double>(numFmt1) / (numFmt1 + numFmt2))
                  : 0.0;
    }

    oa(cereal::make_nvp("num_frags_with_concordant_consistent_mappings",
                        numAgree));
    oa(cereal::make_nvp("num_frags_with_inconsistent_or_orphan_mappings",
                        numDisagree));
    oa(cereal::make_nvp("strand_mapping_bias", ratio));

    double disagreeRatio = 1.0 - compatFragmentRatio;
    if (disagreeRatio > 0.05) {
      errstr << "NOTE: Read Lib [" << rl.readFilesAsString() << "] :\n";
      errstr << "\nGreater than 5\% of the fragments "
             << "disagreed with the provided library type; "
             << "check the file: " << opath.string() << " for details\n";

      log->warn(errstr.str());
      errstr.clear();
    }

    for (size_t i = 0; i < counts.size(); ++i) {
      std::string desc = LibraryFormat::formatFromID(i).toString();
      if (!desc.empty()) {
        oa(cereal::make_nvp(desc, counts[i].load()));
      }
    }
  }
}

template <typename EQBuilderT>
void ReadExperiment<EQBuilderT>::setTranscriptLengthClasses_(
    std::vector<uint32_t>& lengths, size_t nbins) {
  auto n = lengths.size();
  if (n > nbins) {
    lengthQuantiles_.clear();
    lengthQuantiles_.reserve(nbins);

    size_t step = lengths.size() / nbins;
    size_t cumStep = 0;
    for (size_t i = 0; i < nbins; ++i) {
      cumStep += step;
      size_t ind = std::min(cumStep, n - 1);
      std::nth_element(lengths.begin(), lengths.begin() + ind, lengths.end());
      lengthQuantiles_.push_back(*(lengths.begin() + ind));
    }
  } else {
    lengthQuantiles_.clear();
    lengthQuantiles_.reserve(n);
    std::sort(lengths.begin(), lengths.end());
    for (auto l : lengths) {
      lengthQuantiles_.push_back(l);
    }
    state_.posBiasFW.resize(n);
    state_.posBiasRC.resize(n);
    state_.posBiasExpectFW.resize(n);
    state_.posBiasExpectRC.resize(n);
  }

  auto qb = lengthQuantiles_.begin();
  auto qe = lengthQuantiles_.end();
  auto maxQuant = std::distance(qb, qe) - 1;
  for (auto& t : transcripts_) {
    auto ind = std::min(
        maxQuant, std::distance(qb, std::upper_bound(qb, qe, t.RefLength)));
    t.lengthClassIndex(ind);
  }
}
