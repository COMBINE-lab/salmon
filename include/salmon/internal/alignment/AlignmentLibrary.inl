template <typename FragT, typename EQBuilderT, typename AlignModelT>
AlignmentLibrary<FragT, EQBuilderT, AlignModelT>::AlignmentLibrary(
    std::vector<boost::filesystem::path>& alnFiles,
    const boost::filesystem::path& transcriptFile, LibraryFormat libFmt,
    SalmonOpts& salmonOpts)
    : alignmentFiles_(alnFiles),
      transcriptFile_(transcriptFile),
      state_(libFmt, salmonOpts.numConditionalGCBins, salmonOpts.numFragGCBins),
      transcripts_(std::vector<Transcript>()),
      eqBuilder_(salmonOpts.jointLog, salmonOpts.maxHashResizeThreads),
      quantificationPasses_(0) {
  namespace bfs = boost::filesystem;

  for (auto& alignmentFile : alignmentFiles_) {
    if (!bfs::exists(alignmentFile)) {
      std::stringstream ss;
      ss << "The provided alignment file: " << alignmentFile
         << " does not exist!\n";
      throw std::invalid_argument(ss.str());
    }
  }

  if (!bfs::exists(transcriptFile_)) {
    std::stringstream ss;
    ss << "The provided transcript file: " << transcriptFile_
       << " does not exist!\n";
    throw std::invalid_argument(ss.str());
  }

  size_t numParseThreads = salmonOpts.numParseThreads;
  std::cerr << "parseThreads = " << numParseThreads << "\n";
  bq = std::unique_ptr<BAMQueue<FragT>>(new BAMQueue<FragT>(
      alnFiles, state_.library, numParseThreads,
      salmonOpts.mappingCacheMemoryLimit));

  std::cerr << "Checking that provided alignment files have consistent "
               "headers . . . ";
  if (!salmon::utils::headersAreConsistent(bq->headers())) {
    std::stringstream ss;
    ss << "\nThe multiple alignment files provided had inconsistent "
          "headers.\n";
    ss << "Currently, we require that if multiple SAM/BAM files are "
          "provided,\n";
    ss << "they must have identical @SQ records.\n";
    throw std::invalid_argument(ss.str());
  }
  std::cerr << "done\n";

  AlignmentHeader* header = bq->header();
  aligner_ = salmon::bam_utils::inferAlignerFromHeader(header);

  phmap::flat_hash_set<std::string> decoys;
  if (aligner_ == salmon::bam_utils::AlignerDetails::PUFFERFISH) {
    for (decltype(header->nref) i = 0; i < header->nref; ++i) {
      AlignmentHeaderTag* tag;
      for (tag = header->ref[i].tag; tag; tag = tag->next) {
        if ((tag->len == 4) and (std::strncmp(tag->str, "DS:D", 4) == 0)) {
          decoys.insert(header->ref[i].name);
          break;
        }
      }
    }
  }

  if (!decoys.empty()) {
    bq->forceEndParsing();
    bq.reset();
    salmonOpts.jointLog->error(
        "Salmon is being run in alignment-mode with a SAM/BAM file that contains decoy\n"
        "sequences (marked as such during salmon indexing). This SAM/BAM file had {}\n"
        "such sequences tagged in the header. Since alignments to decoys are not\n"
        "intended for decoy-level quantification, this functionality is not currently\n"
        "supported.  If you wish to run salmon with this SAM/BAM file, please \n"
        "filter out / remove decoy transcripts (those tagged with `DS:D`) from the \n"
        "header, and all SAM/BAM records that represent alignments to decoys \n"
        "(those tagged with `XT:A:D`). If you believe you are receiving this message\n"
        "in error, please report this issue on GitHub.",
        decoys.size());
    salmonOpts.jointLog->flush();
    std::stringstream ss;
    ss << "\nCannot quantify from SAM/BAM file containing decoy transcripts or alignment records!\n";
    throw std::runtime_error(ss.str());
  }

  double alpha = 0.005;
  transcripts_.reserve(header->nref);
  for (decltype(header->nref) i = 0; i < header->nref; ++i) {
    transcripts_.emplace_back(i, header->ref[i].name, header->ref[i].len,
                              alpha);
  }

  FASTAParser fp(transcriptFile.string());
  fmt::print(stderr, "Populating targets from aln = {}, fasta = {} . . .",
             alnFiles.front().string(), transcriptFile_.string());
  fp.populateTargets(transcripts_, salmonOpts);

  std::vector<uint32_t> lengths;
  lengths.reserve(transcripts_.size());
  for (auto& txp : transcripts_) {
    lengths.push_back(txp.RefLength);
  }
  setTranscriptLengthClasses_(lengths, state_.posBiasFW.size());

  fmt::print(stderr, "done\n");
  clusters_.reset(new ClusterForest(transcripts_.size(), transcripts_));

  size_t maxFragLen = salmonOpts.fragLenDistMax;
  double meanFragLen = salmonOpts.fragLenDistPriorMean;
  double fragLenStd = salmonOpts.fragLenDistPriorSD;
  size_t fragLenKernelN = 4;
  double fragLenKernelP = 0.5;
  flDist_.reset(new FragmentLengthDistribution(1.0, maxFragLen, meanFragLen,
                                               fragLenStd, fragLenKernelN,
                                               fragLenKernelP, 1));

  alnMod_.reset(new AlignModelT(1.0, salmonOpts.numErrorBins));
  alnMod_->setLogger(salmonOpts.jointLog);

  if (libFmt.type == ReadType::SINGLE_END) {
    std::vector<double> logPMF;
    size_t minVal;
    size_t maxVal;
    flDist_->dumpPMF(logPMF, minVal, maxVal);
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

  salmon::utils::markAuxiliaryTargets(salmonOpts, transcripts_);
  NullFragmentFilter<FragT>* nff = nullptr;
  bool onlyProcessAmbiguousAlignments = false;
  bq->start(nff, onlyProcessAmbiguousAlignments);
}

template <typename FragT, typename EQBuilderT, typename AlignModelT>
AlignmentLibrary<FragT, EQBuilderT, AlignModelT>::AlignmentLibrary(
    std::vector<boost::filesystem::path>& alnFiles, LibraryFormat libFmt,
    SalmonOpts& salmonOpts, bool /*eqClassMode_*/,
    std::vector<std::string>& tnames, std::vector<double>& tefflens)
    : alignmentFiles_(alnFiles),
      state_(libFmt, salmonOpts.numConditionalGCBins, salmonOpts.numFragGCBins),
      transcripts_(std::vector<Transcript>()),
      eqBuilder_(salmonOpts.jointLog, salmonOpts.maxHashResizeThreads),
      quantificationPasses_(0) {
  namespace bfs = boost::filesystem;

  for (auto& alignmentFile : alignmentFiles_) {
    if (!bfs::exists(alignmentFile)) {
      std::stringstream ss;
      ss << "The provided eqClass file: " << alignmentFile
         << " does not exist!\n";
      throw std::invalid_argument(ss.str());
    }
  }

  double alpha = 0.005;
  for (size_t i = 0; i < tnames.size(); ++i) {
    transcripts_.emplace_back(i, tnames[i].c_str(), tefflens[i], true, alpha);
  }

  size_t maxFragLen = salmonOpts.fragLenDistMax;
  double meanFragLen = salmonOpts.fragLenDistPriorMean;
  double fragLenStd = salmonOpts.fragLenDistPriorSD;
  size_t fragLenKernelN = 4;
  double fragLenKernelP = 0.5;
  flDist_.reset(new FragmentLengthDistribution(1.0, maxFragLen, meanFragLen,
                                               fragLenStd, fragLenKernelN,
                                               fragLenKernelP, 1));

  alnMod_.reset(new AlignModelT(1.0, salmonOpts.numErrorBins));
  alnMod_->setLogger(salmonOpts.jointLog);
  salmon::utils::markAuxiliaryTargets(salmonOpts, transcripts_);
}

template <typename FragT, typename EQBuilderT, typename AlignModelT>
void AlignmentLibrary<FragT, EQBuilderT, AlignModelT>::updateTranscriptLengthsAtomic(
    std::atomic<bool>& done) {
  if (sl_.try_lock()) {
    if (!done) {
      auto fld = flDist_.get();
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

template <typename FragT, typename EQBuilderT, typename AlignModelT>
void AlignmentLibrary<FragT, EQBuilderT, AlignModelT>::setTranscriptLengthClasses_(
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
