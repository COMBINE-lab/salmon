#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

// Our includes
#include "ClusterForest.hpp"
#include "DistributionUtils.hpp"
#include "EquivalenceClassBuilder.hpp"
#include "FragmentLengthDistribution.hpp"
#include "FragmentStartPositionDistribution.hpp"
#include "GCFragModel.hpp"
#include "ReadKmerDist.hpp"
#include "ReadLibrary.hpp"
#include "SBModel.hpp"
#include "SalmonIndex.hpp"
#include "SalmonOpts.hpp"
#include "SalmonUtils.hpp"
// #include "SequenceBiasModel.hpp"
#include "SimplePosBias.hpp"
#include "SpinLock.hpp" // RapMap's with try_lock
#include "Transcript.hpp"
#include "UtilityFunctions.hpp"

// Logger includes
#include "spdlog/spdlog.h"

// Boost includes
#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

// Cereal includes
#include "cereal/archives/json.hpp"

// Standard includes
#include <fstream>
#include <memory>
#include <vector>

/**
 *  This class represents a library of alignments used to quantify
 *  a set of target transcripts.  The AlignmentLibrary contains info
 *  about both the alignment file and the target sequence (transcripts).
 *  It is used to group them together and track information about them
 *  during the quantification procedure.
 */
template <typename EQBuilderT>
class ReadExperiment {

public:
  ReadExperiment(std::vector<ReadLibrary>& readLibraries,
                 // const boost::filesystem::path& transcriptFile,
                 const boost::filesystem::path& indexDirectory,
                 SalmonOpts& sopt)
      : readLibraries_(readLibraries),
        // transcriptFile_(transcriptFile),
        transcripts_(std::vector<Transcript>()), totalAssignedFragments_(0),
        fragStartDists_(5), posBiasFW_(5), posBiasRC_(5), posBiasExpectFW_(5),
        posBiasExpectRC_(5), /*seqBiasModel_(1.0),*/ eqBuilder_(sopt.jointLog, sopt.maxHashResizeThreads),
        expectedBias_(constExprPow(4, readBias_[0].getK()), 1.0),
        expectedGC_(sopt.numConditionalGCBins, sopt.numFragGCBins,
                    distribution_utils::DistributionSpace::LOG),
        observedGC_(sopt.numConditionalGCBins, sopt.numFragGCBins,
                    distribution_utils::DistributionSpace::LOG) {
    namespace bfs = boost::filesystem;

    // Make sure the read libraries are valid.
    for (auto& rl : readLibraries_) {
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

    if (readLibraries_.front().getFormat().type == ReadType::SINGLE_END) {
      // Convert the PMF to non-log scale
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

      // Create the non-logged distribution.
      // Here, we multiply by 100 to discourage small
      // numbers in the correctionFactorsfromCounts call
      // below.
      std::vector<double> pmf(maxVal + 1, 0.0);
      for (size_t i = minVal; i < maxVal; ++i) {
        pmf[i] = 100.0 * std::exp(logPMF[i - minVal]);
      }

      using distribution_utils::DistributionSpace;
      // We compute the factors in linear space (since we've de-logged the pmf)
      conditionalMeans_ = distribution_utils::correctionFactorsFromMass(
          pmf, DistributionSpace::LINEAR);
    }

    // Make sure the transcript file exists.
    /*
    if (!bfs::exists(transcriptFile_)) {
        std::stringstream ss;
        ss << "The provided transcript file: " << transcriptFile_ <<
            " does not exist!\n";
        throw std::invalid_argument(ss.str());
    }
    */

    // ==== Figure out the index type
    boost::filesystem::path versionPath = indexDirectory / "versionInfo.json";
    SalmonIndexVersionInfo versionInfo;
    versionInfo.load(versionPath);
    if (versionInfo.indexVersion() == 0) {
      fmt::MemoryWriter infostr;
      infostr << "Error: The index version file " << versionPath.string()
              << " doesn't seem to exist.  Please try re-building the salmon "
                 "index.";
      throw std::invalid_argument(infostr.str());
    }
    // Check index version compatibility here
    auto indexType = versionInfo.indexType();
    // ==== Figure out the index type

    salmonIndex_.reset(new SalmonIndex(sopt.jointLog, indexType));
    salmonIndex_->load(indexDirectory);

    // Now we'll have either an FMD-based index or a QUASI index
    // dispatch on the correct type.
    fmt::MemoryWriter infostr;

    switch (salmonIndex_->indexType()) {
    case SalmonIndexType::QUASI:
      infostr << "Error: This version of salmon does not support a RapMap-based index.";
      throw std::invalid_argument(infostr.str());
      break;
    case SalmonIndexType::FMD:
      infostr << "Error: This version of salmon does not support the FMD index mode.";
      throw std::invalid_argument(infostr.str());
      break;
    case SalmonIndexType::PUFF:
      if (salmonIndex_->isSparse()){
        loadTranscriptsFromPuff(salmonIndex_->puffSparseIndex(), sopt);
      } else {
        loadTranscriptsFromPuff(salmonIndex_->puffIndex(), sopt);
      }
    }

    salmon::utils::markAuxiliaryTargets(sopt, transcripts_);

    // Create the cluster forest for this set of transcripts
    clusters_.reset(new ClusterForest(transcripts_.size(), transcripts_));
  }

  EQBuilderT& equivalenceClassBuilder() { return eqBuilder_; }

  std::string getIndexSeqHash256() const { return salmonIndex_->seqHash256(); }
  std::string getIndexNameHash256() const { return salmonIndex_->nameHash256(); }
  std::string getIndexSeqHash512() const { return salmonIndex_->seqHash512(); }
  std::string getIndexNameHash512() const { return salmonIndex_->nameHash512(); }
  std::string getIndexDecoySeqHash256() const { return salmonIndex_->decoySeqHash256(); }
  std::string getIndexDecoyNameHash256() const { return salmonIndex_->decoyNameHash256(); }

  std::vector<Transcript>& transcripts() { return transcripts_; }
  const std::vector<Transcript>& transcripts() const { return transcripts_; }

  const std::vector<double>& condMeans() const { return conditionalMeans_; }

  void updateTranscriptLengthsAtomic(std::atomic<bool>& done) {
    if (sl_.try_lock()) {
      if (!done) {
        auto fld = fragLengthDist_.get();
        // Convert the PMF to non-log scale
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

        // Create the non-logged distribution.
        // Here, we multiply by 100 to discourage small
        // numbers in the correctionFactorsfromCounts call
        // below.
        std::vector<double> pmf(maxVal + 1, 0.0);
        for (size_t i = minVal; i < maxVal; ++i) {
          pmf[i] = 100.0 * std::exp(logPMF[i - minVal]);
        }

        using distribution_utils::DistributionSpace;
        // We compute the factors in linear space (since we've de-logged the
        // pmf)
        auto correctionFactors = distribution_utils::correctionFactorsFromMass(
            pmf, DistributionSpace::LINEAR);
        // Since we'll continue treating effective lengths in log space,
        // populate them as such
        distribution_utils::computeSmoothedEffectiveLengths(
            pmf.size(), transcripts_, correctionFactors,
            DistributionSpace::LOG);

        /*
        // Update the effective length of *every* transcript
        for( auto& t : transcripts_ ) {
            t.updateEffectiveLength(logPMF, logFLDMean, minVal, maxVal);
        }
        */
        // then declare that we are done
        done = true;
      }
      sl_.unlock();
    }
  }

  uint64_t numAssignedFragments() { return numAssignedFragments_; }
  uint64_t numMappedFragments() const { return numAssignedFragments_; }

  uint64_t upperBoundHits() { return upperBoundHits_; }
  void setUpperBoundHits(uint64_t ubh) { upperBoundHits_ = ubh; }

  std::atomic<uint64_t>& numAssignedFragmentsAtomic() {
    return numAssignedFragments_;
  }

  void setNumObservedFragments(uint64_t numObserved) {
    numObservedFragments_ = numObserved;
  }

  void updateShortFrags(salmon::utils::ShortFragStats& fs) {
    sl_.lock();
    shortFragStats_.numTooShort += fs.numTooShort;
    shortFragStats_.shortest = (fs.shortest < shortFragStats_.shortest)
                                   ? fs.shortest
                                   : shortFragStats_.shortest;
    sl_.unlock();
  }

  salmon::utils::ShortFragStats getShortFragStats() const {
    return shortFragStats_;
  }

  uint64_t numObservedFragments() const { return numObservedFragments_; }

  double mappingRate() {
    if (quantificationPasses_ > 0) {
      return static_cast<double>(numAssignedFragsInFirstPass_) /
             numObservedFragsInFirstPass_;
    } else {
      return static_cast<double>(numAssignedFragments_) / numObservedFragments_;
    }
  }

  SalmonIndex* getIndex() { return salmonIndex_.get(); }

  template <typename PuffIndexT>
  void loadTranscriptsFromPuff(PuffIndexT* idx_, const SalmonOpts& sopt) {
    // Get the list of reference names
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
    size_t numIndexedRefs = idx_->getIndexedRefCount();
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
      // copy over the length, then we're done.
      transcripts_.emplace_back(id, name, len, alpha);
      auto& txp = transcripts_.back();
      txp.setCompleteLength(completeRefLengths[i]);

      // TODO : PF_INTEGRATION
      // We won't every have the sequence for a decoy
      // NOTE : The below is a huge hack right now.  Since we have
      // the whole reference sequence in 2-bit in the index, there is
      // really no reason to convert to ASCII and keep a duplicate
      // copy around.
      if (!isDecoy and !isShort and (sopt.biasCorrect or sopt.gcBiasCorrect)) {
        auto tid = i - numShort;
        int64_t refAccPos = tid > 0 ? refAccumLengths[tid - 1] : 0;
        int64_t refTotalLength = refAccumLengths[tid] - refAccPos;
        if (len != refTotalLength) {
          log->warn("len : {:n}, but txp.RefLength : {:n} :: refTotalLength : {:n}", len, txp.RefLength, refTotalLength);
        }
        char* tseq = pufferfish::util::getRefSeqOwned(allRefSeq, refAccPos, refTotalLength);
        txp.setSequenceOwned(tseq, sopt.gcBiasCorrect, sopt.reduceGCMemory);
      }
      txp.setDecoy(isDecoy);
      numShort += isShort ? 1 : 0;
      if (isDecoy) { 
        ++numDecoys_; 
      } else { // only use this reference for length class computation if not a decoy
        lengths.push_back(txp.RefLength);
      }
    }
    sopt.jointLog->info("Number of decoys : {:n}", numDecoys_);
    auto firstDecoyIndex = idx_->firstDecoyIndex();
    if (firstDecoyIndex < numRecords) {
      sopt.jointLog->info("First decoy index : {:n} ", firstDecoyIndex);
    }
    //std::exit(1);
    // ====== Done loading the transcripts from file
    setTranscriptLengthClasses_(lengths, posBiasFW_.size());
  }

  template <typename QuasiIndexT>
  void loadTranscriptsFromQuasi(QuasiIndexT* idx_, const SalmonOpts& sopt) {
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
      // copy over the length, then we're done.
      transcripts_.emplace_back(id, name, len, alpha);
      auto& txp = transcripts_.back();
      txp.setCompleteLength(idx_->txpCompleteLens[i]);
      // The transcript sequence
      // auto txpSeq = idx_->seq.substr(idx_->txpOffsets[i], len);
      // Set the transcript sequence
      txp.setSequenceBorrowed(idx_->seq.c_str() + idx_->txpOffsets[i],
                              sopt.gcBiasCorrect, sopt.reduceGCMemory);
      txp.setDecoy(isDecoy);
      if (isDecoy) { ++numDecoys_; }

      lengths.push_back(txp.RefLength);
      /*
      // Length classes taken from
      //
      https://github.com/cole-trapnell-lab/cufflinks/blob/master/src/biascorrection.cpp
      // ======
      // Roberts, Adam, et al.
      // "Improving RNA-Seq expression estimates by correcting for fragment
      bias."
      // Genome Biol 12.3 (2011): R22.
      // ======
      // perhaps, define these in a more data-driven way
      if (txp.RefLength <= 791) {
      txp.lengthClassIndex(0);
      } else if (txp.RefLength <= 1265) {
      txp.lengthClassIndex(1);
      } else if (txp.RefLength <= 1707) {
      txp.lengthClassIndex(2);
      } else if (txp.RefLength <= 2433) {
      txp.lengthClassIndex(3);
      } else {
      txp.lengthClassIndex(4);
      }
      */
    }
    // ====== Done loading the transcripts from file
    setTranscriptLengthClasses_(lengths, posBiasFW_.size());
  }

  void dropDecoyTranscripts() {
    if (numDecoys_ > 0) {
      size_t numValidTargets = transcripts_.size() - numDecoys_;
      transcripts_.resize(numValidTargets);
    }
  }

  template <typename CallbackT>
  bool processReads(const uint32_t& numThreads, const SalmonOpts& sopt,
                    CallbackT& processReadLibrary) {
    std::atomic<bool> burnedIn{
        totalAssignedFragments_ + numAssignedFragments_ >= sopt.numBurninFrags};
    for (auto& rl : readLibraries_) {
      processReadLibrary(rl, salmonIndex_.get(), transcripts_, clusterForest(),
                         *(fragLengthDist_.get()), numAssignedFragments_,
                         numThreads, burnedIn);
    }
    return true;
  }

  ~ReadExperiment() {
    // ---- Get rid of things we no longer need --------
  }

  ClusterForest& clusterForest() { return *clusters_.get(); }

  std::string readFilesAsString() {
    std::stringstream sstr;
    size_t ln{0};
    size_t numReadLibraries{readLibraries_.size()};

    for (auto& rl : readLibraries_) {
      sstr << rl.readFilesAsString();
      if (ln++ < numReadLibraries) {
        sstr << "; ";
      }
    }
    return sstr.str();
  }

  uint64_t numAssignedFragsInFirstPass() {
    return numAssignedFragsInFirstPass_;
  }

  uint64_t numObservedFragsInFirstPass() {
    return numObservedFragsInFirstPass_;
  }

  double effectiveMappingRate() const { return effectiveMappingRate_; }

  void setEffectiveMappingRate(double emr) { effectiveMappingRate_ = emr; }

  std::vector<FragmentStartPositionDistribution>&
  fragmentStartPositionDistributions() {
    return fragStartDists_;
  }

  void computePolyAPositions() {
    for (auto& t : transcripts_) {
      t.computePolyAPositions();
    }
  }

  // SequenceBiasModel& sequenceBiasModel() { return seqBiasModel_; }

  bool softReset() {
    if (quantificationPasses_ == 0) {
      numAssignedFragsInFirstPass_ = numAssignedFragments_;
      numObservedFragsInFirstPass_ = numObservedFragments_;
    }
    numObservedFragments_ = 0;
    totalAssignedFragments_ += numAssignedFragments_;
    numAssignedFragments_ = 0;
    quantificationPasses_++;
    return true;
  }

  bool reset() {
    namespace bfs = boost::filesystem;
    for (auto& rl : readLibraries_) {
      if (!rl.isRegularFile()) {
        return false;
      }
    }

    if (quantificationPasses_ == 0) {
      numAssignedFragsInFirstPass_ = numAssignedFragments_;
      numObservedFragsInFirstPass_ = numObservedFragments_;
    }

    numObservedFragments_ = 0;
    totalAssignedFragments_ += numAssignedFragments_;
    numAssignedFragments_ = 0;
    quantificationPasses_++;
    return true;
  }

  void summarizeLibraryTypeCounts(boost::filesystem::path& opath) {
    LibraryFormat fmt1(ReadType::SINGLE_END, ReadOrientation::NONE,
                       ReadStrandedness::U);
    LibraryFormat fmt2(ReadType::SINGLE_END, ReadOrientation::NONE,
                       ReadStrandedness::U);

    std::ofstream os(opath.string());
    cereal::JSONOutputArchive oa(os);

    // std::ofstream ofile(opath.string());

    fmt::MemoryWriter errstr;

    auto log = spdlog::get("jointLog");

    uint64_t numFmt1{0};
    uint64_t numFmt2{0};
    uint64_t numAgree{0};
    uint64_t numDisagree{0};
    uint64_t numDisagreeUnstranded{0};
    uint64_t numDisagreeStranded{0};

    for (auto& rl : readLibraries_) {
      auto fmt = rl.format();
      auto& counts = rl.libTypeCounts();

      oa(cereal::make_nvp("read_files", rl.readFilesAsString()));
      std::string expectedFormat = rl.format().toString();
      oa(cereal::make_nvp("expected_format", expectedFormat));

      double compatFragmentRatio =
          rl.numCompat() / static_cast<double>(numAssignedFragments_);
      oa(cereal::make_nvp("compatible_fragment_ratio", compatFragmentRatio));
      oa(cereal::make_nvp("num_compatible_fragments", rl.numCompat()));
      oa(cereal::make_nvp("num_assigned_fragments",
                          numAssignedFragments_.load()));

      std::vector<ReadStrandedness> strands;
      // How we should define sense and antisense 
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

      // fmt1 is :
      // reads that arise from sense, or 
      // fragments whose first read arises from sense
      fmt1.type = fmt.type;
      fmt1.orientation = fmt.orientation;
      fmt1.strandedness = strands[0];
      // fmt 2 is :
      // reads that arise from antisense, or
      // fragments whose first read arises from antisense
      fmt2.type = fmt.type;
      fmt2.orientation = fmt.orientation;
      fmt2.strandedness = strands[1];

      numFmt1 = 0;
      numFmt2 = 0;

      for (size_t i = 0; i < counts.size(); ++i) {
        // This counts agree or disagree when the 
        // provided or detected library type is *not*
        // unstranded
        if (i == fmt.formatID()) {
          numAgree = counts[i];
        } else {
          numDisagreeStranded += counts[i];
        }

        // Regardless of if the provided or detected
        // library type is unstranded, count the 
        // sense and antisense compatible reads
        if (i == fmt1.formatID()) {
          numFmt1 = counts[i];
        } else if (i == fmt2.formatID()) {
          numFmt2 = counts[i];
        } else {
          // collect this number for if
          // we have an unstranded library type
          numDisagreeUnstranded += counts[i];
        }
      }

      double ratio = 0.0;

      // If the provided or detected library type is *unstranded*
      if (fmt.strandedness == ReadStrandedness::U) {
        // overwrite numAgree, since that was computed for 
        // a stranded library 
        numAgree = numFmt1 + numFmt2;
        numDisagree = numDisagreeUnstranded;
        ratio = numAgree > 0 ? (static_cast<double>(numFmt1) / (numFmt1 + numFmt2)) : 0.0;

        if (numAgree == 0) {
          errstr << "NOTE: Read Lib [" << rl.readFilesAsString() << "] :\n";
          errstr << "\nFound no concordant and consistent mappings. "
                    "If this is a paired-end library, are you sure the reads are properly paired? "
                    "check the file: " << opath.string() << " for details\n";
          log->warn(errstr.str());
          errstr.clear();
        } else if (std::abs(ratio - 0.5) > 0.01) {
          // check that we have a similar number of mappings in both
          // directions and then aggregate the forward and
          // reverse counts.
          errstr << "NOTE: Read Lib [" << rl.readFilesAsString() << "] :\n";
          errstr << "\nDetected a *potential* strand bias > 1\% in an "
                    "unstranded protocol "
                 << "check the file: " << opath.string() << " for details\n";

          log->warn(errstr.str());
          errstr.clear();
        }
      } else {
        // If the provided or detected library type is *stranded*
        numDisagree = numDisagreeStranded;
        ratio = numAgree > 0 ? (static_cast<double>(numFmt1) / (numFmt1 + numFmt2)) : 0.0;
      } // end else

      oa(cereal::make_nvp("num_frags_with_concordant_consistent_mappings", numAgree));
      oa(cereal::make_nvp("num_frags_with_inconsistent_or_orphan_mappings", numDisagree));
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

  std::vector<ReadLibrary>& readLibraries() { return readLibraries_; }
  const std::vector<ReadLibrary>& readLibraries() const {
    return readLibraries_;
  }
  FragmentLengthDistribution* fragmentLengthDistribution() const {
    return fragLengthDist_.get();
  }

  void setGCFracForward(double fracForward) { gcFracFwd_ = fracForward; }

  double gcFracFwd() const { return gcFracFwd_; }
  double gcFracRC() const { return 1.0 - gcFracFwd_; }

  std::vector<double>& expectedSeqBias() { return expectedBias_; }

  const std::vector<double>& expectedSeqBias() const { return expectedBias_; }

  void setExpectedGCBias(const GCFragModel& expectedBiasIn) {
    expectedGC_ = expectedBiasIn;
  }

  GCFragModel& expectedGCBias() { return expectedGC_; }

  const GCFragModel& expectedGCBias() const { return expectedGC_; }

  const GCFragModel& observedGC() const { return observedGC_; }

  GCFragModel& observedGC() { return observedGC_; }

  std::vector<SimplePosBias>& posBias(salmon::utils::Direction dir) {
    return (dir == salmon::utils::Direction::FORWARD) ? posBiasFW_ : posBiasRC_;
  }
  const std::vector<SimplePosBias>&
  posBias(salmon::utils::Direction dir) const {
    return (dir == salmon::utils::Direction::FORWARD) ? posBiasFW_ : posBiasRC_;
  }

  std::vector<SimplePosBias>& posBiasExpected(salmon::utils::Direction dir) {
    return (dir == salmon::utils::Direction::FORWARD) ? posBiasExpectFW_
                                                      : posBiasExpectRC_;
  }
  const std::vector<SimplePosBias>&
  posBiasExpected(salmon::utils::Direction dir) const {
    return (dir == salmon::utils::Direction::FORWARD) ? posBiasExpectFW_
                                                      : posBiasExpectRC_;
  }

  ReadKmerDist<6, std::atomic<uint32_t>>&
  readBias(salmon::utils::Direction dir) {
    return (dir == salmon::utils::Direction::FORWARD) ? readBias_[0]
                                                      : readBias_[1];
  }
  const ReadKmerDist<6, std::atomic<uint32_t>>&
  readBias(salmon::utils::Direction dir) const {
    return (dir == salmon::utils::Direction::FORWARD) ? readBias_[0]
                                                      : readBias_[1];
  }

  SBModel& readBiasModelObserved(salmon::utils::Direction dir) {
    return (dir == salmon::utils::Direction::FORWARD)
               ? readBiasModelObserved_[0]
               : readBiasModelObserved_[1];
  }
  const SBModel& readBiasModelObserved(salmon::utils::Direction dir) const {
    return (dir == salmon::utils::Direction::FORWARD)
               ? readBiasModelObserved_[0]
               : readBiasModelObserved_[1];
  }

  SBModel& readBiasModelExpected(salmon::utils::Direction dir) {
    return (dir == salmon::utils::Direction::FORWARD)
               ? readBiasModelExpected_[0]
               : readBiasModelExpected_[1];
  }
  const SBModel& readBiasModelExpected(salmon::utils::Direction dir) const {
    return (dir == salmon::utils::Direction::FORWARD)
               ? readBiasModelExpected_[0]
               : readBiasModelExpected_[1];
  }
  void setReadBiasModelExpected(SBModel&& model, salmon::utils::Direction dir) {
    size_t idx = (dir == salmon::utils::Direction::FORWARD) ? 0 : 1;
    readBiasModelExpected_[idx] = std::move(model);
  }

  const std::vector<uint32_t>& getLengthQuantiles() const {
    return lengthQuantiles_;
  }

  uint64_t getNumDecoys() const {
    return numDecoys_;
  }

  salmon::utils::DuplicateTargetStatus index_retains_duplicates() const { 
    return salmonIndex_->index_retains_duplicates(); 
  }

private:
  void setTranscriptLengthClasses_(std::vector<uint32_t>& lengths,
                                   size_t nbins) {
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
        // Find the proper quantile
        lengthQuantiles_.push_back(*(lengths.begin() + ind));
      }
    } else {
      lengthQuantiles_.clear();
      lengthQuantiles_.reserve(n);
      std::sort(lengths.begin(), lengths.end());
      for (auto l : lengths) {
        lengthQuantiles_.push_back(l);
      }
      posBiasFW_.resize(n);
      posBiasRC_.resize(n);
      posBiasExpectFW_.resize(n);
      posBiasExpectRC_.resize(n);
    }

    auto qb = lengthQuantiles_.begin();
    auto qe = lengthQuantiles_.end();
    auto maxQuant = std::distance(qb, qe) - 1;
    for (auto& t : transcripts_) {
      auto ind = std::min(
          maxQuant, std::distance(qb, std::upper_bound(qb, qe, t.RefLength)));
      // the index is the smallest quantile longer than this length
      t.lengthClassIndex(ind);
    }
  }

  /**
   * The file from which the alignments will be read.
   * This can be a SAM or BAM file, and can be a regular
   * file or a fifo.
   */
  std::vector<ReadLibrary> readLibraries_;
  /**
   * The file from which the transcripts are read.
   * This is expected to be a FASTA format file.
   */
  // boost::filesystem::path transcriptFile_;
  /**
   * The targets (transcripts) to be quantified.
   */
  std::vector<Transcript> transcripts_;
  /**
   * The index we've built on the set of transcripts.
   */
  std::unique_ptr<SalmonIndex> salmonIndex_{nullptr};
  /**
   * The cluster forest maintains the dynamic relationship
   * defined by transcripts and reads --- if two transcripts
   * share an ambiguously mapped read, then they are placed
   * in the same cluster.
   */
  std::unique_ptr<ClusterForest> clusters_;
  /**
   *
   *
   */
  std::vector<FragmentStartPositionDistribution> fragStartDists_;

  // SequenceBiasModel seqBiasModel_;

  /** Keeps track of the number of passes that have been
   *  made through the alignment file.
   */
  salmon::utils::ShortFragStats shortFragStats_;
  std::atomic<uint64_t> numObservedFragments_{0};
  std::atomic<uint64_t> numAssignedFragments_{0};
  uint64_t totalAssignedFragments_{0};
  size_t quantificationPasses_{0};
  uint64_t numAssignedFragsInFirstPass_{0};
  uint64_t numObservedFragsInFirstPass_{0};
  uint64_t upperBoundHits_{0};
  double effectiveMappingRate_{0.0};
  SpinLock sl_;
  std::unique_ptr<FragmentLengthDistribution> fragLengthDist_;
  EQBuilderT eqBuilder_;

  /** Positional bias things**/
  std::vector<uint32_t> lengthQuantiles_;
  std::vector<SimplePosBias> posBiasFW_;
  std::vector<SimplePosBias> posBiasRC_;
  std::vector<SimplePosBias> posBiasExpectFW_;
  std::vector<SimplePosBias> posBiasExpectRC_;

  /** GC-fragment bias things **/
  // One bin for each percentage GC content
  double gcFracFwd_{-1.0};
  GCFragModel observedGC_;
  GCFragModel expectedGC_;

  /** Sequence specific bias things **/
  // Since multiple threads can touch this dist, we
  // need atomic counters.
  std::array<ReadKmerDist<6, std::atomic<uint32_t>>, 2> readBias_;
  std::array<SBModel, 2> readBiasModelObserved_;
  std::array<SBModel, 2> readBiasModelExpected_;
  // std::array<std::vector<double>, 2> expectedBias_;
  std::vector<double> expectedBias_;
  std::vector<double> conditionalMeans_;

  uint64_t numDecoys_{0};
};

#endif // EXPERIMENT_HPP
