#ifndef ALIGNMENT_LIBRARY_HPP
#define ALIGNMENT_LIBRARY_HPP

extern "C" {
#include "io_lib/os.h"
#include "io_lib/scram.h"
#undef max
#undef min
}

// Our includes
#include "AlignmentGroup.hpp"
#include "AlignmentModel.hpp"
#include "BAMQueue.hpp"
#include "BAMUtils.hpp"
#include "ClusterForest.hpp"
#include "DistributionUtils.hpp"
#include "EquivalenceClassBuilder.hpp"
#include "FASTAParser.hpp"
#include "FragmentLengthDistribution.hpp"
#include "FragmentStartPositionDistribution.hpp"
#include "GCFragModel.hpp"
#include "LibraryFormat.hpp"
#include "LibraryTypeDetector.hpp"
#include "ReadKmerDist.hpp"
#include "SBModel.hpp"
#include "SalmonOpts.hpp"
#include "SalmonUtils.hpp"
#include "SimplePosBias.hpp"
#include "SpinLock.hpp" // From pufferfish, with try_lock
#include "Transcript.hpp"
#include "concurrentqueue.h"
#include "parallel_hashmap/phmap.h"

// Boost includes
#include <boost/filesystem.hpp>

// Standard includes
#include <functional>
#include <memory>
#include <vector>
#include <stdexcept>

template <typename T> class NullFragmentFilter;

/**
 *  This class represents a library of alignments used to quantify
 *  a set of target transcripts.  The AlignmentLibrary contains info
 *  about both the alignment file and the target sequence (transcripts).
 *  It is used to group them together and track information about them
 *  during the quantification procedure.
 */
template <typename FragT, typename EQBuilderT> class AlignmentLibrary {

public:
  AlignmentLibrary(std::vector<boost::filesystem::path>& alnFiles,
                   const boost::filesystem::path& transcriptFile,
                   LibraryFormat libFmt, SalmonOpts& salmonOpts)
      : alignmentFiles_(alnFiles), transcriptFile_(transcriptFile),
        libFmt_(libFmt), transcripts_(std::vector<Transcript>()),
        fragStartDists_(5), posBiasFW_(5), posBiasRC_(5), posBiasExpectFW_(5),
        posBiasExpectRC_(5), /*seqBiasModel_(1.0),*/
        eqBuilder_(salmonOpts.jointLog, salmonOpts.maxHashResizeThreads), quantificationPasses_(0),
        expectedBias_(constExprPow(4, readBias_[0].getK()), 1.0),
        expectedGC_(salmonOpts.numConditionalGCBins, salmonOpts.numFragGCBins,
                    distribution_utils::DistributionSpace::LOG),
        observedGC_(salmonOpts.numConditionalGCBins, salmonOpts.numFragGCBins,
                    distribution_utils::DistributionSpace::LOG) {
    namespace bfs = boost::filesystem;

    // Make sure the alignment file exists.
    for (auto& alignmentFile : alignmentFiles_) {
      if (!bfs::exists(alignmentFile)) {
        std::stringstream ss;
        ss << "The provided alignment file: " << alignmentFile
           << " does not exist!\n";
        throw std::invalid_argument(ss.str());
      }
    }

    // Make sure the transcript file exists.
    if (!bfs::exists(transcriptFile_)) {
      std::stringstream ss;
      ss << "The provided transcript file: " << transcriptFile_
         << " does not exist!\n";
      throw std::invalid_argument(ss.str());
    }

    // The alignment file existed, so create the alignment queue
    size_t numParseThreads = salmonOpts.numParseThreads;
    std::cerr << "parseThreads = " << numParseThreads << "\n";
    bq = std::unique_ptr<BAMQueue<FragT>>(
        new BAMQueue<FragT>(alnFiles, libFmt_, numParseThreads,
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

    SAM_hdr* header = bq->header();

    // Figure out aligner information from the header if we can
    aligner_ = salmon::bam_utils::inferAlignerFromHeader(header);

    // in this case check for decoys and make a list of their names
    phmap::flat_hash_set<std::string> decoys;
    if (aligner_ == salmon::bam_utils::AlignerDetails::PUFFERFISH) {
     // for each reference
     for (decltype(header->nref) i = 0; i < header->nref; ++i) {
       // for each tag 
       SAM_hdr_tag *tag;
	     for (tag = header->ref[i].tag; tag; tag = tag->next) {
         // if this tag marks it as a decoy
         if ((tag->len == 4) and (std::strncmp(tag->str, "DS:D", 4) == 0)) {
           decoys.insert(header->ref[i].name);
           break;
         } // end if decoy tag

       } // end for each tag
      } // end for each referecne
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
        "in error, please report this issue on GitHub.", decoys.size());
        salmonOpts.jointLog->flush();
        std::stringstream ss;
        ss << "\nCannot quantify from SAM/BAM file containing decoy transcripts or alignment records!\n";
        throw std::runtime_error(ss.str());
    }

    // The transcript file existed, so load up the transcripts
    double alpha = 0.005;
    // we know how many we will have, so reserve the space for 
    // them.
    transcripts_.reserve(header->nref);
    for (decltype(header->nref) i = 0; i < header->nref; ++i) {
      transcripts_.emplace_back(i, header->ref[i].name, header->ref[i].len,
                                alpha);
    }

    FASTAParser fp(transcriptFile.string());

    fmt::print(stderr, "Populating targets from aln = {}, fasta = {} . . .",
               alnFiles.front(), transcriptFile_);
    fp.populateTargets(transcripts_, salmonOpts);
    /*
for (auto& txp : transcripts_) {
    // Length classes taken from
    // ======
    // Roberts, Adam, et al.
    // "Improving RNA-Seq expression estimates by correcting for fragment bias."
    // Genome Biol 12.3 (2011): R22.
    // ======
    // perhaps, define these in a more data-driven way
    if (txp.RefLength <= 1334) {
        txp.lengthClassIndex(0);
    } else if (txp.RefLength <= 2104) {
        txp.lengthClassIndex(0);
    } else if (txp.RefLength <= 2988) {
        txp.lengthClassIndex(0);
    } else if (txp.RefLength <= 4389) {
        txp.lengthClassIndex(0);
    } else {
        txp.lengthClassIndex(0);
    }
}
*/

    std::vector<uint32_t> lengths;
    lengths.reserve(transcripts_.size());
    for (auto& txp : transcripts_) {
      lengths.push_back(txp.RefLength);
    }
    setTranscriptLengthClasses_(lengths, posBiasFW_.size());

    fmt::print(stderr, "done\n");

    // Create the cluster forest for this set of transcripts
    clusters_.reset(new ClusterForest(transcripts_.size(), transcripts_));

    // Initialize the fragment length distribution
    size_t maxFragLen = salmonOpts.fragLenDistMax;
    double meanFragLen = salmonOpts.fragLenDistPriorMean;
    double fragLenStd = salmonOpts.fragLenDistPriorSD;
    size_t fragLenKernelN = 4;
    double fragLenKernelP = 0.5;
    flDist_.reset(new FragmentLengthDistribution(1.0, maxFragLen, meanFragLen,
                                                 fragLenStd, fragLenKernelN,
                                                 fragLenKernelP, 1));

    alnMod_.reset(new AlignmentModel(1.0, salmonOpts.numErrorBins));
    alnMod_->setLogger(salmonOpts.jointLog);

    if (libFmt.type == ReadType::SINGLE_END) {
      // Convert the PMF to non-log scale
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

    salmon::utils::markAuxiliaryTargets(salmonOpts, transcripts_);

    // Start parsing the alignments
    NullFragmentFilter<FragT>* nff = nullptr;
    bool onlyProcessAmbiguousAlignments = false;
    bq->start(nff, onlyProcessAmbiguousAlignments);
  }

  AlignmentLibrary(std::vector<boost::filesystem::path>& alnFiles,
                   LibraryFormat libFmt, SalmonOpts& salmonOpts,
                   bool /*eqClassMode_*/, std::vector<std::string>& tnames,
                   std::vector<double>& tefflens)
      : alignmentFiles_(alnFiles),
        libFmt_(libFmt), transcripts_(std::vector<Transcript>()),
        fragStartDists_(5), posBiasFW_(5), posBiasRC_(5), posBiasExpectFW_(5),
        posBiasExpectRC_(5), /*seqBiasModel_(1.0),*/
        eqBuilder_(salmonOpts.jointLog, salmonOpts.maxHashResizeThreads), quantificationPasses_(0),
        expectedBias_(constExprPow(4, readBias_[0].getK()), 1.0),
        expectedGC_(salmonOpts.numConditionalGCBins, salmonOpts.numFragGCBins,
                    distribution_utils::DistributionSpace::LOG),
        observedGC_(salmonOpts.numConditionalGCBins, salmonOpts.numFragGCBins,
                    distribution_utils::DistributionSpace::LOG) {
    namespace bfs = boost::filesystem;

    // Make sure the alignment file exists.
    for (auto& alignmentFile : alignmentFiles_) {
      if (!bfs::exists(alignmentFile)) {
        std::stringstream ss;
        ss << "The provided eqClass file: " << alignmentFile
           << " does not exist!\n";
        throw std::invalid_argument(ss.str());
      }
    }

    // The transcript file existed, so load up the transcripts
    double alpha = 0.005;
    for (size_t i = 0; i < tnames.size(); ++i) {
      transcripts_.emplace_back(i, tnames[i].c_str(), tefflens[i], true, alpha);
    }

    // Initialize the fragment length distribution
    size_t maxFragLen = salmonOpts.fragLenDistMax;
    double meanFragLen = salmonOpts.fragLenDistPriorMean;
    double fragLenStd = salmonOpts.fragLenDistPriorSD;
    size_t fragLenKernelN = 4;
    double fragLenKernelP = 0.5;
    flDist_.reset(new FragmentLengthDistribution(1.0, maxFragLen, meanFragLen,
                                                 fragLenStd, fragLenKernelN,
                                                 fragLenKernelP, 1));

    alnMod_.reset(new AlignmentModel(1.0, salmonOpts.numErrorBins));
    alnMod_->setLogger(salmonOpts.jointLog);
    salmon::utils::markAuxiliaryTargets(salmonOpts, transcripts_);
  }

  EQBuilderT& equivalenceClassBuilder() { return eqBuilder_; }

  std::string getIndexSeqHash256() const { return ""; }
  std::string getIndexNameHash256() const { return ""; }
  std::string getIndexSeqHash512() const { return ""; }
  std::string getIndexNameHash512() const { return ""; }
  std::string getIndexDecoySeqHash256() const { return ""; }
  std::string getIndexDecoyNameHash256() const { return ""; }

  /**
   * Return true if this read library is for paired-end reads and false
   * otherwise.
   */
  bool isPairedEnd() { return (libFmt_.type == ReadType::PAIRED_END); }

  salmon::bam_utils::AlignerDetails getAlignerType() const { return aligner_; }

  // TODO: Make same as mapping-based
  void updateTranscriptLengthsAtomic(std::atomic<bool>& done) {
    if (sl_.try_lock()) {
      if (!done) {

        auto fld = flDist_.get();
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

  const std::vector<double>& condMeans() const { return conditionalMeans_; }

  std::vector<Transcript>& transcripts() { return transcripts_; }
  const std::vector<Transcript>& transcripts() const { return transcripts_; }

  inline bool getAlignmentGroup(AlignmentGroup<FragT>*& ag) {
    return bq->getAlignmentGroup(ag);
  }

  // inline t_pool* threadPool() { return threadPool_.get(); }

  inline SAM_hdr* header() { return bq->header(); }

  inline std::vector<FragmentStartPositionDistribution>&
  fragmentStartPositionDistributions() {
    return fragStartDists_;
  }

  inline FragmentLengthDistribution* fragmentLengthDistribution() const {
    return flDist_.get();
  }

  inline AlignmentModel& alignmentModel() { return *alnMod_.get(); }

  // SequenceBiasModel& sequenceBiasModel() { return seqBiasModel_; }

  //    inline tbb::concurrent_queue<FragT*>& fragmentQueue() {
  inline tbb::concurrent_queue<FragT*>& fragmentQueue() {
    return bq->getFragmentQueue();
  }

  //    inline tbb::concurrent_bounded_queue<AlignmentGroup<FragT*>*>&
  //    alignmentGroupQueue() {
  inline moodycamel::ConcurrentQueue<AlignmentGroup<FragT*>*>&
  alignmentGroupQueue() {
    return bq->getAlignmentGroupQueue();
  }

  inline BAMQueue<FragT>& getAlignmentGroupQueue() { return *bq.get(); }

  inline size_t upperBoundHits() { return bq->numMappedFragments(); }
  inline size_t numObservedFragments() const {
    return bq->numObservedFragments();
  }
  inline size_t numMappedFragments() const { return bq->numMappedFragments(); }
  inline size_t numUniquelyMappedFragments() {
    return bq->numUniquelyMappedFragments();
  }
  inline double effectiveMappingRate() const {
    return static_cast<double>(numMappedFragments()) / numObservedFragments();
  }

  // const boost::filesystem::path& alignmentFile() { return alignmentFile_; }

  ClusterForest& clusterForest() { return *clusters_.get(); }

  template <typename FilterT>
  bool reset(bool incPasses = true, FilterT filter = nullptr,
             bool onlyProcessAmbiguousAlignments = false) {
    namespace bfs = boost::filesystem;

    for (auto& alignmentFile : alignmentFiles_) {
      if (!bfs::is_regular_file(alignmentFile)) {
        return false;
      }
    }

    bq->reset();
    bq->start(filter, onlyProcessAmbiguousAlignments);
    if (incPasses) {
      quantificationPasses_++;
      fmt::print(stderr, "Current iteration = {}\n", quantificationPasses_);
    }
    return true;
  }

  inline LibraryFormat format() { return libFmt_; }
  inline const LibraryFormat format() const { return libFmt_; }

  /**
   * If this is set, attempt to automatically detect this library's type
   */
  void enableAutodetect() {
    // if auto detection is not already enabled, and we're enabling it
    if (!detector_) {
      detector_.reset(new LibraryTypeDetector(libFmt_.type));
    }
  }

  bool autoDetect() const { return (detector_.get() != nullptr); }

  LibraryTypeDetector* getDetector() { return detector_.get(); }

  LibraryFormat& getFormat() { return libFmt_; }
  const LibraryFormat& getFormat() const { return libFmt_; }

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
    return salmon::utils::DuplicateTargetStatus::UNKNOWN; 
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
  std::vector<boost::filesystem::path> alignmentFiles_;
  /**
   * The file from which the transcripts are read.
   * This is expected to be a FASTA format file.
   */
  boost::filesystem::path transcriptFile_;
  /**
   * Describes the expected format of the sequencing
   * fragment library.
   */
  LibraryFormat libFmt_;
  /**
   * The targets (transcripts) to be quantified.
   */
  std::vector<Transcript> transcripts_;
  /**
   * A pointer to the queue from which the fragments
   * will be read.
   */
  // std::unique_ptr<t_pool, std::function<void(t_pool*)>> threadPool_;
  std::unique_ptr<BAMQueue<FragT>> bq;

  // SequenceBiasModel seqBiasModel_;

  /**
   * The cluster forest maintains the dynamic relationship
   * defined by transcripts and reads --- if two transcripts
   * share an ambiguously mapped read, then they are placed
   * in the same cluster.
   */
  std::unique_ptr<ClusterForest> clusters_;

  /**
   * The emperical fragment start position distribution
   */
  std::vector<FragmentStartPositionDistribution> fragStartDists_;

  /**
   * The emperical fragment length distribution.
   *
   */
  std::unique_ptr<FragmentLengthDistribution> flDist_;
  /**
   * The emperical error model
   */
  std::unique_ptr<AlignmentModel> alnMod_;

  /** Keeps track of the number of passes that have been
   *  made through the alignment file.
   */
  size_t quantificationPasses_;
  SpinLock sl_;
  EQBuilderT eqBuilder_;

  /** Positional bias things**/
  std::vector<uint32_t> lengthQuantiles_;
  std::vector<SimplePosBias> posBiasFW_;
  std::vector<SimplePosBias> posBiasRC_;
  std::vector<SimplePosBias> posBiasExpectFW_;
  std::vector<SimplePosBias> posBiasExpectRC_;

  /** GC-fragment bias things **/
  // One bin for each percentage GC content
  double gcFracFwd_;
  GCFragModel observedGC_;
  GCFragModel expectedGC_;

  // Since multiple threads can touch this dist, we
  // need atomic counters.
  std::array<ReadKmerDist<6, std::atomic<uint32_t>>, 2> readBias_;
  std::array<SBModel, 2> readBiasModelObserved_;
  std::array<SBModel, 2> readBiasModelExpected_;

  // ReadKmerDist<6, std::atomic<uint32_t>> readBias_;
  std::vector<double> expectedBias_;
  std::unique_ptr<LibraryTypeDetector> detector_{nullptr};
  std::vector<double> conditionalMeans_;

  salmon::bam_utils::AlignerDetails aligner_{salmon::bam_utils::AlignerDetails::UNKNOWN};
  uint64_t numDecoys_{0};
};

#endif // ALIGNMENT_LIBRARY_HPP
