#ifndef ALIGNMENT_LIBRARY_HPP
#define ALIGNMENT_LIBRARY_HPP

#include "salmon/internal/io/AlignmentIO.hpp"

// Our includes
#include "salmon/internal/alignment/AlignmentGroup.hpp"
#include "salmon/internal/alignment/AlignmentModel.hpp"
#include "salmon/internal/alignment/ONTAlignmentModel.hpp"
#include "salmon/internal/alignment/BAMQueue.hpp"
#include "salmon/internal/alignment/BAMUtils.hpp"
#include "salmon/internal/alignment/NullFragmentFilter.hpp"
#include "salmon/internal/quant/ClusterForest.hpp"
#include "salmon/internal/util/DistributionUtils.hpp"
#include "EquivalenceClassBuilder.hpp"
#include "salmon/internal/io/FASTAParser.hpp"
#include "salmon/internal/model/FragmentLengthDistribution.hpp"
#include "salmon/internal/model/FragmentStartPositionDistribution.hpp"
#include "salmon/internal/model/GCFragModel.hpp"
#include "salmon/internal/model/LibraryFormat.hpp"
#include "salmon/internal/model/LibraryTypeDetector.hpp"
#include "salmon/internal/model/ReadKmerDist.hpp"
#include "salmon/internal/model/SBModel.hpp"
#include "salmon/internal/quant/BiasLibraryState.hpp"
#include "SalmonOpts.hpp"
#include "salmon/internal/util/SalmonUtils.hpp"
#include "salmon/internal/model/SimplePosBias.hpp"
#include "SpinLock.hpp" // From pufferfish, with try_lock
#include "salmon/internal/model/Transcript.hpp"
#include "salmon/internal/alignment/ReadPair.hpp"
#include "salmon/internal/alignment/UnpairedRead.hpp"
#include <concurrentqueue.h>
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
template <typename FragT, typename EQBuilderT, typename AlignModelT> class AlignmentLibrary {

public:
  AlignmentLibrary(std::vector<boost::filesystem::path>& alnFiles,
                   const boost::filesystem::path& transcriptFile,
                   LibraryFormat libFmt, SalmonOpts& salmonOpts);

  AlignmentLibrary(std::vector<boost::filesystem::path>& alnFiles,
                   LibraryFormat libFmt, SalmonOpts& salmonOpts,
                   bool /*eqClassMode_*/, std::vector<std::string>& tnames,
                   std::vector<double>& tefflens);

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
  bool isPairedEnd() { return (state_.library.type == ReadType::PAIRED_END); }

  salmon::bam_utils::AlignerDetails getAlignerType() const { return aligner_; }

  // TODO: Make same as mapping-based
  void updateTranscriptLengthsAtomic(std::atomic<bool>& done);

  const std::vector<double>& condMeans() const { return state_.conditionalMeans; }

  std::vector<Transcript>& transcripts() { return transcripts_; }
  const std::vector<Transcript>& transcripts() const { return transcripts_; }

  inline bool getAlignmentGroup(AlignmentGroup<FragT*>*& ag) {
    return bq->getAlignmentGroup(ag);
  }

  // inline t_pool* threadPool() { return threadPool_.get(); }

  inline SAM_hdr* header() { return bq->header(); }

  inline std::vector<FragmentStartPositionDistribution>&
  fragmentStartPositionDistributions() {
    return state_.fragmentStartDists;
  }

  inline FragmentLengthDistribution* fragmentLengthDistribution() const {
    return flDist_.get();
  }

  inline AlignModelT& alignmentModel() { return *alnMod_.get(); }

  // SequenceBiasModel& sequenceBiasModel() { return seqBiasModel_; }

  //    inline oneapi::tbb::concurrent_queue<FragT*>& fragmentQueue() {
  inline oneapi::tbb::concurrent_queue<FragT*>& fragmentQueue() {
    return bq->getFragmentQueue();
  }

  //    inline oneapi::tbb::concurrent_bounded_queue<AlignmentGroup<FragT*>*>&
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

  inline LibraryFormat format() { return state_.library; }
  inline const LibraryFormat format() const { return state_.library; }

  /**
   * If this is set, attempt to automatically detect this library's type
   */
  void enableAutodetect() {
    // if auto detection is not already enabled, and we're enabling it
    if (!state_.detector) {
      state_.detector.reset(new LibraryTypeDetector(state_.library.type));
    }
  }

  bool autoDetect() const { return (state_.detector.get() != nullptr); }

  LibraryTypeDetector* getDetector() { return state_.detector.get(); }

  LibraryFormat& getFormat() { return state_.library; }
  const LibraryFormat& getFormat() const { return state_.library; }

  void setGCFracForward(double fracForward) { state_.gcFracFwd = fracForward; }

  double gcFracFwd() const { return state_.gcFracFwd; }
  double gcFracRC() const { return 1.0 - state_.gcFracFwd; }

  std::vector<double>& expectedSeqBias() { return state_.expectedSeqBias; }

  const std::vector<double>& expectedSeqBias() const { return state_.expectedSeqBias; }

  void setExpectedGCBias(const GCFragModel& expectedBiasIn) {
    state_.expectedGC = expectedBiasIn;
  }

  GCFragModel& expectedGCBias() { return state_.expectedGC; }

  const GCFragModel& expectedGCBias() const { return state_.expectedGC; }

  const GCFragModel& observedGC() const { return state_.observedGC; }

  GCFragModel& observedGC() { return state_.observedGC; }

  std::vector<SimplePosBias>& posBias(salmon::utils::Direction dir) {
    return (dir == salmon::utils::Direction::FORWARD) ? state_.posBiasFW : state_.posBiasRC;
  }
  const std::vector<SimplePosBias>&
  posBias(salmon::utils::Direction dir) const {
    return (dir == salmon::utils::Direction::FORWARD) ? state_.posBiasFW : state_.posBiasRC;
  }

  std::vector<SimplePosBias>& posBiasExpected(salmon::utils::Direction dir) {
    return (dir == salmon::utils::Direction::FORWARD) ? state_.posBiasExpectFW
                                                      : state_.posBiasExpectRC;
  }

  const std::vector<SimplePosBias>&
  posBiasExpected(salmon::utils::Direction dir) const {
    return (dir == salmon::utils::Direction::FORWARD) ? state_.posBiasExpectFW
                                                      : state_.posBiasExpectRC;
  }

  ReadKmerDist<6, std::atomic<uint32_t>>&
  readBias(salmon::utils::Direction dir) {
    return (dir == salmon::utils::Direction::FORWARD) ? state_.readBias[0]
                                                      : state_.readBias[1];
  }
  const ReadKmerDist<6, std::atomic<uint32_t>>&
  readBias(salmon::utils::Direction dir) const {
    return (dir == salmon::utils::Direction::FORWARD) ? state_.readBias[0]
                                                      : state_.readBias[1];
  }

  SBModel& readBiasModelObserved(salmon::utils::Direction dir) {
    return (dir == salmon::utils::Direction::FORWARD)
               ? state_.readBiasModelObserved[0]
               : state_.readBiasModelObserved[1];
  }
  const SBModel& readBiasModelObserved(salmon::utils::Direction dir) const {
    return (dir == salmon::utils::Direction::FORWARD)
               ? state_.readBiasModelObserved[0]
               : state_.readBiasModelObserved[1];
  }

  SBModel& readBiasModelExpected(salmon::utils::Direction dir) {
    return (dir == salmon::utils::Direction::FORWARD)
               ? state_.readBiasModelExpected[0]
               : state_.readBiasModelExpected[1];
  }
  const SBModel& readBiasModelExpected(salmon::utils::Direction dir) const {
    return (dir == salmon::utils::Direction::FORWARD)
               ? state_.readBiasModelExpected[0]
               : state_.readBiasModelExpected[1];
  }

  void setReadBiasModelExpected(SBModel&& model, salmon::utils::Direction dir) {
    size_t idx = (dir == salmon::utils::Direction::FORWARD) ? 0 : 1;
    state_.readBiasModelExpected[idx] = std::move(model);
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
                                   size_t nbins);

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
  salmon::quant::BiasLibraryState<LibraryFormat> state_;
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
   * The emperical fragment length distribution.
   *
   */
  std::unique_ptr<FragmentLengthDistribution> flDist_;
  /**
   * The emperical error model
   */
  std::unique_ptr<AlignModelT> alnMod_;

  /** Keeps track of the number of passes that have been
   *  made through the alignment file.
   */
  size_t quantificationPasses_;
  SpinLock sl_;
  EQBuilderT eqBuilder_;

  std::vector<uint32_t> lengthQuantiles_;

  salmon::bam_utils::AlignerDetails aligner_{salmon::bam_utils::AlignerDetails::UNKNOWN};
  uint64_t numDecoys_{0};
};

#include "salmon/internal/alignment/AlignmentLibrary.inl"

extern template class AlignmentLibrary<ReadPair, EquivalenceClassBuilder<TGValue>, AlignmentModel>;
extern template class AlignmentLibrary<UnpairedRead, EquivalenceClassBuilder<TGValue>, AlignmentModel>;
extern template class AlignmentLibrary<UnpairedRead, EquivalenceClassBuilder<TGValue>, ONTAlignmentModel>;

#endif // ALIGNMENT_LIBRARY_HPP
