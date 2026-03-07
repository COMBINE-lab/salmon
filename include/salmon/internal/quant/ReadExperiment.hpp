#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

// Our includes
#include "salmon/internal/quant/ClusterForest.hpp"
#include "salmon/internal/util/DistributionUtils.hpp"
#include "EquivalenceClassBuilder.hpp"
#include "salmon/internal/model/FragmentLengthDistribution.hpp"
#include "salmon/internal/model/FragmentStartPositionDistribution.hpp"
#include "salmon/internal/model/GCFragModel.hpp"
#include "salmon/internal/model/ReadKmerDist.hpp"
#include "salmon/internal/quant/ReadLibrary.hpp"
#include "salmon/internal/quant/BiasLibraryState.hpp"
#include "salmon/internal/model/SBModel.hpp"
#include "salmon/internal/index/SalmonIndex.hpp"
#include "SalmonOpts.hpp"
#include "salmon/internal/util/SalmonUtils.hpp"
#include "salmon/internal/util/FmtCompat.hpp"
// #include "salmon/internal/model/SequenceBiasModel.hpp"
#include "salmon/internal/model/SimplePosBias.hpp"
#include "SpinLock.hpp" // RapMap's with try_lock
#include "salmon/internal/model/Transcript.hpp"
#include "salmon/internal/util/UtilityFunctions.hpp"

// Logger includes
#include <spdlog/spdlog.h>

// Boost includes
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
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
                 SalmonIndex* salmonIndex,
                 // const boost::filesystem::path& indexDirectory,
                 SalmonOpts& sopt);
  
  EQBuilderT& equivalenceClassBuilder() { return eqBuilder_; }

  std::string getIndexSeqHash256() const { return salmonIndex_->seqHash256(); }
  std::string getIndexNameHash256() const { return salmonIndex_->nameHash256(); }
  std::string getIndexSeqHash512() const { return salmonIndex_->seqHash512(); }
  std::string getIndexNameHash512() const { return salmonIndex_->nameHash512(); }
  std::string getIndexDecoySeqHash256() const { return salmonIndex_->decoySeqHash256(); }
  std::string getIndexDecoyNameHash256() const { return salmonIndex_->decoyNameHash256(); }

  std::vector<Transcript>& transcripts() { return transcripts_; }
  const std::vector<Transcript>& transcripts() const { return transcripts_; }

  const std::vector<double>& condMeans() const { return state_.conditionalMeans; }

  void updateTranscriptLengthsAtomic(std::atomic<bool>& done);

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

  SalmonIndex* getIndex() { return salmonIndex_; }

  template <typename PuffIndexT>
  void loadTranscriptsFromPuff(PuffIndexT* idx_, const SalmonOpts& sopt);

  template <typename QuasiIndexT>
  void loadTranscriptsFromQuasi(QuasiIndexT* idx_, const SalmonOpts& sopt);

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
    for (auto& rl : state_.library) {
      processReadLibrary(rl, salmonIndex_, transcripts_, clusterForest(),
                         *(fragLengthDist_.get()), numAssignedFragments_,
                         numThreads, burnedIn);
    }
    return true;
  }

  ~ReadExperiment() {
    // ---- Get rid of things we no longer need --------
  }

  ClusterForest& clusterForest() { return *clusters_.get(); }

  std::string readFilesAsString();

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
    return state_.fragmentStartDists;
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
    for (auto& rl : state_.library) {
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

  void summarizeLibraryTypeCounts(boost::filesystem::path& opath);

  std::vector<ReadLibrary>& readLibraries() { return state_.library; }
  const std::vector<ReadLibrary>& readLibraries() const {
    return state_.library;
  }
  FragmentLengthDistribution* fragmentLengthDistribution() const {
    return fragLengthDist_.get();
  }

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
    return salmonIndex_->index_retains_duplicates(); 
  }

private:
  void setTranscriptLengthClasses_(std::vector<uint32_t>& lengths,
                                   size_t nbins);

  /**
   * The file from which the alignments will be read.
   * This can be a SAM or BAM file, and can be a regular
   * file or a fifo.
   */
  salmon::quant::BiasLibraryState<std::vector<ReadLibrary>> state_;
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
  SalmonIndex* salmonIndex_{nullptr};
  /**
   * The cluster forest maintains the dynamic relationship
   * defined by transcripts and reads --- if two transcripts
   * share an ambiguously mapped read, then they are placed
   * in the same cluster.
   */
  std::unique_ptr<ClusterForest> clusters_;
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

  std::vector<uint32_t> lengthQuantiles_;

  uint64_t numDecoys_{0};
};

#include "salmon/internal/quant/ReadExperiment.inl"

extern template class ReadExperiment<EquivalenceClassBuilder<TGValue>>;
extern template class ReadExperiment<EquivalenceClassBuilder<SCTGValue>>;

#endif // EXPERIMENT_HPP
