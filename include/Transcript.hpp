#ifndef TRANSCRIPT
#define TRANSCRIPT

#include "FragmentLengthDistribution.hpp"
#include "GCFragModel.hpp"
#include "SalmonMath.hpp"
#include "SalmonStringUtils.hpp"
#include "SalmonUtils.hpp"
#include "SequenceBiasModel.hpp"
#include "stx/string_view.hpp"
#include "IOUtils.hpp"
#include <atomic>
#include <cmath>
#include <limits>
#include <memory>

//#include "rapmap/bit_array.h"
#include "pufferfish/compact_vector/compact_vector.hpp"
#include "pufferfish/rank9b.hpp"

class Transcript {
  static constexpr const uint32_t adapterBindingLength{5};
public:
  /*
  struct BitArrayDeleter {
    void operator()(BIT_ARRAY* b) {
      if (b != nullptr) {
        bit_array_free(b);
      }
    }
  };
  */
  using BitArray = compact::vector<uint64_t,1>;
  //using BitArrayPointer = std::unique_ptr<BIT_ARRAY, BitArrayDeleter>;
  using Rank9bPointer = std::unique_ptr<rank9b>;

  Transcript()
      : RefName(nullptr), RefLength(std::numeric_limits<uint32_t>::max()),
        CompleteLength(std::numeric_limits<uint32_t>::max()),
        EffectiveLength(-1.0), id(std::numeric_limits<uint32_t>::max()),
        logPerBasePrior_(salmon::math::LOG_0), priorMass_(salmon::math::LOG_0),
        mass_(salmon::math::LOG_0), sharedCount_(0.0),
        avgMassBias_(salmon::math::LOG_0), active_(false), skipBiasCorrection_(false) {
    uniqueCount_.store(0);
    totalCount_.store(0); // thanks @come-raczy
    cachedEffectiveLength_.store(salmon::math::LOG_0);
  }

  Transcript(size_t idIn, const char* name, uint32_t len, double alpha = 0.05)
      : RefName(name), RefLength(len), CompleteLength(len),
        EffectiveLength(-1.0), id(idIn), logPerBasePrior_(std::log(alpha)),
        priorMass_(std::log(alpha * len)), mass_(salmon::math::LOG_0),
        sharedCount_(0.0), avgMassBias_(salmon::math::LOG_0), active_(false), skipBiasCorrection_(false) {
    uniqueCount_.store(0);
    totalCount_.store(0); // thanks @come-raczy
    cachedEffectiveLength_.store(std::log(static_cast<double>(RefLength)));
  }

  Transcript(size_t idIn, const char* name, double len, bool /*updateEffLength*/, double alpha = 0.05)
    : RefName(name), RefLength(len), CompleteLength(len),
      EffectiveLength(len), id(idIn), logPerBasePrior_(std::log(alpha)),
      priorMass_(std::log(alpha * len)), mass_(salmon::math::LOG_0),
      sharedCount_(0.0), avgMassBias_(salmon::math::LOG_0), active_(false), skipBiasCorrection_(false) {
    uniqueCount_.store(0);
    totalCount_.store(0); // thanks @come-raczy
    cachedEffectiveLength_.store(std::log(len));
  }

  // We cannot copy; only move
  Transcript(Transcript& other) = delete;
  Transcript& operator=(Transcript& other) = delete;

  Transcript(Transcript&& other) {
    id = other.id;

    RefName = std::move(other.RefName);
    RefLength = other.RefLength;
    CompleteLength = other.CompleteLength;
    EffectiveLength = other.EffectiveLength;

    SAMSequence_ = std::move(other.SAMSequence_);
    Sequence_ = std::move(other.Sequence_);
    GCCount_ = std::move(other.GCCount_);
    reduceGCMemory_ = other.reduceGCMemory_;
    gcFracLen_ = other.gcFracLen_;
    lastRegularSample_ = other.lastRegularSample_;
    gcBitArray_ = std::move(other.gcBitArray_);
    gcRank_ = std::move(other.gcRank_);

    uniqueCount_.store(other.uniqueCount_);
    totalCount_.store(other.totalCount_.load());
    sharedCount_.store(other.sharedCount_.load());
    mass_.store(other.mass_.load());
    cachedEffectiveLength_.store(other.cachedEffectiveLength_.load());
    lengthClassIndex_ = other.lengthClassIndex_;
    logPerBasePrior_ = other.logPerBasePrior_;
    priorMass_ = other.priorMass_;
    avgMassBias_.store(other.avgMassBias_.load());
    hasAnchorFragment_.store(other.hasAnchorFragment_.load());
    active_ = other.active_;
    isDecoy_ = other.isDecoy_;
    skipBiasCorrection_ = other.skipBiasCorrection_;
  }

  Transcript& operator=(Transcript&& other) {
    id = other.id;

    RefName = std::move(other.RefName);
    RefLength = other.RefLength;
    CompleteLength = other.CompleteLength;
    EffectiveLength = other.EffectiveLength;
    SAMSequence_ = std::move(other.SAMSequence_);
    Sequence_ = std::move(other.Sequence_);
    GCCount_ = std::move(other.GCCount_);
    reduceGCMemory_ = other.reduceGCMemory_;
    gcFracLen_ = other.gcFracLen_;
    lastRegularSample_ = other.lastRegularSample_;
    gcBitArray_ = std::move(other.gcBitArray_);
    gcRank_ = std::move(other.gcRank_);

    uniqueCount_.store(other.uniqueCount_);
    totalCount_.store(other.totalCount_.load());
    sharedCount_.store(other.sharedCount_.load());
    mass_.store(other.mass_.load());
    cachedEffectiveLength_.store(other.cachedEffectiveLength_.load());
    lengthClassIndex_ = other.lengthClassIndex_;
    logPerBasePrior_ = other.logPerBasePrior_;
    priorMass_ = other.priorMass_;
    avgMassBias_.store(other.avgMassBias_.load());
    hasAnchorFragment_.store(other.hasAnchorFragment_.load());
    active_ = other.active_;
    isDecoy_ = other.isDecoy_;
    skipBiasCorrection_ = other.skipBiasCorrection_;
    return *this;
  }

  inline double sharedCount() { return sharedCount_.load(); }
  inline size_t uniqueCount() { return uniqueCount_.load(); }
  inline size_t totalCount() { return totalCount_.load(); }

  inline void addUniqueCount(size_t newCount) { uniqueCount_ += newCount; }
  inline void addTotalCount(size_t newCount) { totalCount_ += newCount; }

  inline double uniqueUpdateFraction() const {
    double ambigCount = static_cast<double>(totalCount_ - uniqueCount_);
    return uniqueCount_ / ambigCount;
  }

  inline char charBaseAt(size_t idx, salmon::stringtools::strand dir =
                                         salmon::stringtools::strand::forward) {
    return salmon::stringtools::samCodeToChar[baseAt(idx, dir)];
  }

  inline uint8_t baseAt(size_t idx, salmon::stringtools::strand dir =
                                        salmon::stringtools::strand::forward) {
    using salmon::stringtools::strand;
    using salmon::stringtools::encodedRevComp;
    size_t byte = idx >> 1;
    size_t nibble = idx & 0x1;

    // NOTE 10.2
    auto& sseq = SAMSequence_;
    //if (byte >= sseq.size()) { std::cerr << "requested index " << byte << " for vector of size " << sseq.size() << " for reference " << RefName << std::endl; return 0;}


    switch (dir) {
    case strand::forward:
      if (nibble) {
        return sseq[byte] & 0x0F;
      } else {
        return ((sseq[byte] & 0xF0) >> 4) & 0x0F;
      }
      break;
    case strand::reverse:
      if (nibble) {
        return encodedRevComp[sseq[byte] & 0x0F];
      } else {
        return encodedRevComp[((sseq[byte] & 0xF0) >> 4) & 0x0F];
      }
      break;
    }

    return std::numeric_limits<uint8_t>::max();
  }

  inline void setSharedCount(double sc) { sharedCount_.store(sc); }

  inline void addSharedCount(double sc) {
    salmon::utils::incLoop(sharedCount_, sc);
  }

  inline void addBias(double bias) {
    salmon::utils::incLoopLog(avgMassBias_, bias);
  }

  inline void addMass(double mass) { salmon::utils::incLoopLog(mass_, mass); }

  inline void setMass(double mass) { mass_.store(mass); }

  inline double mass(bool withPrior = true) {
    return (withPrior) ? salmon::math::logAdd(priorMass_, mass_.load())
                       : mass_.load();
  }

  void setActive() { active_ = true; }
  bool getActive() { return active_; }

  inline double bias() {
    return (totalCount_.load() > 0)
               ? avgMassBias_ - std::log(totalCount_.load())
               : salmon::math::LOG_1;
  }

  /*
  double getAverageSequenceBias(SequenceBiasModel& m) {
      double bias = salmon::math::LOG_0;
      for (int32_t i = 0; i < RefLength; ++i) {
          bias = salmon::math::logAdd(bias, m.biasFactor(*this, i));
      }
      return bias - std::log(RefLength);
  }
  */

  /**
   *  NOTE: Adopted from "est_effective_length" at
   * (https://github.com/adarob/eXpress/blob/master/src/targets.cpp) originally
   * written by Adam Roberts.
   *
   *
   */
  double computeLogEffectiveLength(std::vector<double>& logPMF,
                                   size_t minVal,
                                   size_t maxVal) {

    double effectiveLength = salmon::math::LOG_0;
    double refLen = static_cast<double>(RefLength);
    double logRefLength = std::log(refLen);

    uint32_t mval = maxVal;
    size_t clen = minVal;
    size_t maxLen = std::min(RefLength, mval);
    while (clen <= maxLen) {
      size_t i = clen - minVal;
      effectiveLength = salmon::math::logAdd(
          effectiveLength, logPMF[i] + std::log(refLen - clen + 1));
      ++clen;
    }
    if (salmon::math::isLog0(effectiveLength) or
        std::exp(effectiveLength) < 1.0) {
      effectiveLength = logRefLength;
    }

    return effectiveLength;
  }

  /**
   * Return the cached value for the log of the effective length.
   */
  double getCachedLogEffectiveLength() { return cachedEffectiveLength_.load(); }

  void setCachedLogEffectiveLength(double l) {
    cachedEffectiveLength_.store(l);
  }

  void updateEffectiveLength(std::vector<double>& logPMF, double /*logFLDMean*/,
                             size_t minVal, size_t maxVal) {
    double cel = computeLogEffectiveLength(logPMF, minVal, maxVal);
    cachedEffectiveLength_.store(cel);
  }

  double perBasePrior() { return std::exp(logPerBasePrior_); }

  void lengthClassIndex(uint32_t ind) { lengthClassIndex_ = ind; }
  uint32_t lengthClassIndex() const { return lengthClassIndex_; }

  void setAnchorFragment() { hasAnchorFragment_.store(true); }

  bool hasAnchorFragment() { return hasAnchorFragment_.load(); }

  inline GCDesc gcDesc(int32_t s, int32_t e, bool& valid) const {
    int outsideContext{3};
    int insideContext{2};

    int outside5p = outsideContext + 1;
    int outside3p = outsideContext;

    int inside5p = insideContext - 1;
    int inside3p = insideContext;

    double contextSize = outsideContext + insideContext;
    int lastPos = RefLength - 1;
    if (!reduceGCMemory_) {
      auto cs = (s > 0) ? GCCount_[s - 1] : 0;
      auto ce = GCCount_[e];

      int fs = s - outside5p;
      int fe = s + inside5p;
      int ts = e - inside3p;
      int te = e + outside3p;

      bool fpLeftExists = (fs >= 0);
      bool fpRightExists = (fe <= lastPos);
      bool tpLeftExists = (ts >= 0);
      bool tpRightExists = (te <= lastPos);

      auto fps = (fpLeftExists) ? GCCount_[fs] : 0;
      auto fpe = (fpRightExists) ? GCCount_[fe] : ce;
      auto tps = (tpLeftExists) ? GCCount_[ts] : 0;
      auto tpe = (tpRightExists) ? GCCount_[te] : ce;

      // now, clamp to actual bounds
      fs = (fs < 0) ? 0 : fs;
      fe = (fe > lastPos) ? lastPos : fe;
      ts = (ts < 0) ? 0 : ts;
      te = (te > lastPos) ? lastPos : te;
      int fpContextSize = (!fpLeftExists) ? (fe + 1) : (fe - fs);
      int tpContextSize = (!tpLeftExists) ? (te + 1) : (te - ts);
      contextSize = static_cast<double>(fpContextSize + tpContextSize);
      if (contextSize == 0) {
        // std::cerr << "" << std::endl;
        return GCDesc();
      }
      valid = true;

      int32_t fragFrac = std::lrint((100.0 * (ce - cs)) / (e - s + 1));
      int32_t contextFrac =
          std::lrint((100.0 * (((fpe - fps) + (tpe - tps)) / (contextSize))));
      /*
      if (contextFrac > 100) {
        std::cerr << "NOTE : 5' count = " << (fpeCount - fpsCount) << ", 3'
      count =" << (tpeCount - tpsCount) << ", context size = " << contextSize <<
      std::endl; std::cerr << "s = " << s << ", e = " << e << ", l = " <<
      RefLength << ", fs = " << fs  << ", fe =  " << fe << ", ts = " << ts << ",
      te = " << te << std::endl; std::cerr << "fpsCount = " << fpsCount <<
          ", fpeCount = " << fpeCount <<
          ", tpsCount = " << tpsCount <<
          ", tpeCount = " << tpeCount  <<
          ", fpContextSize = " << fpContextSize  <<
          ", tpContextSize = " << tpContextSize << std::endl;

      }
      */
      GCDesc desc = {fragFrac, contextFrac};
      return desc;
    } else {
      auto cs = (s > 0) ? gcCountInterp_(s - 1) : 0;
      auto ce = gcCountInterp_(e);

      int fs = s - outside5p;
      int fe = s + inside5p;
      int ts = e - inside3p;
      int te = e + outside3p;

      bool fpLeftExists = (fs >= 0);
      bool fpRightExists = (fe <= lastPos);
      bool tpLeftExists = (ts >= 0);
      bool tpRightExists = (te <= lastPos);

      auto fps = (fpLeftExists) ? gcCountInterp_(fs) : 0;
      auto fpe = (fpRightExists) ? gcCountInterp_(fe) : ce;
      auto tps = (tpLeftExists) ? gcCountInterp_(ts) : 0;
      auto tpe = (tpRightExists) ? gcCountInterp_(te) : ce;

      // now, clamp to actual bounds
      fs = (fs < 0) ? 0 : fs;
      fe = (fe > lastPos) ? lastPos : fe;
      ts = (ts < 0) ? 0 : ts;
      te = (te > lastPos) ? lastPos : te;
      int fpContextSize = (!fpLeftExists) ? (fe + 1) : (fe - fs);
      int tpContextSize = (!tpLeftExists) ? (te + 1) : (te - ts);
      contextSize = static_cast<double>(fpContextSize + tpContextSize);
      if (contextSize == 0) {
        return GCDesc();
      }
      valid = true;

      int32_t fragFrac = std::lrint((100.0 * (ce - cs)) / (e - s + 1));
      int32_t contextFrac =
          std::lrint((100.0 * (((fpe - fps) + (tpe - tps)) / (contextSize))));
      GCDesc desc = {fragFrac, contextFrac};
      return desc;
      /** PREVIOUS SAMPLED IMPL (before May 23, 2017) **/
      /*
        auto cs = gcCountInterp_(s);
        auto ce = gcCountInterp_(e);

        valid = true;
        auto fps = (s >= outside5p) ? gcCountInterp_(s-outside5p) : 0;
        auto fpe = (inside5p > 0) ? gcCountInterp_(std::min(s+inside5p,
        lastPos)) : cs; auto tps = (inside3p > 0) ?
        ((e >= inside3p) ? gcCountInterp_(e-inside3p) : 0) : ce;
        auto tpe = gcCountInterp_(std::min(e+outside3p, lastPos));

        int32_t fragFrac = std::lrint((100.0 * (ce - cs)) / (e - s + 1));
        int32_t contextFrac = std::lrint((100.0 * (((fpe - fps) + (tpe - tps)) /
        (2.0 * contextSize)))); GCDesc desc = {fragFrac, contextFrac}; return
        desc;
      */
    }
  }
  inline double gcAt(int32_t s) const {
    int32_t sRefLength = static_cast<int32_t>(RefLength);
    return (s < 0) ? 0.0
                   : ((s >= sRefLength) ? gcCount_(sRefLength - 1) : gcCount_(s));
  }

  // Return the fractional GC content along this transcript
  // in the interval [s,e] (note; this interval is closed on both sides).
  inline int32_t gcFrac(int32_t s, int32_t e) const {
    if (!reduceGCMemory_) {
      auto cs = (s > 0) ? GCCount_[s - 1] : 0;
      auto ce = GCCount_[e];
      return std::lrint((100.0 * (ce - cs)) / (e - s + 1));
    } else {
      auto cs = (s > 0) ? gcCountInterp_(s - 1) : 0;
      auto ce = gcCountInterp_(e);
      return std::lrint((100.0 * (ce - cs)) / (e - s + 1));
    }
    return 0;
  }

  /**
   * Return the next polyA site that occurs in this transcript
   * after position p
   **/
  /*
  inline int32_t getNextPolyA(int32_t p) {
    if (p+1 >= static_cast<int32_t>(RefLength)) { return RefLength; }
    auto r = polyARank_->rank(p+1);
    return polyAPos_[r];
  }
  */

  void setDecoy(bool isDecoy) {
    isDecoy_ = isDecoy;
  }

  bool isDecoy() const { return isDecoy_; }

  // Will *not* delete seq on destruction
  void setSequenceBorrowed(const char* seq, bool needGC = false,
                           bool reduceGCMemory = false) {
    Sequence_ = std::unique_ptr<const char, void (*)(const char*)>(
        seq,                 // store seq
        [](const char* /*p*/) {} // do nothing deleter
    );
    if (needGC) {
      computeGCContent_(reduceGCMemory);
    }
  }

  // Will delete seq on destruction
  void setSequenceOwned(const char* seq, bool needGC = false,
                        bool reduceGCMemory = false) {
    Sequence_ = std::unique_ptr<const char, void (*)(const char*)>(
        seq,                              // store seq
        [](const char* p) { delete[] p; } // do nothing deleter
    );
    if (needGC) {
      computeGCContent_(reduceGCMemory);
    }
  }

  // Will delete seq on destruction
  void setSAMSequenceOwned(std::vector<uint8_t>&& seq, bool needGC = false,
                           bool reduceGCMemory = false) {

    if ((2*seq.size() < RefLength) or (2*seq.size() > RefLength + 1)) {
      std::stringstream errstream;
      errstream << "\n\nSAM file says target " << RefName << " has length " << RefLength
                << ", but the FASTA file contains a sequence of length [" << seq.size() * 2 << " or " << seq.size() * 2 - 1 << "]\n\n";
      std::cerr << ioutils::SET_RED << errstream.str();
      std::exit(1);
    }

    SAMSequence_ = std::move(seq);
    if (needGC) {
      computeGCContent_(reduceGCMemory);
    }
  }

  bool have_sequence() const { if (Sequence_) { return true; } else { return false; } }
  
  const char* Sequence() const { return Sequence_.get(); }

  uint8_t* SAMSequence() const { return const_cast<uint8_t*>(SAMSequence_.data()); }

  void setCompleteLength(uint32_t completeLengthIn) {
    CompleteLength = completeLengthIn;
  }

  void computePolyAPositions() { computePolyAPositions_(); }

  void setSkipBiasCorrection(bool skip) { skipBiasCorrection_ = skip; }
  bool skipBiasCorrection() const { return skipBiasCorrection_; }

  std::string RefName;
  uint32_t RefLength;
  uint32_t CompleteLength;
  double EffectiveLength;
  uint32_t id;

  double uniqueCounts{0.0};
  double totalCounts{0.0};
  double projectedCounts{0.0};
  double sharedCounts{0.0};

private:
  // NOTE: Is it worth it to check if we have GC here?
  // we should never access these without bias correction.
  inline double gcCount_(int32_t p) {
    return (!reduceGCMemory_) ? static_cast<double>(GCCount_[p])
                              : gcCountInterp_(p);
  }
  inline double gcCount_(int32_t p) const {
    return (!reduceGCMemory_) ? static_cast<double>(GCCount_[p])
                              : gcCountInterp_(p);
  }

  /*
    inline int32_t closestBin_(int32_t p) const {
      return static_cast<int32_t>(std::round( static_cast<double>(p) / gcStep_
    ));
    }
  */

  inline double gcCountInterp_(int32_t p) const {
    int32_t sRefLength = static_cast<int32_t>(RefLength);
    uint64_t extraCount = 0;
    if ((p+1) >= sRefLength) {
      extraCount = gcBitArray_[(RefLength-1)];
      p = sRefLength - 2;
    }
    return static_cast<double>(gcRank_->rank(p + 1) + extraCount);
  }

  /** Previous GC count interp implementation (May 23, 2017) **/
  /*
    inline double gcCountInterp_(int32_t p) const {
      //std::cerr << "in gcCountInterp\n";
      if (p == RefLength - 1) {
        // If p is the last position, just return the last value
        return static_cast<double>(GCCount_.back());
      }

      // The index of the closest bin
      auto cb = closestBin_(p);
      // The actual position to which this bin corresponds
      int32_t binPos = cb * gcStep_;
      // Can't go past the end
      if (binPos > RefLength - 1) {
        binPos = RefLength - 1;
        cb = GCCount_.size() - 1;
      }

      // The count of {G,C} at the checkpoint
      auto binCount = GCCount_[cb];
      // The count before or after the bin, until p
      int32_t count{0};
      const char* seq = Sequence_.get();

      // we hit a sampled position
      if (binPos == p) {
      } else if (binPos > p) {
        for (size_t i = binPos; i > p; --i) {
          auto c = seq[i];
          // If the character is a G or C, we subtract 1
          count -= (c == 'G' or c == 'C') ? 1 : 0;
        }
      } else {
        for (size_t i = binPos + 1; i <= p; ++i) {
          auto c = seq[i];
          // If the character is a G or C, we add 1
          count += (c == 'G' or c == 'C') ? 1 : 0;
        }
      }
      return  binCount + count;
    }
  */

  /*
    void computeGCContentSampled_(uint32_t step) {
        gcStep_ = step;
        const char* seq = Sequence_.get();
        size_t nsamp = std::ceil(static_cast<double>(RefLength) / step);
        GCCount_.reserve(nsamp + 2);

        size_t lastSamp{0};
        size_t totGC{0};
        for (size_t i = 0; i < RefLength; ++i) {
            auto c = std::toupper(seq[i]);
            if (c == 'G' or c == 'C') {
                totGC++;
            }
            if (i % step == 0) {
                GCCount_.push_back(totGC);
                lastSamp = i;
            }
        }

        if (lastSamp < RefLength - 1) {
            GCCount_.push_back(totGC);
        }

        gcFracLen_ = static_cast<double>(RefLength - 1) / gcStep_;
        lastRegularSample_ = std::ceil(gcFracLen_);
    }
  */

  void computeGCContent_(bool reduceGCMemory) {
    reduceGCMemory_ = reduceGCMemory;
    const char* seq = Sequence_.get();
    GCCount_.clear();
    if (!reduceGCMemory) {
      GCCount_.resize(RefLength, 0);
      size_t totGC{0};
      for (size_t i = 0; i < RefLength; ++i) {
        auto c = std::toupper(seq[i]);
        if (c == 'G' or c == 'C') {
          totGC++;
        }
        GCCount_[i] = totGC;
      }
    } else {
      gcBitArray_.resize(RefLength);
      gcBitArray_.clear_mem();
      //BIT_ARRAY* rawArray = bit_array_create(RefLength);
      for (size_t i = 0; i < RefLength; ++i) {
        auto c = std::toupper(seq[i]);
        if (c == 'G' or c == 'C') {
          gcBitArray_[i] = 1;
          //bit_array_set_bit(rawArray, i);
        }
      }
      //gcBitArray_.reset(rawArray);
      gcRank_.reset(new rank9b(gcBitArray_.get(), RefLength));
    }
  }

  void computePolyAPositions_() {
    /*
    polyAPos_.clear();
    BIT_ARRAY* rawArray = bit_array_create(RefLength);
    std::string polyA(adapterBindingLength, 'A');
    stx::string_view polyAView(polyA);
    stx::string_view seq(Sequence_.get(), RefLength);
    auto occIt = seq.find(polyAView);
    auto prev = occIt;
    auto end = stx::string_view::npos;
    while (occIt != end) {
      auto d = occIt;
      bit_array_set_bit(rawArray, d);
      polyAPos_.push_back(d);
      prev = occIt;
      occIt = seq.find(polyAView, d+polyAView.length());
      // if this is the same stretch of polyA, skip again
      if (occIt != end and (static_cast<int64_t>(occIt) - prev) < (polyAView.length() + 1)) {
        occIt = seq.find_first_not_of('A', d+polyAView.length());
        occIt = seq.find(polyAView, occIt);
      }
    }
    polyAPos_.push_back(RefLength);
    polyABitArray_.reset(rawArray);
    polyARank_.reset(new rank9b(polyABitArray_->words, RefLength));
    */
  }


  // NOTE 10.2
  std::vector<uint8_t> SAMSequence_;

  std::unique_ptr<const char, void (*)(const char*)> Sequence_ =
      std::unique_ptr<const char, void (*)(const char*)>(nullptr,
                                                         [](const char*) {});

  std::atomic<size_t> uniqueCount_;
  std::atomic<size_t> totalCount_;
  double priorMass_;
  std::atomic<double> mass_;
  std::atomic<double> sharedCount_;
  std::atomic<double> cachedEffectiveLength_;
  std::atomic<double> avgMassBias_;
  uint32_t lengthClassIndex_;
  double logPerBasePrior_;
  // In a paired-end protocol, a transcript has
  // an "anchor" fragment if it has a proper
  // pair of reads mapping to it.
  std::atomic<bool> hasAnchorFragment_{false};
  bool active_;
  bool isDecoy_{false};
  bool skipBiasCorrection_{false};

  bool reduceGCMemory_{false};
  double gcFracLen_{0.0};
  uint32_t lastRegularSample_{0};
  std::vector<uint32_t> GCCount_;
  //BitArrayPointer gcBitArray_{nullptr};
  BitArray gcBitArray_;
  Rank9bPointer gcRank_{nullptr};
  //BitArrayPointer polyABitArray_{nullptr};
  //Rank9bPointer polyARank_{nullptr};
  //std::vector<int32_t> polyAPos_;
};

#endif // TRANSCRIPT
