#ifndef __DISTRIBUTION_UTILS__
#define __DISTRIBUTION_UTILS__

#include <cstddef>
#include <cstdint>
#include <vector>

// Do things involving distributions
class FragmentLengthDistribution;
class Transcript;

namespace distribution_utils {
enum class DistributionSpace : uint8_t { LOG = 0, LINEAR = 1 };

class DistSummary {
public:
DistSummary(double mean_in, double sd_in, uint32_t support) :
  mean(mean_in), sd(sd_in), samples(support,0) {}

double mean;
double sd;
std::vector<int32_t> samples;
};

/**
 *  Draw samples from the provided fragment length distribution.
 *  \param fld A pointer to the FragmentLengthDistribution from which
 *             samples will be drawn.
 *  \param numSamples  The number of samples to draw.
 */
DistSummary samplesFromLogPMF(FragmentLengthDistribution* fld,
                                       int32_t numSamples);

/**
 * The following two functions compute conditional means of the empirical
 *fragment length distribution, and apply them to the transcripts to compute the
 *effective lengths. To the best of our knowledge, this particular approach to
 *effective length correction was first introduced in Kallisto[1]. [1] Bray,
 *Nicolas L., et al. "Near-optimal probabilistic RNA-seq quantification."
 *Nature biotechnology 34.5 (2016): 525-527.
 **/

/**
 * Compute the conditional mean fragment length for every length
 * in the input fragment length distribution.  For each length i,
 * the conditional mean assumes that it is not possible to sample fragments
 * of length > i, and so the probability mass is normalized by the
 * probability of all lengths <= i.
 *
 * \param mass The input fragment length distribution.  This should contain a
 * number for each fragment length that is proportional to the probability of
 *             drawing a fragment of that length.  The input need not be
 * normalized. \param inputSpace A DistributionSpace parameter that determines
 * whether mass should be interpreted as exisitng in log space or linear space.
 * \returns The conditional means for each fragment length.
 */
std::vector<double> correctionFactorsFromMass(std::vector<double>& mass,
                                              DistributionSpace inputSpace);

/**
 * Populate the effective lengths of the input transcripts based on the
 * conditional means.
 *
 * \sa correctionFactorsFromMass()
 * \param maxLength The maximum fragment length.
 * \param transcripts The transcripts whose lengths are to be corrected.
 * \param correctionFactors The conditional means (computed with
 * correctionFactorsFromMass())
 */
void computeSmoothedEffectiveLengths(size_t maxLength,
                                     std::vector<Transcript>& transcripts,
                                     std::vector<double>& correctionFactors,
                                     DistributionSpace outputSpace);

/**
 * Return the logCMF from a "snapshot" of the input distribution at the current
 * point in time.  This provides the ability to quickly recall approximate cumulative
 * probabilities quickly, versus having to compute them each time.
 *
 * \param fld A pointer to the fragment length distribution.
 */
std::vector<double> evaluateLogCMF(FragmentLengthDistribution* fld);

class LogCMFCache {
public:
  LogCMFCache(FragmentLengthDistribution* fld,
              bool singleEndLib, size_t fragUpdateThresh=100000);
  void refresh(size_t processedReads, bool burnedIn);
  double getAmbigFragLengthProb(bool fwd, int32_t pos, int32_t rlen, int32_t tlen, bool burnedIn) const;
private:
  inline double cmfValue_(size_t len, bool useFLD) const;
  FragmentLengthDistribution* fld_{nullptr};
  bool singleEndLib_{false};
  size_t fragUpdateThresh_{100000};
  size_t prevProcessedReads_{0};
  std::vector<double> cachedCMF_;
};

template <typename T>
class VersionedValue {
  public:
    T val{T()};
    uint64_t gen{0};
};

// A simple cache where values are retreived by their index
// and the cached items are versioned
template <typename T>
class IndexedVersionedCache {
  public:

  IndexedVersionedCache(size_t max_index) :
    cache_(max_index+1), max_index_(max_index), current_gen_(0) {}

  inline void increment_generation() { ++current_gen_; }

  inline bool get_value(size_t index, T& v) {
    size_t idx = (index > max_index_) ? max_index_ : index;
    const VersionedValue<T>& vv = cache_[idx];
    bool is_stale = vv.gen < current_gen_;
    v = vv.val;
    return is_stale;
  }

  inline void update_value(size_t index, T v) {
    size_t idx = (index > max_index_) ? max_index_ : index;
    VersionedValue<T>& vv = cache_[idx];
    vv.val = v;
    vv.gen = current_gen_;
  }

  template <typename F>
  inline T get_or_update(size_t index, F& gen_value) {
    size_t idx = (index > max_index_) ? max_index_ : index;
    VersionedValue<T>& vv = cache_[idx];
    // if the current value is stale, compute a new one and cache it
    if (vv.gen < current_gen_) {
      vv.val = gen_value(idx);
      vv.gen = current_gen_;
    }
    // return the (possibly newly) cached value
    return vv.val;
  }

  private:
  std::vector<VersionedValue<T>> cache_;
  size_t max_index_{0};
  uint64_t current_gen_{0};
};

} // namespace distribution_utils

#endif // __DISTRIBUTION_UTILS__
