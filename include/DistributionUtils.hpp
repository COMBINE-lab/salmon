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

/**
 *  Draw samples from the provided fragment length distribution.
 *  \param fld A pointer to the FragmentLengthDistribution from which
 *             samples will be drawn.
 *  \param numSamples  The number of samples to draw.
 */
std::vector<int32_t> samplesFromLogPMF(FragmentLengthDistribution* fld,
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

} // namespace distribution_utils

#endif // __DISTRIBUTION_UTILS__
