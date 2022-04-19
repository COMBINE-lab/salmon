#include <atomic>
#include <unordered_map>
#include <vector>
#include <exception>
#include <cmath>

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_reduce.h"
#include "tbb/partitioner.h"
// <-- deprecated in TBB --> #include "tbb/task_scheduler_init.h"
#include "tbb/global_control.h"

//#include "fastapprox.h"
#include <boost/math/special_functions/digamma.hpp>

// C++ string formatting library
#include "spdlog/fmt/fmt.h"

#include "Eigen/Dense"
#include "cuckoohash_map.hh"

#include "AlignmentLibrary.hpp"
#include "BootstrapWriter.hpp"
#include "CollapsedEMOptimizer.hpp"
#include "MultinomialSampler.hpp"
#include "ReadExperiment.hpp"
#include "ReadPair.hpp"
#include "SalmonMath.hpp"
#include "Transcript.hpp"
#include "TranscriptGroup.hpp"
#include "UnpairedRead.hpp"
#include "EMUtils.hpp"

using BlockedIndexRange = tbb::blocked_range<size_t>;

// intelligently chosen value originally adopted from
// https://github.com/pachterlab/kallisto/blob/master/src/EMAlgorithm.h#L18
// later modified since denorm_min seems to be too permissive.
constexpr double minEQClassWeight = std::numeric_limits<double>::min();
constexpr double minWeight = std::numeric_limits<double>::min();
// A bit more conservative of a minimum as an argument to the digamma function.
constexpr double digammaMin = 1e-10;

double normalize(std::vector<std::atomic<double>>& vec) {
  double sum{0.0};
  for (auto& v : vec) {
    sum += v;
  }

  // too small!
  if (sum < ::minWeight) {
    return sum;
  }

  double invSum = 1.0 / sum;
  for (auto& v : vec) {
    v.store(v.load() * invSum);
  }

  return sum;
}

template <typename VecT>
double truncateCountVector(VecT& alphas, std::vector<double>& cutoff) {
  // Truncate tiny expression values
  double alphaSum = 0.0;

  for (size_t i = 0; i < alphas.size(); ++i) {
    if (alphas[i] <= cutoff[i]) {
      alphas[i] = 0.0;
    }
    alphaSum += alphas[i];
  }
  return alphaSum;
}

/**
 *  Populate the prior parameters for the VBEM
 *  Note: effLens *must* be valid before calling this function.
 */
std::vector<double> populatePriorAlphas_(
    std::vector<Transcript>& transcripts, // transcripts
    Eigen::VectorXd& effLens,             // current effective length estimate
    double priorValue,      // the per-nucleotide prior value to use
    bool perTranscriptPrior // true if prior is per-txp, else per-nucleotide
) {
  // start out with the per-txp prior
  std::vector<double> priorAlphas(transcripts.size(), priorValue);

  // If the prior is per-nucleotide (default, then we need a potentially
  // different value for each transcript based on its length).
  if (!perTranscriptPrior) {
    for (size_t i = 0; i < transcripts.size(); ++i) {
      priorAlphas[i] = priorValue * effLens(i);
    }
  }
  return priorAlphas;
}

/**
 * Single-threaded VBEM-update routine for use in bootstrapping
 */
template <typename VecT>
void VBEMUpdate_(std::vector<std::vector<uint32_t>>& txpGroupLabels,
                 std::vector<std::vector<double>>& txpGroupCombinedWeights,
                 const std::vector<uint64_t>& txpGroupCounts,
                 std::vector<double>& priorAlphas, 
                 const VecT& alphaIn, VecT& alphaOut, VecT& expTheta) {

  assert(alphaIn.size() == alphaOut.size());
  size_t M = alphaIn.size();
  size_t numEQClasses = txpGroupLabels.size();
  double alphaSum = {0.0};
  for (size_t i = 0; i < M; ++i) {
    alphaSum += alphaIn[i] + priorAlphas[i];
  }

  double logNorm = boost::math::digamma(alphaSum);

  // double prior = priorAlpha;

  for (size_t i = 0; i < M; ++i) {
    auto ap = alphaIn[i] + priorAlphas[i];
    if (ap > ::digammaMin) {
        expTheta[i] =
          std::exp(boost::math::digamma(ap) - logNorm);
    } else {
      expTheta[i] = 0.0;
    }
    alphaOut[i] = 0.0; // priorAlphas[i];
  }

  for (size_t eqID = 0; eqID < numEQClasses; ++eqID) {
    uint64_t count = txpGroupCounts[eqID];
    const std::vector<uint32_t>& txps = txpGroupLabels[eqID];
    const auto& auxs = txpGroupCombinedWeights[eqID];

    size_t groupSize = txpGroupCombinedWeights[eqID].size(); // txps.size();
    // If this is a single-transcript group,
    // then it gets the full count.  Otherwise,
    // update according to our VBEM rule.
    if (BOOST_LIKELY(groupSize > 1)) {
      double denom = 0.0;
      for (size_t i = 0; i < groupSize; ++i) {
        auto tid = txps[i];
        auto aux = auxs[i];
        if (expTheta[tid] > 0.0) {
          double v = expTheta[tid] * aux;
          denom += v;
        }
      }
      if (denom <= ::minEQClassWeight) {
        // tgroup.setValid(false);
      } else {
        double invDenom = count / denom;
        for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          auto aux = auxs[i];
          if (expTheta[tid] > 0.0) {
            double v = expTheta[tid] * aux;
            salmon::utils::incLoop(alphaOut[tid], v * invDenom);
          }
        }
      }

    } else {
      salmon::utils::incLoop(alphaOut[txps.front()], count);
    }
  }
}

template <typename VecT>
void VBEMUpdate_augmented(std::vector<std::vector<uint32_t>>& txpGroupLabels,
                 std::vector<std::vector<double>>& txpGroupCombinedWeights,
                 const std::vector<uint64_t>& txpGroupCounts,
                 std::vector<double>& priorAlphas, 
                 const VecT& alphaIn, VecT& alphaOut, VecT& expTheta,
                 std::vector<uint32_t>& sampled_txps_counts,
                 std::vector<bool>& eqClass_augmentable) {


  assert(alphaIn.size() == alphaOut.size());
  size_t M = alphaIn.size();
  size_t numEQClasses = txpGroupLabels.size();
  double alphaSum = {0.0};
  for (size_t i = 0; i < M; ++i) {
    alphaSum += alphaIn[i] + priorAlphas[i];
  }

  double logNorm = boost::math::digamma(alphaSum);

  // double prior = priorAlpha;

  for (size_t i = 0; i < M; ++i) {
    auto ap = alphaIn[i] + priorAlphas[i];
    if (ap > ::digammaMin) {
        expTheta[i] =
          std::exp(boost::math::digamma(ap) - logNorm);
    } else {
      expTheta[i] = 0.0;
    }
    alphaOut[i] = 0.0; // priorAlphas[i];
  }

  for (size_t eqID = 0; eqID < numEQClasses; ++eqID) {
    uint64_t count = txpGroupCounts[eqID];
    const std::vector<uint32_t>& txps = txpGroupLabels[eqID];
    const auto& auxs = txpGroupCombinedWeights[eqID];

    size_t groupSize = txpGroupCombinedWeights[eqID].size(); // txps.size();
    // If this is a single-transcript group,
    // then it gets the full count.  Otherwise,
    // update according to our VBEM rule.
    if (BOOST_LIKELY(groupSize > 1)) {
      double denom = 0.0;
      for (size_t i = 0; i < groupSize; ++i) {
        auto tid = txps[i];
        auto aux = auxs[i];
        auto expTheta_i = expTheta[tid];
        if (!eqClass_augmentable[eqID]){
          auto ap = alphaIn[tid] + priorAlphas[tid] - static_cast<double>(sampled_txps_counts[tid]);
          if (ap > ::digammaMin) {
              expTheta_i =
                std::exp(boost::math::digamma(ap) - logNorm);
          } else {
            expTheta_i = 0.0;
          }          
        }
        if (expTheta_i > 0.0) {
          double v = expTheta_i * aux;
          denom += v;
        }
      }
      if (denom <= ::minEQClassWeight) {
        // tgroup.setValid(false);
      } else {
        double invDenom = count / denom;
        for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          auto aux = auxs[i];
          auto expTheta_i = expTheta[tid];
          if (!eqClass_augmentable[eqID]){
            auto ap = alphaIn[tid] + priorAlphas[tid] - static_cast<double>(sampled_txps_counts[tid]);
            if (ap > ::digammaMin) {
                expTheta_i =
                  std::exp(boost::math::digamma(ap) - logNorm);
            } else {
              expTheta_i = 0.0;
            }          
          }
          if (expTheta_i > 0.0) {
            double v = expTheta_i * aux;
            salmon::utils::incLoop(alphaOut[tid], v * invDenom);
          }
        }
      }

    } else {
      salmon::utils::incLoop(alphaOut[txps.front()], count);
    }
  }
}

template <typename VecT>
void EMUpdate_forAugmented(std::vector<std::vector<uint32_t>>& txpGroupLabels,
               std::vector<std::vector<double>>& txpGroupCombinedWeights,
               const std::vector<double>& txpGroupCounts,
               // std::unordered_set<uint32_t>& augmented_tids,
               std::vector<bool>& augmented_tids,
               std::vector<bool>& augmented_eqIDs,
               const VecT& alphaIn,
               VecT& alphaOut) {

  assert(alphaIn.size() == alphaOut.size());

  size_t numEqClasses = txpGroupLabels.size();
  for (size_t eqID = 0; eqID < numEqClasses; ++eqID) {
    if (!augmented_eqIDs[eqID]) {
      continue;
    }
 
    uint64_t count = txpGroupCounts[eqID];
    // for each transcript in this class
    const std::vector<uint32_t>& txps = txpGroupLabels[eqID];
    const auto& auxs = txpGroupCombinedWeights[eqID];

    double denom = 0.0;
    size_t groupSize = txpGroupCombinedWeights[eqID].size(); // txps.size();
    // If this is a single-transcript group,
    // then it gets the full count.  Otherwise,
    // update according to our VBEM rule.
    if (BOOST_LIKELY(groupSize > 1)) {
      for (size_t i = 0; i < groupSize; ++i) {
        auto tid = txps[i];
        if (augmented_tids[tid]) {
          auto aux = auxs[i];
          double v = alphaIn[tid] * aux;
          denom += v;
        }
      }

      if (denom <=  std::numeric_limits<double>::denorm_min()) {
        // tgroup.setValid(false);
      } else {
        double invDenom = count / denom;
        for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          if (augmented_tids[tid]) {
            auto aux = auxs[i];
            double v = alphaIn[tid] * aux;
            if (!std::isnan(v) && augmented_tids[tid]) {
              salmon::utils::incLoop(alphaOut[tid], v * invDenom);
            }
          }
        }
      }
    } else {
      salmon::utils::incLoop(alphaOut[txps.front()], count);
    }
  }
}


/*
 * Use the "standard" EM algorithm over equivalence
 * classes to estimate the latent variables (alphaOut)
 * given the current estimates (alphaIn).
 */
template <typename EQVecT>
void EMUpdate_(EQVecT& eqVec,
               std::vector<double>& priorAlphas,
               const CollapsedEMOptimizer::VecType& alphaIn,
               CollapsedEMOptimizer::VecType& alphaOut) {

  assert(alphaIn.size() == alphaOut.size());

  tbb::parallel_for(
      BlockedIndexRange(size_t(0), size_t(eqVec.size())),
      [&eqVec, &priorAlphas, &alphaIn, &alphaOut](const BlockedIndexRange& range) -> void {
        for (auto eqID : boost::irange(range.begin(), range.end())) {
          auto& kv = eqVec[eqID];

          uint64_t count = kv.second.count;
          // for each transcript in this class
          const TranscriptGroup& tgroup = kv.first;
          if (tgroup.valid) {
            const std::vector<uint32_t>& txps = tgroup.txps;
            const auto& auxs = kv.second.combinedWeights;

            size_t groupSize = kv.second.weights.size(); // txps.size();
            // If this is a single-transcript group,
            // then it gets the full count.  Otherwise,
            // update according to our VBEM rule.
            if (BOOST_LIKELY(groupSize > 1)) {
              double denom = 0.0;
              for (size_t i = 0; i < groupSize; ++i) {
                auto tid = txps[i];
                auto aux = auxs[i];
                double v = (alphaIn[tid]) * aux;
                denom += v;
              }

              if (denom <= ::minEQClassWeight) {
                // tgroup.setValid(false);
              } else {
                double invDenom = count / denom;
                for (size_t i = 0; i < groupSize; ++i) {
                  auto tid = txps[i];
                  auto aux = auxs[i];
                  double v = (alphaIn[tid]) * aux;
                  if (!std::isnan(v)) {
                    salmon::utils::incLoop(alphaOut[tid], v * invDenom);
                  }
                }
              }
            } else {
              salmon::utils::incLoop(alphaOut[txps.front()], count);
            }
          }
        }
      });
}

template <typename EQVecT>
void EMUpdate_augmented(EQVecT& eqVec,
               std::vector<double>& priorAlphas,
               const CollapsedEMOptimizer::VecType& alphaIn,
               CollapsedEMOptimizer::VecType& alphaOut,
               std::vector<uint32_t>& sampled_txps_counts,
               std::vector<bool>& eqClass_augmentable) {

  assert(alphaIn.size() == alphaOut.size());

  tbb::parallel_for(
      BlockedIndexRange(size_t(0), size_t(eqVec.size())),
      [&eqVec, &priorAlphas, &alphaIn, &alphaOut,
      &sampled_txps_counts, &eqClass_augmentable](const BlockedIndexRange& range) -> void {
        for (auto eqID : boost::irange(range.begin(), range.end())) {
          auto& kv = eqVec[eqID];

          uint64_t count = kv.second.count;
          // for each transcript in this class
          const TranscriptGroup& tgroup = kv.first;
          if (tgroup.valid) {
            const std::vector<uint32_t>& txps = tgroup.txps;
            const auto& auxs = kv.second.combinedWeights;

            size_t groupSize = kv.second.weights.size(); // txps.size();
            // If this is a single-transcript group,
            // then it gets the full count.  Otherwise,
            // update according to our VBEM rule.
            if (BOOST_LIKELY(groupSize > 1)) {
              double denom = 0.0;
              for (size_t i = 0; i < groupSize; ++i) {
                auto tid = txps[i];
                auto aux = auxs[i];
                // double augmentedCount = eqClass_augmentable[eqID] ? 0 : static_cast<double>(sampled_txps_counts[tid]);
                /*if (augmentedCount > 0) {
                  std::cerr<< i << ":" << tid << "\t" 
                    << eqClass_augmentable[i] << " " 
                    << augmentedCount << " " 
                    << alphaIn[tid] << "\n";
                }*/
                double v = alphaIn[tid] * aux;
                denom += v;
              }

              if (denom <= ::minEQClassWeight) {
                // tgroup.setValid(false);
              } else {
                double invDenom = count / denom;
                for (size_t i = 0; i < groupSize; ++i) {
                  auto tid = txps[i];
                  auto aux = auxs[i];
                  // double augmentedCount = eqClass_augmentable[i] ? 0 : static_cast<double>(sampled_txps_counts[tid]);
                  double v = alphaIn[tid] * aux;
                  if (!std::isnan(v)) {
                    salmon::utils::incLoop(alphaOut[tid], v * invDenom);
                  }
                }
              }
            } else {
              salmon::utils::incLoop(alphaOut[txps.front()], count);
            }
          }
        }
      });
}

/*
 * Use the Variational Bayesian EM algorithm over equivalence
 * classes to estimate the latent variables (alphaOut)
 * given the current estimates (alphaIn).
 */
template <typename EQVecT>
void VBEMUpdate_(EQVecT& eqVec,
                 std::vector<double>& priorAlphas, 
                 const CollapsedEMOptimizer::VecType& alphaIn,
                 CollapsedEMOptimizer::VecType& alphaOut,
                 CollapsedEMOptimizer::VecType& expTheta) {

  assert(alphaIn.size() == alphaOut.size());
  size_t M = alphaIn.size();
  double alphaSum = {0.0};
  for (size_t i = 0; i < M; ++i) {
    alphaSum += alphaIn[i] + priorAlphas[i];
  }

  double logNorm = boost::math::digamma(alphaSum);

  tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(priorAlphas.size())),
                    [logNorm, &priorAlphas, &alphaIn, &alphaOut,
                     &expTheta](const BlockedIndexRange& range) -> void {

                      // double prior = priorAlpha;

                      for (auto i : boost::irange(range.begin(), range.end())) {
                        auto ap = alphaIn[i].load() + priorAlphas[i];
                        if (ap > ::digammaMin) {
                          expTheta[i] =
                              std::exp(boost::math::digamma(ap) - logNorm);
                        } else {
                          expTheta[i] = 0.0;
                        }
                        // alphaOut[i] = prior * transcripts[i].RefLength;
                        alphaOut[i] = 0.0;
                      }
                    });

  tbb::parallel_for(
      BlockedIndexRange(size_t(0), size_t(eqVec.size())),
      [&eqVec, &alphaIn, &alphaOut,
       &expTheta](const BlockedIndexRange& range) -> void {
        for (auto eqID : boost::irange(range.begin(), range.end())) {
          auto& kv = eqVec[eqID];

          uint64_t count = kv.second.count;
          // for each transcript in this class
          const TranscriptGroup& tgroup = kv.first;
          if (tgroup.valid) {
            const std::vector<uint32_t>& txps = tgroup.txps;
            const auto& auxs = kv.second.combinedWeights;

            size_t groupSize = kv.second.weights.size(); // txps.size();
            // If this is a single-transcript group,
            // then it gets the full count.  Otherwise,
            // update according to our VBEM rule.
            if (BOOST_LIKELY(groupSize > 1)) {
              double denom = 0.0;
              for (size_t i = 0; i < groupSize; ++i) {
                auto tid = txps[i];
                auto aux = auxs[i];
                if (expTheta[tid] > 0.0) {
                  double v = expTheta[tid] * aux;
                  denom += v;
                }
              }
              if (denom <= ::minEQClassWeight) {
                // tgroup.setValid(false);
              } else {
                double invDenom = count / denom;
                for (size_t i = 0; i < groupSize; ++i) {
                  auto tid = txps[i];
                  auto aux = auxs[i];
                  if (expTheta[tid] > 0.0) {
                    double v = expTheta[tid] * aux;
                    salmon::utils::incLoop(alphaOut[tid], v * invDenom);
                  }
                }
              }

            } else {
              salmon::utils::incLoop(alphaOut[txps.front()], count);
            }
          }
        }
      });
}

template <typename EQVecT>
void VBEMUpdate_augmented(EQVecT& eqVec,
                 std::vector<double>& priorAlphas, 
                 const CollapsedEMOptimizer::VecType& alphaIn,
                 CollapsedEMOptimizer::VecType& alphaOut,
                 CollapsedEMOptimizer::VecType& expTheta,
                 std::vector<uint32_t>& sampled_txps_counts,
                 std::vector<bool>& eqClass_augmentable) {

  assert(alphaIn.size() == alphaOut.size());
  size_t M = alphaIn.size();
  double alphaSum = {0.0};
  for (size_t i = 0; i < M; ++i) {
    alphaSum += alphaIn[i] + priorAlphas[i];
  }

  double logNorm = boost::math::digamma(alphaSum);

  tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(priorAlphas.size())),
                    [logNorm, &priorAlphas, &alphaIn, &alphaOut,
                     &expTheta](const BlockedIndexRange& range) -> void {

                      // double prior = priorAlpha;

                      for (auto i : boost::irange(range.begin(), range.end())) {
                        auto ap = alphaIn[i].load() + priorAlphas[i];
                        if (ap > ::digammaMin) {
                          expTheta[i] =
                              std::exp(boost::math::digamma(ap) - logNorm);
                        } else {
                          expTheta[i] = 0.0;
                        }
                        // alphaOut[i] = prior * transcripts[i].RefLength;
                        alphaOut[i] = 0.0;
                      }
                    });

  tbb::parallel_for(
      BlockedIndexRange(size_t(0), size_t(eqVec.size())),
      [&eqVec, &alphaIn, &alphaOut, &priorAlphas, &logNorm, &eqClass_augmentable, &sampled_txps_counts,
       &expTheta](const BlockedIndexRange& range) -> void {
        for (auto eqID : boost::irange(range.begin(), range.end())) {
          auto& kv = eqVec[eqID];

          uint64_t count = kv.second.count;
          // for each transcript in this class
          const TranscriptGroup& tgroup = kv.first;
          if (tgroup.valid) {
            const std::vector<uint32_t>& txps = tgroup.txps;
            const auto& auxs = kv.second.combinedWeights;

            size_t groupSize = kv.second.weights.size(); // txps.size();
            // If this is a single-transcript group,
            // then it gets the full count.  Otherwise,
            // update according to our VBEM rule.
            if (BOOST_LIKELY(groupSize > 1)) {
              double denom = 0.0;
              for (size_t i = 0; i < groupSize; ++i) {
                auto tid = txps[i];
                auto aux = auxs[i];
                auto expTheta_i = expTheta[tid];
                if (!eqClass_augmentable[eqID]) {
                  auto ap = alphaIn[tid].load() + priorAlphas[tid] - static_cast<double>(sampled_txps_counts[tid]);
                  if (ap > ::digammaMin) {
                    expTheta_i =
                        std::exp(boost::math::digamma(ap) - logNorm);
                  } else {
                    expTheta_i = 0.0;
                  }
                }
                if ( expTheta_i > 0.0) {
                  double v = expTheta_i * aux;
                  denom += v;
                }
              }
              if (denom <= ::minEQClassWeight) {
                // tgroup.setValid(false);
              } else {
                double invDenom = count / denom;
                for (size_t i = 0; i < groupSize; ++i) {
                  auto tid = txps[i];
                  auto aux = auxs[i];
                  auto expTheta_i = expTheta[tid];
                  if (!eqClass_augmentable[eqID]) {
                    auto ap = alphaIn[tid].load() + priorAlphas[tid] - static_cast<double>(sampled_txps_counts[tid]);
                    if (ap > ::digammaMin) {
                      expTheta_i =
                          std::exp(boost::math::digamma(ap) - logNorm);
                    } else {
                      expTheta_i = 0.0;
                    }
                  }
                  if (expTheta_i > 0.0) {
                    double v = expTheta_i * aux;
                    salmon::utils::incLoop(alphaOut[tid], v * invDenom);
                  }
                }
              }

            } else {
              salmon::utils::incLoop(alphaOut[txps.front()], count);
            }
          }
        }
      });
}


template <typename VecT, typename EQVecT>
size_t markDegenerateClasses(
    EQVecT& eqVec,
    VecT& alphaIn, std::vector<bool>& available,
    std::shared_ptr<spdlog::logger> jointLog, bool verbose = false) {

  size_t numDropped{0};
  for (auto& kv : eqVec) {
    uint64_t count = kv.second.count;
    // for each transcript in this class
    const TranscriptGroup& tgroup = kv.first;
    const std::vector<uint32_t>& txps = tgroup.txps;
    const auto& auxs = kv.second.combinedWeights;

    double denom = 0.0;
    size_t groupSize = kv.second.weights.size();
    for (size_t i = 0; i < groupSize; ++i) {
      auto tid = txps[i];
      auto aux = auxs[i];
      double v = alphaIn[tid] * aux;
      if (!std::isnan(v)) {
        denom += v;
      } else {
        std::cerr << "val is NAN; alpha( " << tid << " ) = " << alphaIn[tid]
                  << ", aux = " << aux << "\n";
      }
    }
    if (denom <= ::minEQClassWeight) {
      if (verbose) {
        fmt::MemoryWriter errstream;

        errstream << "\nDropping weighted eq class\n";
        errstream << "============================\n";

        errstream << "denom = 0, count = " << count << "\n";
        errstream << "class = { ";
        for (size_t tn = 0; tn < groupSize; ++tn) {
          errstream << txps[tn] << " ";
        }
        errstream << "}\n";
        errstream << "alphas = { ";
        for (size_t tn = 0; tn < groupSize; ++tn) {
          errstream << alphaIn[txps[tn]] << " ";
        }
        errstream << "}\n";
        errstream << "weights = { ";
        for (size_t tn = 0; tn < groupSize; ++tn) {
          errstream << auxs[tn] << " ";
        }
        errstream << "}\n";
        errstream << "============================\n\n";

        jointLog->info(errstream.str());
      }
      ++numDropped;
      kv.first.setValid(false);
    } else {
      for (size_t i = 0; i < groupSize; ++i) {
        auto tid = txps[i];
        available[tid] = true;
      }
    }
  }
  return numDropped;
}

CollapsedEMOptimizer::CollapsedEMOptimizer() {}

bool doBootstrap(
    std::vector<std::vector<uint32_t>>& txpGroups,
    std::vector<std::vector<double>>& txpGroupCombinedWeights,
    std::vector<Transcript>& transcripts, Eigen::VectorXd& effLens,
    const std::vector<double>& sampleWeights, std::vector<uint64_t>& origCounts,
    uint64_t totalNumFrags,
    uint64_t numMappedFrags, double uniformTxpWeight,
    std::atomic<uint32_t>& bsNum, SalmonOpts& sopt,
    std::vector<double>& priorAlphas,
    std::function<bool(const std::vector<double>&)>& writeBootstrap,
    double relDiffTolerance, uint32_t maxIter,
    uint32_t single_count_class_offset,
    // transcripts that have an eq id collision with some other
    // txp; useful in the augmented setting. This vector is
    // sorted in ascending order.
    std::vector<uint32_t>& eq_identical_txps,
    std::atomic<uint32_t>& sampled_txps_total) {

  // An EM termination criterion, adopted from Bray et al. 2016
  uint32_t minIter = 50;

  // Determine up front if we're going to use scaled counts.
  bool useScaledCounts = !(sopt.useQuasi or sopt.allowOrphans);
  bool useVBEM{sopt.useVBOpt};
  size_t numClasses = txpGroups.size();
  CollapsedEMOptimizer::SerialVecType alphas(transcripts.size(), 0.0);
  CollapsedEMOptimizer::SerialVecType alphasPrime(transcripts.size(), 0.0);
  CollapsedEMOptimizer::SerialVecType expTheta(transcripts.size(), 0.0);
  std::vector<uint64_t> sampCounts(numClasses, 0);
  std::vector<double> newCounts(numClasses, 0);

  uint32_t numBootstraps = sopt.numBootstraps;
  bool perTranscriptPrior{sopt.perTranscriptPrior};

  auto& jointLog = sopt.jointLog;

  #if defined(__linux) && defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128
    std::random_device rd("/dev/urandom");
  #else
    std::random_device rd;
  #endif  // defined(__GLIBCXX__) && __GLIBCXX__ >= 2020012

  std::mt19937 gen(rd());
  // MultinomialSampler msamp(rd);
  std::discrete_distribution<uint64_t> csamp_orig(sampleWeights.begin(),
                                             sampleWeights.begin() + single_count_class_offset);
  std::discrete_distribution<uint64_t> csamp_aug(sampleWeights.begin() + single_count_class_offset,
                                             sampleWeights.end());

  std::vector<bool> augmented_tids(transcripts.size(), false);
  std::vector<bool> augmented_eqIDs(txpGroups.size(), false);

  for (size_t eqID = 0; eqID < txpGroups.size(); ++eqID) {
    if (eqID >= single_count_class_offset) {
      augmented_eqIDs[eqID] = true;
      continue;
    }
    const auto& txps = txpGroups[eqID];
    size_t groupSize = txps.size();
    for (size_t i = 0; i < groupSize; ++i) {
      auto tid = txps[i];
      if (std::binary_search(eq_identical_txps.begin(), eq_identical_txps.end(), tid)) {
        augmented_eqIDs[eqID] = true;
        augmented_tids[tid] = true;
      }
    }
  }

  while (bsNum++ < numBootstraps) {
    csamp_orig.reset();
    csamp_aug.reset();

    for (size_t sc = 0; sc < sampCounts.size(); ++sc) {
      sampCounts[sc] = 0;
    }
    std::vector<uint32_t> sampled_txps(transcripts.size(), 0);
    std::vector<uint32_t> sampled_txps_counts(transcripts.size(), 0);
    sampled_txps_total = 0;
    for (size_t fn = 0; fn < totalNumFrags; ++fn) {
      auto class_number = csamp_orig(gen);
      /*if (class_number >= single_count_class_offset) {
        std::cerr<<"Should not happen1!\n";
        exit(0);
        // auto tid = txpGroups[class_number].front();
        // sampled_txps[tid] = 1;
        // sampled_txps_counts[tid] += 1;
        // sampled_txps_total += 1;
      }*/
      ++sampCounts[class_number];
    }

    /* size_t aug_total = static_cast<size_t>(std::ceil(totalNumFrags * sopt.augmented_bootstrap_weight));
    for (size_t fn = 0; fn < aug_total; ++fn) {
      auto class_number = csamp_aug(gen);
      if (class_number >= sampleWeights.size() - single_count_class_offset) {
        std::cerr<<"Should not happen2!\n";
        exit(0);
      }
      auto global_eq_class_number = single_count_class_offset + class_number;
      auto tid = txpGroups[global_eq_class_number].front();
      sampled_txps_counts[tid] += 1;
      ++sampCounts[global_eq_class_number];
    } */

    // Do a new bootstrap
    double totalLen{0.0};
    for (size_t i = 0; i < transcripts.size(); ++i) {
      alphas[i] =
          transcripts[i].getActive() ? uniformTxpWeight * totalNumFrags : 0.0;
      totalLen += effLens(i);
    }

    bool converged{false};
    double maxRelDiff = -std::numeric_limits<double>::max();
    size_t itNum = 0;

    // If we use VBEM, we'll need the prior parameters
    // double priorAlpha = 1.00;

    // EM termination criteria, adopted from Bray et al. 2016
    double minAlpha = 1e-8;
    double alphaCheckCutoff = 1e-2;
    double cutoff = minAlpha;

    /*std::vector<bool> eqClass_augmentable(txpGroups.size(), false);
    for (size_t eqID = 0; eqID < txpGroups.size(); ++eqID) {
      const auto& txps = txpGroups[eqID];
      if (eqID >= single_count_class_offset) {
        // eqClass_augmentable[eqID] = true;
        continue;
      }
      size_t groupSize = txps.size();
      for (size_t i = 0; i < groupSize; ++i) {
        auto tid = txps[i];
        if (std::binary_search(eq_identical_txps.begin(), eq_identical_txps.end(), tid)) {
          eqClass_augmentable[eqID] = true;
          break;
        }
      }
      // original approach
      /*
      // for each transcript in this class
      // const TranscriptGroup& tgroup = kv.first;
      // const std::vector<uint32_t>& txps = tgroup.txps;
      size_t groupSize = txps.size();
      for (size_t i = 0; i < groupSize; ++i) {
      auto tid = txps[i];
      if (transcripts[tid].uniqueCount() > 0) {
      eqClass_augmentable[eqID] = false;
      break;
      }
      }
      */
    //}

    size_t num_augmentable = 0;
    for (auto a : augmented_eqIDs) {      
      num_augmentable += a;
    }
    size_t num_augmentable_txps = 0;
    for (auto a : augmented_tids) {      
      num_augmentable_txps += a;
    }

    auto augmented_class_count = txpGroups.size() - single_count_class_offset;
    num_augmentable -= augmented_class_count;
    std::cerr << "added class count (number of augmentable txps) = "
              << augmented_class_count << " " << num_augmentable_txps << "\n";
    std::cerr << "number of augmentable eq classes = " << num_augmentable << "\n";
    std::cerr << "eq identical txps = " << eq_identical_txps.size() << "\n";
    std::cerr << "single_count_class_offset = " << single_count_class_offset << "\n";
    std::cerr << "num eq classes = " << txpGroups.size() << "\n";
    std::cerr << (num_augmentable) << " out of " 
              << single_count_class_offset << " equivalence classes were augmentable\n";

    while (itNum < minIter or (itNum < maxIter and !converged)) {
      if (useVBEM) {
        if (sopt.eqClassBasedAugmentation)
          VBEMUpdate_augmented(txpGroups, txpGroupCombinedWeights, sampCounts, 
                    priorAlphas, alphas, alphasPrime, expTheta, sampled_txps_counts, augmented_eqIDs);
        else
          VBEMUpdate_(txpGroups, txpGroupCombinedWeights, sampCounts, 
                    priorAlphas, alphas, alphasPrime, expTheta);
      } else {
        if (false and sopt.eqClassBasedAugmentation and (itNum > 0)) {
          EMUpdate_augmented(txpGroups, txpGroupCombinedWeights, sampCounts,
                alphas, alphasPrime, sampled_txps_counts, augmented_eqIDs);
        } else {
          EMUpdate_(txpGroups, txpGroupCombinedWeights, sampCounts, 
                alphas, alphasPrime);
        }
      }

      converged = true;
      maxRelDiff = -std::numeric_limits<double>::max();
      for (size_t i = 0; i < transcripts.size(); ++i) {
        if (alphasPrime[i] > alphaCheckCutoff) {
          double relDiff =
              std::abs(alphas[i] - alphasPrime[i]) / alphasPrime[i];
          maxRelDiff = (relDiff > maxRelDiff) ? relDiff : maxRelDiff;
          if (relDiff > relDiffTolerance) {
            converged = false;
          }
        }
        alphas[i] = alphasPrime[i];
        alphasPrime[i] = 0.0;
      }

      ++itNum;
    }

    // Perform augmentation for augmentable transcripts
    size_t aug_total = static_cast<size_t>(std::ceil(totalNumFrags * sopt.augmented_bootstrap_weight));
    for (size_t fn = 0; fn < aug_total; ++fn) {
      auto class_number = csamp_aug(gen);
      if (class_number >= sampleWeights.size() - single_count_class_offset) {
        std::cerr<<"Should not happen2!\n";
        std::exit(0);
      }
      auto global_eq_class_number = single_count_class_offset + class_number;
      auto tid = txpGroups[global_eq_class_number].front();
      sampled_txps_counts[tid] += 1;
      ++newCounts[global_eq_class_number];
    }
    // Create new set of eq classes with augmentable txps
    size_t numEqClasses = txpGroups.size();
    for (size_t eqID = 0; eqID < numEqClasses; ++eqID) {
      if (!augmented_eqIDs[eqID]) {
        continue;
      }
      size_t groupSize = txpGroupCombinedWeights[eqID].size();
      const std::vector<uint32_t>& txps = txpGroups[eqID];
      const auto& auxs = txpGroupCombinedWeights[eqID];
      uint64_t count = sampCounts[eqID];
      if (BOOST_LIKELY(groupSize > 1)) {
        double new_count = 0;
        double denom = 0.0;

        for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          auto aux = auxs[i];
          double v = alphas[tid] * aux;
          denom += v;
        }
        if (denom <=  std::numeric_limits<double>::denorm_min()) {
          
        } else {
          double invDenom = count / denom;
          for (size_t i = 0; i < groupSize; ++i) {
            auto tid = txps[i];
            auto aux = auxs[i];
            double v = alphas[tid] * aux;
            if (augmented_tids[tid]) {
              salmon::utils::incLoop(alphasPrime[tid], v * invDenom);
              salmon::utils::incLoop(new_count, v * invDenom);          
            }
          }
          newCounts[eqID] = new_count;
        }      
      } else {
        salmon::utils::incLoop(alphasPrime[txps.front()], count);
        newCounts[eqID] = static_cast<double>(count);
      }
    }

    for (size_t tid = 0; tid < alphas.size(); tid++) {
      if (augmented_tids[tid]) {
        alphas[tid] = alphasPrime[tid];
        alphasPrime[tid] = 0;
      }
    }

    // EM2 after augmentation with augmented transcripts only being able to change
    while (itNum < minIter or (itNum < maxIter and !converged)) {
      EMUpdate_forAugmented(txpGroups, txpGroupCombinedWeights, newCounts, 
            augmented_tids, augmented_eqIDs, alphas, alphasPrime);

      converged = true;
      maxRelDiff = -std::numeric_limits<double>::max();
      for (size_t i = 0; i < transcripts.size(); ++i) {
        if (augmented_tids[i]) {
          if (alphasPrime[i] > alphaCheckCutoff) {
            double relDiff =
                std::abs(alphas[i] - alphasPrime[i]) / alphasPrime[i];
            maxRelDiff = (relDiff > maxRelDiff) ? relDiff : maxRelDiff;
            if (relDiff > relDiffTolerance) {
              converged = false;
            }
          }
          alphas[i] = alphasPrime[i];
          alphasPrime[i] = 0.0;
        }
      }

      ++itNum;
    }

    /*std::vector<bool> augmented_tids(transcripts.size(), false);
    for (auto eqID : eqClass_augmentable) {
      if (eqClass_augmentable[eqID]) {
        for (auto tid: txpGroups[eqID])
          if (sampled_txps_counts[tid] > 0)
            augmented_tids[tid] = true;
      }
    }*/

    for (size_t i = 0; i < sampled_txps_counts.size(); ++i) { //single_count_class_offset; i < txpGroups.size(); ++i) {
      double adjusted_count = alphas[i] - sampled_txps_counts[i];
      // if (adjusted_count < 0.0 ) { std::cerr << "adjusted count " << i << " = " << adjusted_count; }
      alphas[i] = std::max(0.0, adjusted_count);
    }

    // Consider the projection of the abundances onto the *original* equivalence class
    // counts
    if (sopt.bootstrapReproject) {
      if (useVBEM) {
        VBEMUpdate_(txpGroups, txpGroupCombinedWeights, origCounts, 
                    priorAlphas, alphas, alphasPrime, expTheta);
      } else {
        EMUpdate_(txpGroups, txpGroupCombinedWeights, origCounts, 
                  alphas, alphasPrime);
      }
    }

    // Truncate tiny expression values
    double alphaSum = 0.0;
    if (useVBEM and !perTranscriptPrior) {
      std::vector<double> cutoffs(transcripts.size(), 0.0);
      for (size_t i = 0; i < transcripts.size(); ++i) {
        cutoffs[i] = minAlpha;
      }
      // alphaSum = truncateCountVector(alphas, cutoffs);
      alphaSum = truncateCountVector(alphas, cutoffs);
    } else {
      // Truncate tiny expression values
      alphaSum = truncateCountVector(alphas, cutoff);
    }

    if (alphaSum < ::minWeight) {
      jointLog->error("Total alpha weight was too small! "
                      "Make sure you ran salmon correctly.");
      return false;
    }

    if (useScaledCounts) {
      double mappedFragsDouble = static_cast<double>(numMappedFrags);
      double alphaSum = 0.0;
      for (auto a : alphas) {
        alphaSum += a;
      }
      if (alphaSum > ::minWeight) {
        double scaleFrac = 1.0 / alphaSum;
        // scaleFrac converts alpha to nucleotide fraction,
        // and multiplying by numMappedFrags scales by the total
        // number of mapped fragments to provide an estimated count.
        for (auto& a : alphas) {
          a = mappedFragsDouble * (a * scaleFrac);
        }
      } else { // This shouldn't happen!
        sopt.jointLog->error(
            "Bootstrap had insufficient number of fragments!"
            "Something is probably wrong; please check that you "
            "have run salmon correctly and report this to GitHub.");
      }
    }

    writeBootstrap(alphas);
  }
  return true;
}

std::vector<uint32_t> get_eq_identical_txps(std::vector<Transcript>& transcripts,
                                            std::vector<std::vector<uint32_t>>& txpGroups,
                                            std::vector<uint32_t>& aug_collision_txps) {

  // build a map from transcripts to the equivalence class IDs in which they appear.
  std::vector<std::vector<uint32_t>> txp_to_eqid;
  txp_to_eqid.resize(transcripts.size(), std::vector<uint32_t>());

  std::cerr << "building txp to eqid map\n";
  // here the idea is, for each transcript, we push back onto a vector its
  // equiavlence class "signature"; the set of equivalence classes to which
  // it belongs.
  for (size_t i = 0; i < txpGroups.size(); ++i) {
    auto& tg = txpGroups[i];
    for (auto t : tg) {
      txp_to_eqid[t].push_back(static_cast<uint32_t>(i));
    }
  }

  struct EqIdVecHasher {
    std::size_t operator()(const std::vector<uint32_t>& k) const {
      size_t seed{0};
      return XXH64(static_cast<void*>(const_cast<uint32_t*>(k.data())), k.size() * sizeof(uint32_t),
                   seed);
    }
  };

  std::cerr << "building eq_id_vec to txp map\n";
  // map from the equivalence-class-id set to the transcripts having this eq-class signature
  // here the idea is to use hashing to efficiently determing the set of transcripts with
  // exactly the same equivalence class signature.
  std::unordered_map<std::vector<uint32_t>, std::vector<uint32_t>, EqIdVecHasher> eq_id_to_txps;
  for (size_t i = 0; i < txp_to_eqid.size(); ++i) {
    auto& eqvec = txp_to_eqid[i];
    if (!eqvec.empty()) {
      eq_id_to_txps[eqvec].push_back(static_cast<uint32_t>(i));
    }
  }

  std::vector<uint32_t> eq_sig_collision_txps;

  std::cerr << "finding target transcripts\n";
  for (auto kv = eq_id_to_txps.begin(); kv != eq_id_to_txps.end(); ++kv) {
    auto& v = kv->second;
    if (v.size() > 1) {
      for (auto& t : v ) { eq_sig_collision_txps.push_back(t); }
      auto& eq_ids = txp_to_eqid[v[0]];
      for (auto& eq: eq_ids) {
        auto& tg = txpGroups[eq];
        for (auto& tid: tg) { aug_collision_txps.push_back(tid); }
      }
    }
  }

  std::sort(eq_sig_collision_txps.begin(), eq_sig_collision_txps.end());
  eq_sig_collision_txps.erase(std::unique( eq_sig_collision_txps.begin(),
                                           eq_sig_collision_txps.end() ), eq_sig_collision_txps.end() );

  std::sort(aug_collision_txps.begin(), aug_collision_txps.end());
  aug_collision_txps.erase(std::unique( aug_collision_txps.begin(),
                                        aug_collision_txps.end() ), aug_collision_txps.end() );

  std::cerr << "There were " << eq_sig_collision_txps.size() << " signature colliding transcripts\n";
  std::cerr << "There were " << aug_collision_txps.size() << " transcripts in the corresponding eq classes\n";
  return eq_sig_collision_txps;
}

template <typename ExpT>
bool CollapsedEMOptimizer::gatherBootstraps(
    ExpT& readExp, SalmonOpts& sopt,
    std::function<bool(const std::vector<double>&)>& writeBootstrap,
    double relDiffTolerance, uint32_t maxIter) {

  std::vector<Transcript>& transcripts = readExp.transcripts();
  std::vector<bool> available(transcripts.size(), false);
  using VecT = CollapsedEMOptimizer::SerialVecType;
  // With atomics
  VecT alphas(transcripts.size(), 0.0);
  VecT alphasPrime(transcripts.size(), 0.0);
  VecT expTheta(transcripts.size());
  Eigen::VectorXd effLens(transcripts.size());
  double minAlpha = 1e-8;

  bool scaleCounts = (!sopt.useQuasi and !sopt.allowOrphans);

  uint64_t numMappedFrags =
      scaleCounts ? readExp.upperBoundHits() : readExp.numMappedFragments();

  uint32_t numBootstraps = sopt.numBootstraps;

  auto& eqBuilder = readExp.equivalenceClassBuilder();
  auto& eqVec = eqBuilder.eqVec();

  std::unordered_set<uint32_t> activeTranscriptIDs;
  const size_t numClasses = eqVec.size();
  for (size_t cid = 0; cid < numClasses; ++cid) { 
    auto nt = eqBuilder.getNumTranscriptsForClass(cid);
    auto& txps = eqVec[cid].first.txps;
    for (size_t tctr = 0; tctr < nt; ++tctr) {
      auto t = txps[tctr];
      transcripts[t].setActive();
      activeTranscriptIDs.insert(t);
    }
  }

  bool perTranscriptPrior{sopt.perTranscriptPrior};
  double priorValue{sopt.vbPrior};

  auto jointLog = sopt.jointLog;

  jointLog->info("Will draw {:n} bootstrap samples", numBootstraps);
  jointLog->info("Optimizing over {:n} equivalence classes", eqVec.size());

  double totalNumFrags{static_cast<double>(numMappedFrags)};
  double totalLen{0.0};

  if (activeTranscriptIDs.size() == 0) {
    jointLog->error("It seems that no transcripts are expressed; something is "
                    "likely wrong!");
    std::exit(1);
  }

  double scale = 1.0 / activeTranscriptIDs.size();
  for (size_t i = 0; i < transcripts.size(); ++i) {
    // double m = transcripts[i].mass(false);
    alphas[i] = transcripts[i].getActive() ? scale * totalNumFrags : 0.0;
    effLens(i) = (sopt.noEffectiveLengthCorrection)
                     ? transcripts[i].RefLength
                     : transcripts[i].EffectiveLength;
    totalLen += effLens(i);
  }

  // If we use VBEM, we'll need the prior parameters
  std::vector<double> priorAlphas = populatePriorAlphas_(
      transcripts, effLens, priorValue, perTranscriptPrior);

  auto numRemoved =
      markDegenerateClasses(eqVec, alphas, available, sopt.jointLog);
  sopt.jointLog->info("Marked {} weighted equivalence classes as degenerate",
                      numRemoved);

  // Since we will use the same weights and transcript groups for each
  // of the bootstrap samples (only the count vector will change), it
  // makes sense to keep only one copy of these.
  using TGroupLabelT = std::vector<uint32_t>;
  using TGroupWeightVec = std::vector<double>;
  std::vector<TGroupLabelT> txpGroups;
  std::vector<TGroupWeightVec> txpGroupCombinedWeights;
  std::vector<uint64_t> origCounts;
  uint64_t totalCount{0};

  // this is redundant with the already existing 
  // activeTranscriptIDs
  // std::vector<bool> txps_presence(transcripts.size(), false);

  for (size_t cid = 0; cid < numClasses; ++cid) { 
    const auto& kv = eqVec[cid];
    uint64_t count = kv.second.count;

    // for each transcript in this class
    const TranscriptGroup& tgroup = kv.first;
    if (tgroup.valid) {
      //const std::vector<uint32_t>& txps = tgroup.txps;
      const auto numTranscripts = eqBuilder.getNumTranscriptsForClass(cid);
      std::vector<uint32_t> txps(tgroup.txps.begin(), tgroup.txps.begin()+numTranscripts);
      const auto& auxs = kv.second.combinedWeights;
      
      if (txps.size() != auxs.size()) {
        sopt.jointLog->critical(
            "# of transcripts ({}) should match length of weight vec. ({})",
            txps.size(), auxs.size());
        sopt.jointLog->flush();
        spdlog::drop_all();
        std::exit(1);
      }
      /*
      for (auto tid : txps) {
        txps_presence[tid] = true;
      }*/
      
      txpGroups.push_back(txps);
      // Convert to non-atomic
      txpGroupCombinedWeights.emplace_back(auxs.begin(), auxs.end());
      origCounts.push_back(count);
      totalCount += count;
    }
  }

  auto single_count_class_offset = txpGroups.size();
  uint32_t new_class_count = 0;

  bool augment_bootstraps = (sopt.augmented_bootstrap_weight > 0.0);
  std::vector<double> samplingWeights(txpGroups.size(), 0.0);

  std::vector<uint32_t> eq_identical_txps;
  std::vector<uint32_t> aug_collision_txps;

  if (augment_bootstraps) {

    // get a vector of transcripts that had an equivalence class
    // signature collision with at least one other transcript
    eq_identical_txps = get_eq_identical_txps(transcripts, txpGroups, aug_collision_txps);

    samplingWeights.resize(txpGroups.size() + activeTranscriptIDs.size(), 0.0); // transcripts.size(), 0.0);
    uint64_t n = totalCount;

    double aug_count = (static_cast<double>(totalCount) * sopt.augmented_bootstrap_weight) / 
                      (static_cast<double>(activeTranscriptIDs.size()));
    // for (auto& tid : activeTranscriptIDs) {
    // for (auto& tid : aug_collision_txps) {
    for (auto& tid : eq_identical_txps) {
      new_class_count += 1;
      txpGroups.push_back({tid});
      std::vector<double> auxs;
      double weight = 1.0;
      auxs.push_back(1.0);
      txpGroupCombinedWeights.emplace_back(auxs.begin(), auxs.end());
      uint64_t count = 1; // aug_count;
      origCounts.push_back(count);
      // totalCount += count;
    }
  }

  double floatCount = totalCount;
  // the original sampling distribution
  for (size_t i = 0; i < single_count_class_offset; ++i) {
    samplingWeights[i] = origCounts[i] / floatCount;
  }
  // the augmented sampling distribution is unfiorm
  for (size_t i = single_count_class_offset; i < origCounts.size(); ++i) {
    samplingWeights[i] = 1;
  }

  size_t numWorkerThreads{1};
  if (sopt.numThreads > 1 and numBootstraps > 1) {
    numWorkerThreads = std::min(sopt.numThreads - 1, numBootstraps - 1);
  }

  std::atomic<uint32_t> sampled_txps_total{0};
  std::atomic<uint32_t> bsCounter{0};
  std::vector<std::thread> workerThreads;
  for (size_t tn = 0; tn < numWorkerThreads; ++tn) {
    workerThreads.emplace_back(
        doBootstrap, std::ref(txpGroups), std::ref(txpGroupCombinedWeights),
        std::ref(transcripts), std::ref(effLens), std::ref(samplingWeights), std::ref(origCounts),
        totalCount, numMappedFrags, scale, std::ref(bsCounter), std::ref(sopt),
        std::ref(priorAlphas), std::ref(writeBootstrap), relDiffTolerance,
        maxIter, single_count_class_offset, std::ref(eq_identical_txps), std::ref(sampled_txps_total));
  }

  for (auto& t : workerThreads) {
    t.join();
  }

  /*std::ofstream outfile;
  outfile.open("/mnt/scratch3/mohsen/boost-boots/noor_data/active.list");
  for (auto &tid: aug_collision_txps) {
    outfile << transcripts[tid].RefName << "\n";
  }
  outfile.close();*/

  /*
  if (augment_bootstraps) {
    double nsampled_txps_rate = static_cast<double>(sampled_txps_total)/static_cast<double>(numBootstraps*new_class_count);
    std::cerr<<"-------------------------------------\n"
            <<"numExpTxps:"<<new_class_count
            <<"\nnumBootstraps:"<<numBootstraps
            <<"\nnumExpTxps * numBootstraps:"<<numBootstraps*new_class_count
            <<"\nsampled_txps_total:" << sampled_txps_total
            <<"\nsampled_txps_rate:"<< nsampled_txps_rate
            <<"\n-------------------------------------\n";
  }
  */

  return true;
}

template <typename EQVecT>
void updateEqClassWeights(
    EQVecT& eqVec,
    Eigen::VectorXd& effLens) {
  tbb::parallel_for(
      BlockedIndexRange(size_t(0), size_t(eqVec.size())),
      [&eqVec, &effLens](const BlockedIndexRange& range) -> void {
        // For each index in the equivalence class vector
        for (auto eqID : boost::irange(range.begin(), range.end())) {
          // The vector entry
          auto& kv = eqVec[eqID];
          // The label of the equivalence class
          const TranscriptGroup& k = kv.first;
          // The size of the label
          size_t classSize = kv.second.weights.size(); // k.txps.size();
          // The weights of the label
          auto& v = kv.second;

          // Iterate over each weight and set it equal to
          // 1 / effLen of the corresponding transcript
          double wsum{0.0};
          for (size_t i = 0; i < classSize; ++i) {
            auto tid = k.txps[i];
            auto probStartPos = 1.0 / effLens(tid);
            v.combinedWeights[i] =
                kv.second.count * (v.weights[i] * probStartPos);
            wsum += v.combinedWeights[i];
          }
          double wnorm = 1.0 / wsum;
          for (size_t i = 0; i < classSize; ++i) {
            v.combinedWeights[i] *= wnorm;
          }
        }
      });
}

template <typename ExpT>
bool CollapsedEMOptimizer::optimize(ExpT& readExp, SalmonOpts& sopt,
                                    double relDiffTolerance, uint32_t maxIter) {

  // <-- deprecated in TBB --> tbb::task_scheduler_init tbbScheduler(sopt.numThreads);
  tbb::global_control c(tbb::global_control::max_allowed_parallelism, sopt.numThreads);
  std::vector<Transcript>& transcripts = readExp.transcripts();
  std::vector<bool> available(transcripts.size(), false);

  // An EM termination criterion, adopted from Bray et al. 2016
  uint32_t minIter = 50;
  bool seqBiasCorrect = sopt.biasCorrect;
  bool gcBiasCorrect = sopt.gcBiasCorrect;
  bool posBiasCorrect = sopt.posBiasCorrect;
  bool doBiasCorrect = seqBiasCorrect or gcBiasCorrect or posBiasCorrect;
  bool metaGenomeMode = sopt.meta;
  bool altInitMode = sopt.alternativeInitMode;

  using VecT = CollapsedEMOptimizer::VecType;
  // With atomics
  VecType alphas(transcripts.size());
  VecType alphasPrime(transcripts.size());
  VecType expTheta(transcripts.size());

  Eigen::VectorXd effLens(transcripts.size());

  auto& eqVec =
      readExp.equivalenceClassBuilder().eqVec();

  bool noRichEq = sopt.noRichEqClasses;

  bool useVBEM{sopt.useVBOpt};
  bool perTranscriptPrior{sopt.perTranscriptPrior};
  double priorValue{sopt.vbPrior};

  auto jointLog = sopt.jointLog;

  auto& fragStartDists = readExp.fragmentStartPositionDistributions();
  double totalNumFrags{static_cast<double>(readExp.numMappedFragments())};
  double totalLen{0.0};

  // If effective length correction isn't turned off, then use effective
  // lengths rather than reference lengths.
  bool useEffectiveLengths = !sopt.noEffectiveLengthCorrection;

  int64_t numActive{0};
  double totalWeight{0.0};

  for (size_t i = 0; i < transcripts.size(); ++i) {
    auto& txp = transcripts[i];
    alphas[i] = txp.projectedCounts;
    totalWeight += alphas[i];
    effLens(i) = useEffectiveLengths
                     ? std::exp(txp.getCachedLogEffectiveLength())
                     : txp.RefLength;
    if (sopt.noLengthCorrection) {
      effLens(i) = 100.0;
    }
    txp.EffectiveLength = effLens(i);

    double uniqueCount = static_cast<double>(txp.uniqueCount() + 0.5);
    auto wi = (sopt.initUniform) ? 100.0 : (uniqueCount * 1e-3 * effLens(i));
    alphasPrime[i] = wi;
    ++numActive;
    totalLen += effLens(i);
  }

  // If we use VBEM, we'll need the prior parameters
  std::vector<double> priorAlphas = populatePriorAlphas_(
      transcripts, effLens, priorValue, perTranscriptPrior);

  // Based on the number of observed reads, use
  // a linear combination of the online estimates
  // and the uniform distribution.
  double uniformPrior = totalWeight / static_cast<double>(numActive);
  double maxFrac = 0.999;
  double fracObserved = std::min(maxFrac, totalWeight / sopt.numRequiredFragments);
  // Above, we placed the uniformative (uniform) initalization into the
  // alphasPrime variables.  If that's what the user requested, then copy those
  // over to the alphas
  if (sopt.initUniform) {
    for (size_t i = 0; i < alphas.size(); ++i) {
      alphas[i].store(alphasPrime[i].load());
      alphasPrime[i] = 1.0;
    }
  } else { // otherwise, initialize with a linear combination of the true and
           // uniform alphas
    for (size_t i = 0; i < alphas.size(); ++i) {
      auto uniAbund = (metaGenomeMode or altInitMode) ? alphasPrime[i].load()
                                                      : uniformPrior;
      alphas[i] =
          (alphas[i] * fracObserved) + (uniAbund * (1.0 - fracObserved));
      alphasPrime[i] = 1.0;
    }
  }

  // If the user requested *not* to use "rich" equivalence classes,
  // then wipe out all of the weight information here and simply replace
  // the weights with the effective length terms (here, the *inverse* of
  // the effective length).  Otherwise, multiply the existing weight terms
  // by the effective length term.
  tbb::parallel_for(
      BlockedIndexRange(size_t(0), size_t(eqVec.size())),
      [&eqVec, &effLens, noRichEq, &sopt](const BlockedIndexRange& range) -> void {
        // For each index in the equivalence class vector
        for (auto eqID : boost::irange(range.begin(), range.end())) {
          // The vector entry
          auto& kv = eqVec[eqID];
          // The label of the equivalence class
          const TranscriptGroup& k = kv.first;
          // The size of the label
          size_t classSize = kv.second.weights.size(); // k.txps.size();
          // The weights of the label
          auto& v = kv.second;

          // Iterate over each weight and set it
          double wsum{0.0};

          for (size_t i = 0; i < classSize; ++i) {
            auto tid = k.txps[i];
            double el = effLens(tid);
            if (el <= 1.0) {
              el = 1.0;
            }
            if (noRichEq) {
              // Keep length factor separate for the time being
              v.weights[i] = 1.0;
            }
            // meaningful values.
            auto probStartPos = 1.0 / el;

            // combined weight
            double wt = sopt.eqClassMode ? v.weights[i] : v.count * v.weights[i] * probStartPos;
            v.combinedWeights.push_back(wt);
            wsum += wt;
          }

          double wnorm = 1.0 / wsum;
          for (size_t i = 0; i < classSize; ++i) {
            v.combinedWeights[i] = v.combinedWeights[i] * wnorm;
          }
        }
      });

  auto numRemoved =
      markDegenerateClasses(eqVec, alphas, available, sopt.jointLog);
  sopt.jointLog->info("Marked {} weighted equivalence classes as degenerate",
                      numRemoved);

  size_t itNum{0};

  // EM termination criteria, adopted from Bray et al. 2016
  double minAlpha = 1e-8;
  double alphaCheckCutoff = 1e-2;
  double cutoff = minAlpha;

  // Iterations in which we will allow re-computing the effective lengths
  // if bias-correction is enabled.
  // std::vector<uint32_t> recomputeIt{100, 500, 1000};
  minIter = 100;

  bool converged{false};
  double maxRelDiff = -std::numeric_limits<double>::max();
  bool needBias = doBiasCorrect;
  size_t targetIt{10};
  /* -- v0.8.x
  double alphaSum = 0.0;
  */

  while (itNum < minIter or (itNum < maxIter and !converged) or needBias) {
    if (needBias and (itNum > targetIt or converged)) {

      jointLog->info(
          "iteration {:n}, adjusting effective lengths to account for biases",
          itNum);
      effLens = salmon::utils::updateEffectiveLengths(sopt, readExp, effLens,
                                                      alphas, available, true);
      // if we're doing the VB optimization, update the priors
      if (useVBEM) {
        priorAlphas = populatePriorAlphas_(transcripts, effLens, priorValue,
                                           perTranscriptPrior);
      }

      // Check for strangeness with the lengths.
      for (int32_t i = 0; i < effLens.size(); ++i) {
        if (effLens(i) <= 0.0) {
          jointLog->warn("Transcript {} had length {}", i, effLens(i));
        }
      }
      updateEqClassWeights(eqVec, effLens);
      needBias = false;

      if ( sopt.eqClassMode ) {
        jointLog->error("Eqclass Mode should not be performing bias correction");
        jointLog->flush();
        exit(1);
      }
    }

    if (useVBEM) {
      VBEMUpdate_(eqVec, priorAlphas, alphas,
                  alphasPrime, expTheta);
    } else {
      /*
      if (itNum > 0 and (itNum % 250 == 0)) {
        for (size_t i = 0; i < transcripts.size(); ++i) {
      	  if (alphas[i] < 1.0) { alphas[i] = 0.0; }
      	}
      }
      */

      EMUpdate_(eqVec, priorAlphas, alphas, alphasPrime);
    }

    converged = true;
    maxRelDiff = -std::numeric_limits<double>::max();
    for (size_t i = 0; i < transcripts.size(); ++i) {
      if (alphasPrime[i] > alphaCheckCutoff) {
        double relDiff = std::abs(alphas[i] - alphasPrime[i]) / alphasPrime[i];
        maxRelDiff = (relDiff > maxRelDiff) ? relDiff : maxRelDiff;
        if (relDiff > relDiffTolerance) {
          converged = false;
        }
      }
      alphas[i].store(alphasPrime[i].load());
      alphasPrime[i].store(0.0);
    }

    /* -- v0.8.x
    if (converged and itNum > minIter and !needBias) {
      if (useVBEM and !perTranscriptPrior) {
        std::vector<double> cutoffs(transcripts.size(), 0.0);
        for (size_t i = 0; i < transcripts.size(); ++i) {
          cutoffs[i] = minAlpha;
        }
        alphaSum = truncateCountVector(alphas, cutoffs);
      } else {
        // Truncate tiny expression values
        alphaSum = truncateCountVector(alphas, cutoff);
      }
      if (useVBEM) {
        VBEMUpdate_(eqVec, priorAlphas, alphas,
    alphasPrime, expTheta); } else { EMUpdate_(eqVec, transcripts, alphas,
    alphasPrime);
      }
      for (size_t i = 0; i < transcripts.size(); ++i) {
        alphas[i] = alphasPrime[i];
        alphasPrime[i] = 0.0;
      }
    }
    */

    if (itNum % 100 == 0) {
      jointLog->info("iteration = {:n} | max rel diff. = {}", itNum, maxRelDiff);
    }

    ++itNum;
  }

  /* -- v0.8.x
  if (alphaSum < ::minWeight) {
    jointLog->error("Total alpha weight was too small! "
                    "Make sure you ran salmon correctly.");
    return false;
  }
  */

  // Reset the original bias correction options
  sopt.gcBiasCorrect = gcBiasCorrect;
  sopt.biasCorrect = seqBiasCorrect;

  jointLog->info("iteration = {:n} | max rel diff. = {}", itNum, maxRelDiff);

  double alphaSum = 0.0;
  if (useVBEM and !perTranscriptPrior) {
    std::vector<double> cutoffs(transcripts.size(), 0.0);
    for (size_t i = 0; i < transcripts.size(); ++i) {
      cutoffs[i] = minAlpha;
    }
    alphaSum = truncateCountVector(alphas, cutoffs);
  } else {
    // Truncate tiny expression values
    alphaSum = truncateCountVector(alphas, cutoff);
  }

  if (alphaSum < ::minWeight) {
    jointLog->error("Total alpha weight was too small! "
                    "Make sure you ran salmon correctly.");
    return false;
  }

  // Set the mass of each transcript using the
  // computed alphas.
  for (size_t i = 0; i < transcripts.size(); ++i) {
    // Set the mass to the normalized (after truncation)
    // relative abundance
    // If we changed the effective lengths, copy them over here
    if (doBiasCorrect) {
      transcripts[i].EffectiveLength = effLens(i);
    }
    transcripts[i].setSharedCount(alphas[i]);
    transcripts[i].setMass(alphas[i] / alphaSum);
  }
  return true;
}

using BulkReadExperimentT = ReadExperiment<EquivalenceClassBuilder<TGValue>>;
template <typename FragT,typename AlignModelT>
using BulkAlnLibT = AlignmentLibrary<FragT, EquivalenceClassBuilder<TGValue>,AlignModelT>;
using SCReadExperimentT = ReadExperiment<EquivalenceClassBuilder<SCTGValue>>;


template bool CollapsedEMOptimizer::optimize<BulkReadExperimentT>(
    BulkReadExperimentT& readExp, SalmonOpts& sopt, double relDiffTolerance,
    uint32_t maxIter);

template bool CollapsedEMOptimizer::optimize<BulkAlnLibT<UnpairedRead,AlignmentModel>>(
    BulkAlnLibT<UnpairedRead,AlignmentModel>& readExp, SalmonOpts& sopt,
    double relDiffTolerance, uint32_t maxIter);

template bool CollapsedEMOptimizer::optimize<BulkAlnLibT<UnpairedRead,ONTAlignmentModel>>(
    BulkAlnLibT<UnpairedRead,ONTAlignmentModel>& readExp, SalmonOpts& sopt,
    double relDiffTolerance, uint32_t maxIter);

template bool CollapsedEMOptimizer::optimize<BulkAlnLibT<ReadPair,AlignmentModel>>(
    BulkAlnLibT<ReadPair,AlignmentModel>& readExp, SalmonOpts& sopt,
    double relDiffTolerance, uint32_t maxIter);


template bool CollapsedEMOptimizer::gatherBootstraps<BulkReadExperimentT>(
    BulkReadExperimentT& readExp, SalmonOpts& sopt,
    std::function<bool(const std::vector<double>&)>& writeBootstrap,
    double relDiffTolerance, uint32_t maxIter);

template bool
CollapsedEMOptimizer::gatherBootstraps<BulkAlnLibT<UnpairedRead,AlignmentModel>>(
    BulkAlnLibT<UnpairedRead,AlignmentModel>& readExp, SalmonOpts& sopt,
    std::function<bool(const std::vector<double>&)>& writeBootstrap,
    double relDiffTolerance, uint32_t maxIter);

template bool
CollapsedEMOptimizer::gatherBootstraps<BulkAlnLibT<UnpairedRead,ONTAlignmentModel>>(
    BulkAlnLibT<UnpairedRead,ONTAlignmentModel>& readExp, SalmonOpts& sopt,
    std::function<bool(const std::vector<double>&)>& writeBootstrap,
    double relDiffTolerance, uint32_t maxIter);

template bool
CollapsedEMOptimizer::gatherBootstraps<BulkAlnLibT<ReadPair,AlignmentModel>>(
    BulkAlnLibT<ReadPair,AlignmentModel>& readExp, SalmonOpts& sopt,
    std::function<bool(const std::vector<double>&)>& writeBootstrap,
    double relDiffTolerance, uint32_t maxIter);
// Unused / old

template
void EMUpdate_forAugmented<std::vector<std::atomic<double>>>(std::vector<std::vector<uint32_t>>& txpGroupLabels,
                                                 std::vector<std::vector<double>>& txpGroupCombinedWeights,
                                                 const std::vector<double>& txpGroupCounts,
                                                 std::vector<bool>& augmented_tids,
                                                 std::vector<bool>& augmented_eqIDs,
                                                 const std::vector<std::atomic<double>>& alphaIn,
                                                 std::vector<std::atomic<double>>& alphaOut);
