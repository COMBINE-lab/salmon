#include <algorithm>
#include <atomic>
#include <cstring>
#include <unordered_map>
#include <vector>
#include <exception>

#include "oneapi/tbb/task_arena.h"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"
#include "oneapi/tbb/parallel_for_each.h"
#include "oneapi/tbb/parallel_reduce.h"
#include "oneapi/tbb/partitioner.h"

//#include "fastapprox.h"
#include <boost/math/special_functions/digamma.hpp>

// C++ string formatting library
#include <spdlog/fmt/fmt.h>

#include "Eigen/Dense"
#include "salmon/vendor/cuckoohash_map.hh"
#include "salmon/vendor/unordered_dense.hpp"

#include "salmon/internal/alignment/AlignmentLibrary.hpp"
#include "salmon/internal/util/FmtCompat.hpp"
#include "salmon/internal/output/BootstrapWriter.hpp"
#include "salmon/internal/inference/CollapsedEMOptimizer.hpp"
#include "salmon/internal/inference/MultinomialSampler.hpp"
#include "salmon/internal/quant/ReadExperiment.hpp"
#include "salmon/internal/alignment/ReadPair.hpp"
#include "salmon/internal/util/SalmonMath.hpp"
#include "salmon/internal/model/Transcript.hpp"
#include "salmon/internal/model/TranscriptGroup.hpp"
#include "salmon/internal/alignment/UnpairedRead.hpp"
#include "salmon/internal/inference/EMUtils.hpp"

using BlockedIndexRange = oneapi::tbb::blocked_range<size_t>;

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

/*
 * Use the "standard" EM algorithm over equivalence
 * classes to estimate the latent variables (alphaOut)
 * given the current estimates (alphaIn).
 */
// Deterministically reduce the per-shard private accumulators into alphaOut.
// Each transcript is summed over the shards in a fixed order and written
// exactly once (no atomics, no contention), so the result is identical
// regardless of thread scheduling -- this mirrors the Rust optimizer's
// single-writer M-step reduction (reduce_shards).
static inline void reduceShards_(oneapi::tbb::task_arena& arena,
                                 std::vector<std::vector<double>>& shards,
                                 CollapsedEMOptimizer::VecType& alphaOut) {
  size_t M = alphaOut.size();
  size_t nshards = shards.size();
  arena.execute([&] {
    oneapi::tbb::parallel_for(
        BlockedIndexRange(size_t(0), M),
        [&shards, &alphaOut, nshards](const BlockedIndexRange& range) -> void {
          for (auto i : boost::irange(range.begin(), range.end())) {
            double acc = 0.0;
            for (size_t s = 0; s < nshards; ++s) {
              acc += shards[s][i];
            }
            alphaOut[i].store(acc);
          }
        });
  });
}

// When the offline optimizer is seeded with the online (blended) estimates,
// mutually non-identifiable transcripts -- those appearing in exactly the same
// equivalence classes with the same combined weights (e.g. exact-duplicate
// transcripts) -- can end up with slightly different seeds purely from
// hog-wild online-phase noise, which the VBEM then amplifies into an arbitrary
// winner-take-all split. This pass groups such transcripts by an
// order-independent signature over their (eqID, combinedWeight) memberships and
// replaces each group's seed with the group mean, so the (deterministic)
// offline update keeps them tied. It is applied ONLY on the --useOnlineSeed
// path; the default uniform seed is already symmetric and needs no equalization.
template <typename EQVecT>
size_t equalizeNonIdentifiableSeed_(EQVecT& eqVec,
                                    CollapsedEMOptimizer::VecType& alphas,
                                    size_t numTxps) {
  constexpr uint64_t kFnvOffset = 1469598103934665603ull;
  constexpr uint64_t kFnvPrime = 1099511628211ull;
  auto mix = [](uint64_t h, uint64_t x) -> uint64_t {
    h ^= x;
    h *= kFnvPrime;
    return h;
  };
  std::vector<uint64_t> sig(numTxps, kFnvOffset);
  for (size_t eqID = 0; eqID < eqVec.size(); ++eqID) {
    auto& kv = eqVec[eqID];
    if (!kv.first.valid) {
      continue;
    }
    const std::vector<uint32_t>& txps = kv.first.txps;
    const auto& auxs = kv.second.combinedWeights;
    for (size_t i = 0; i < txps.size(); ++i) {
      uint64_t wbits = 0;
      double w = auxs[i];
      std::memcpy(&wbits, &w, sizeof(wbits));
      uint64_t h = sig[txps[i]];
      h = mix(h, static_cast<uint64_t>(eqID) + 0x9e3779b97f4a7c15ull);
      h = mix(h, wbits);
      sig[txps[i]] = h;
    }
  }
  ankerl::unordered_dense::map<uint64_t, std::vector<uint32_t>> groups;
  groups.reserve(numTxps);
  for (size_t t = 0; t < numTxps; ++t) {
    groups[sig[t]].push_back(static_cast<uint32_t>(t));
  }
  size_t nGroups = 0;
  for (auto& g : groups) {
    std::vector<uint32_t>& members = g.second;
    if (members.size() < 2) {
      continue;
    }
    double s = 0.0;
    for (auto t : members) {
      s += alphas[t].load();
    }
    double avg = s / static_cast<double>(members.size());
    for (auto t : members) {
      alphas[t].store(avg);
    }
    ++nGroups;
  }
  return nGroups;
}

template <typename EQVecT>
void EMUpdate_(oneapi::tbb::task_arena& arena,
               EQVecT& eqVec,
               std::vector<double>& priorAlphas,
               const CollapsedEMOptimizer::VecType& alphaIn,
               CollapsedEMOptimizer::VecType& alphaOut,
               std::vector<std::vector<double>>& shards) {

  assert(alphaIn.size() == alphaOut.size());

  size_t nshards = shards.size();
  size_t nclasses = eqVec.size();
  size_t chunk = (nclasses + nshards - 1) / nshards;

  // Accumulate into per-shard private buffers using plain (non-atomic) adds.
  // Each shard owns a fixed, contiguous range of equivalence classes, so the
  // work partition -- and the within-shard summation order -- is independent
  // of how TBB schedules the shard tasks across threads.
  arena.execute([&] {
    oneapi::tbb::parallel_for(size_t(0), nshards, [&](size_t s) {
      std::vector<double>& buf = shards[s];
      std::fill(buf.begin(), buf.end(), 0.0);
      size_t start = s * chunk;
      size_t end = std::min(start + chunk, nclasses);
      for (size_t eqID = start; eqID < end; ++eqID) {
        auto& kv = eqVec[eqID];
        uint64_t count = kv.second.count;
        const TranscriptGroup& tgroup = kv.first;
        if (tgroup.valid) {
          const std::vector<uint32_t>& txps = tgroup.txps;
          const auto& auxs = kv.second.combinedWeights;
          size_t groupSize = kv.second.weights.size(); // txps.size();
          if (BOOST_LIKELY(groupSize > 1)) {
            double denom = 0.0;
            for (size_t i = 0; i < groupSize; ++i) {
              denom += (alphaIn[txps[i]]) * auxs[i];
            }
            if (denom > ::minEQClassWeight) {
              double invDenom = count / denom;
              for (size_t i = 0; i < groupSize; ++i) {
                double v = (alphaIn[txps[i]]) * auxs[i];
                if (!std::isnan(v)) {
                  buf[txps[i]] += v * invDenom;
                }
              }
            }
          } else {
            buf[txps.front()] += count;
          }
        }
      }
    });
  });

  reduceShards_(arena, shards, alphaOut);
}

/*
 * Use the Variational Bayesian EM algorithm over equivalence
 * classes to estimate the latent variables (alphaOut)
 * given the current estimates (alphaIn).
 */
template <typename EQVecT>
void VBEMUpdate_(oneapi::tbb::task_arena& arena,
                 EQVecT& eqVec,
                 std::vector<double>& priorAlphas,
                 const CollapsedEMOptimizer::VecType& alphaIn,
                 CollapsedEMOptimizer::VecType& alphaOut,
                 CollapsedEMOptimizer::VecType& expTheta,
                 std::vector<std::vector<double>>& shards) {

  assert(alphaIn.size() == alphaOut.size());
  size_t M = alphaIn.size();
  double alphaSum = {0.0};
  for (size_t i = 0; i < M; ++i) {
    alphaSum += alphaIn[i] + priorAlphas[i];
  }

  double logNorm = boost::math::digamma(alphaSum);

  arena.execute([&]{
  oneapi::tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(priorAlphas.size())),
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
  });

  size_t nshards = shards.size();
  size_t nclasses = eqVec.size();
  size_t chunk = (nclasses + nshards - 1) / nshards;

  // Accumulate into per-shard private buffers (see EMUpdate_ for rationale).
  arena.execute([&] {
    oneapi::tbb::parallel_for(size_t(0), nshards, [&](size_t s) {
      std::vector<double>& buf = shards[s];
      std::fill(buf.begin(), buf.end(), 0.0);
      size_t start = s * chunk;
      size_t end = std::min(start + chunk, nclasses);
      for (size_t eqID = start; eqID < end; ++eqID) {
        auto& kv = eqVec[eqID];
        uint64_t count = kv.second.count;
        const TranscriptGroup& tgroup = kv.first;
        if (tgroup.valid) {
          const std::vector<uint32_t>& txps = tgroup.txps;
          const auto& auxs = kv.second.combinedWeights;
          size_t groupSize = kv.second.weights.size(); // txps.size();
          if (BOOST_LIKELY(groupSize > 1)) {
            double denom = 0.0;
            for (size_t i = 0; i < groupSize; ++i) {
              if (expTheta[txps[i]] > 0.0) {
                denom += expTheta[txps[i]] * auxs[i];
              }
            }
            if (denom > ::minEQClassWeight) {
              double invDenom = count / denom;
              for (size_t i = 0; i < groupSize; ++i) {
                if (expTheta[txps[i]] > 0.0) {
                  buf[txps[i]] += expTheta[txps[i]] * auxs[i] * invDenom;
                }
              }
            }
          } else {
            buf[txps.front()] += count;
          }
        }
      }
    });
  });

  reduceShards_(arena, shards, alphaOut);
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

        jointLog->info("{}", errstream.str());
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
    double relDiffTolerance, uint32_t maxIter) {

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
  std::discrete_distribution<uint64_t> csamp(sampleWeights.begin(),
                                             sampleWeights.end());
  while (bsNum++ < numBootstraps) {
    csamp.reset();

    for (size_t sc = 0; sc < sampCounts.size(); ++sc) {
      sampCounts[sc] = 0;
    }
    for (size_t fn = 0; fn < totalNumFrags; ++fn) {
      ++sampCounts[csamp(gen)];
    }
    // Do a new bootstrap
    // msamp(sampCounts.begin(), totalNumFrags, numClasses,
    // sampleWeights.begin());

    for (size_t i = 0; i < transcripts.size(); ++i) {
      alphas[i] =
          transcripts[i].getActive() ? uniformTxpWeight * totalNumFrags : 0.0;
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

    while (itNum < minIter or (itNum < maxIter and !converged)) {

      if (useVBEM) {
        VBEMUpdate_(txpGroups, txpGroupCombinedWeights, sampCounts, 
                    priorAlphas, alphas, alphasPrime, expTheta);
      } else {
        EMUpdate_(txpGroups, txpGroupCombinedWeights, sampCounts, 
                  alphas, alphasPrime);
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

  jointLog->info("Will draw {} bootstrap samples", numBootstraps);
  jointLog->info("Optimizing over {} equivalence classes", eqVec.size());

  double totalNumFrags{static_cast<double>(numMappedFrags)};
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
      
      txpGroups.push_back(txps);
      // Convert to non-atomic
      txpGroupCombinedWeights.emplace_back(auxs.begin(), auxs.end());
      origCounts.push_back(count);
      totalCount += count;
    }
  }

  double floatCount = totalCount;
  std::vector<double> samplingWeights(txpGroups.size(), 0.0);
  for (size_t i = 0; i < origCounts.size(); ++i) {
    samplingWeights[i] = origCounts[i] / floatCount;
  }

  size_t numWorkerThreads{1};
  if (sopt.numThreads > 1 and numBootstraps > 1) {
    numWorkerThreads = std::min(sopt.numThreads - 1, numBootstraps - 1);
  }

  std::atomic<uint32_t> bsCounter{0};
  std::vector<std::thread> workerThreads;
  for (size_t tn = 0; tn < numWorkerThreads; ++tn) {
    workerThreads.emplace_back(
        doBootstrap, std::ref(txpGroups), std::ref(txpGroupCombinedWeights),
        std::ref(transcripts), std::ref(effLens), std::ref(samplingWeights), std::ref(origCounts),
        totalCount, numMappedFrags, scale, std::ref(bsCounter), std::ref(sopt),
        std::ref(priorAlphas), std::ref(writeBootstrap), relDiffTolerance,
        maxIter);
  }

  for (auto& t : workerThreads) {
    t.join();
  }
  return true;
}

template <typename EQVecT>
void updateEqClassWeights(
    oneapi::tbb::task_arena& arena,
    EQVecT& eqVec,
    Eigen::VectorXd& effLens) {
  
  arena.execute([&]{
  oneapi::tbb::parallel_for(
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
  });
}

template <typename ExpT>
bool CollapsedEMOptimizer::optimize(ExpT& readExp, SalmonOpts& sopt,
                                    double relDiffTolerance, uint32_t maxIter) {

  oneapi::tbb::task_arena arena(sopt.numThreads);

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
  // Seed the offline optimizer. By DEFAULT we initialize *uniformly*: the
  // online (blended) estimates carry hog-wild, multithreaded run-to-run noise
  // that the VBEM update can amplify into arbitrary winner-take-all splits of
  // low-abundance / non-identifiable transcripts (issues #1008, #1011), with
  // no measurable accuracy benefit over a uniform seed. --useOnlineSeed
  // restores the online seed (e.g. for faster convergence); the metagenomic
  // and alternative-init modes also supply their own informative seed.
  // --initUniform always forces uniform.
  bool wantBlended =
      (sopt.useOnlineSeed or metaGenomeMode or altInitMode) and not sopt.initUniform;
  if (not wantBlended) {
    for (size_t i = 0; i < alphas.size(); ++i) {
      alphas[i].store(100.0);
      alphasPrime[i] = 1.0;
    }
  } else { // initialize with a linear combination of the online and uniform alphas
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
  arena.execute([&]{
  oneapi::tbb::parallel_for(
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
  });

  auto numRemoved =
      markDegenerateClasses(eqVec, alphas, available, sopt.jointLog);
  sopt.jointLog->info("Marked {} weighted equivalence classes as degenerate",
                      numRemoved);

  // On the --useOnlineSeed path, neutralize hog-wild seed asymmetry between
  // mutually non-identifiable transcripts (exact duplicates and the like) so
  // the offline optimizer does not split them arbitrarily. The default uniform
  // seed is already symmetric, so this is skipped there.
  if (sopt.useOnlineSeed and not sopt.initUniform) {
    size_t nEqGroups =
        equalizeNonIdentifiableSeed_(eqVec, alphas, transcripts.size());
    sopt.jointLog->info(
        "useOnlineSeed: equalized the online seed across {} non-identifiable "
        "transcript group(s)",
        nEqGroups);
  }

  size_t itNum{0};

  // Per-shard private accumulators for the (V)BEM M-step, allocated once and
  // reused across iterations. The M-step sums each shard's contributions into
  // alphaOut in a fixed order (reduceShards_), making the parallel update
  // deterministic and contention-free instead of using atomic CAS adds.
  size_t numShards = std::max<size_t>(1, static_cast<size_t>(sopt.numThreads));
  std::vector<std::vector<double>> mStepShards(
      numShards, std::vector<double>(transcripts.size(), 0.0));

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
          "iteration {}, adjusting effective lengths to account for biases",
          itNum);
      effLens = salmon::utils::updateEffectiveLengths(arena, sopt, readExp, effLens,
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
      updateEqClassWeights(arena, eqVec, effLens);
      needBias = false;

      if ( sopt.eqClassMode ) {
        jointLog->error("Eqclass Mode should not be performing bias correction");
        jointLog->flush();
        exit(1);
      }
    }

    if (useVBEM) {
      VBEMUpdate_(arena, eqVec, priorAlphas, alphas,
                  alphasPrime, expTheta, mStepShards);
    } else {
      /*
      if (itNum > 0 and (itNum % 250 == 0)) {
        for (size_t i = 0; i < transcripts.size(); ++i) {
      	  if (alphas[i] < 1.0) { alphas[i] = 0.0; }
      	}
      }
      */

      EMUpdate_(arena, eqVec, priorAlphas, alphas, alphasPrime, mStepShards);
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
      jointLog->info("iteration = {} | max rel diff. = {}", itNum, maxRelDiff);
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

  jointLog->info("iteration = {} | max rel diff. = {}", itNum, maxRelDiff);

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
