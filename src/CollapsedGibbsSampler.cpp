#include <atomic>
#include <mutex>
#include <random>
#include <thread>
#include <unordered_map>
#include <vector>

#include "tbb/blocked_range.h"
#include "tbb/combinable.h"
#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_reduce.h"
#include "tbb/partitioner.h"
// <-- deprecated in TBB --> #include "tbb/task_scheduler_init.h"
#include "tbb/global_control.h"

//#include "fastapprox.h"
#include <boost/filesystem.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
// PCG Random number generator
#include "pcg_random.hpp"

// C++ string formatting library
#include "spdlog/fmt/fmt.h"

#include "Eigen/Dense"
#include "cuckoohash_map.hh"

#include "AlignmentLibrary.hpp"
#include "BootstrapWriter.hpp"
#include "CollapsedGibbsSampler.hpp"
#include "MultinomialSampler.hpp"
#include "ReadExperiment.hpp"
#include "ReadPair.hpp"
#include "SalmonMath.hpp"
#include "Transcript.hpp"
#include "TranscriptGroup.hpp"
#include "UnpairedRead.hpp"
#include "ezETAProgressBar.hpp"

using BlockedIndexRange = tbb::blocked_range<size_t>;

// intelligently chosen value adopted from
// https://github.com/pachterlab/kallisto/blob/master/src/EMAlgorithm.h#L18
constexpr double minEQClassWeight = std::numeric_limits<double>::denorm_min();
constexpr double minWeight = std::numeric_limits<double>::denorm_min();

/**
 * http://codereview.stackexchange.com/questions/106773/dividing-a-range-into-n-sub-ranges
 */
template <typename Iterator>
std::vector<std::pair<Iterator, Iterator>>
divide_work(Iterator begin, Iterator end, std::size_t n) {
  std::vector<std::pair<Iterator, Iterator>> ranges;
  if (n == 0)
    return ranges;
  ranges.reserve(n);

  auto dist = std::distance(begin, end);
  n = std::min<size_t>(n, dist);
  auto chunk = dist / n;
  auto remainder = dist % n;

  for (size_t i = 0; i < n - 1; ++i) {
    auto next_end = std::next(begin, chunk + (remainder ? 1 : 0));
    ranges.emplace_back(begin, next_end);

    begin = next_end;
    if (remainder)
      remainder -= 1;
  }

  // last chunk
  ranges.emplace_back(begin, end);
  return ranges;
}

/**
 * This non-collapsed Gibbs step is largely inspired by the method first
 * introduced by  Turro et al. [1].  Given the current estimates `txpCount` of
 *the read count for each transcript,  the mean transcript fractions are sampled
 *from a Gamma distribution ~ Gam( prior[i] + txpCount[i], \Beta + effLens[i]).
 *Then, given these transcript fractions,  The reads are re-assigned within
 *each equivalence class by sampling from a multinomial distributed according
 *to these means.
 *
 * [1] Haplotype and isoform specific expression estimation using multi-mapping
 *RNA-seq reads. Turro E, Su S-Y, Goncalves A, Coin L, Richardson S and Lewin A.
 * Genome Biology, 2011 Feb; 12:R13.  doi: 10.1186/gb-2011-12-2-r13.
 **/
template <typename EQVecT>
void sampleRoundNonCollapsedMultithreaded_(
    EQVecT& eqVec,
    /*std::vector<bool>& active,*/ std::vector<uint32_t>& activeList,
    /*std::vector<uint64_t>& countMap,*/ std::vector<double>& probMap,
    std::vector<double>& muGlobal, Eigen::VectorXd& effLens,
    const std::vector<double>& priorAlphas, std::vector<double>& txpCount,
    std::vector<uint32_t>& offsetMap,
    bool noGammaDraw) {

  // generate coeff for \mu from \alpha and \effLens
  double beta = 0.1;
  double norm = 0.0;

  // Sample the transcript fractions \mu from a gamma distribution, and
  // reset txpCounts to zero for each transcript.
  typedef tbb::enumerable_thread_specific<pcg32_unique> GeneratorType;
  auto getGenerator = []() -> pcg32_unique {
    // why this mess below?  see SalmonUtils.hpp : get_random_device() for more
    // details.
    #if defined(__linux) && defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128
      return pcg32_unique(pcg_extras::seed_seq_from<std::random_device>("/dev/urandom"));
    #else
      return pcg32_unique(pcg_extras::seed_seq_from<std::random_device>());
    #endif  // defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128
  };
  GeneratorType localGenerator(getGenerator);

  // Compute the mu to be used in the equiv class resampling
  // If we are doing a gamma draw (including shot-noise)
  if (noGammaDraw) {
   tbb::parallel_for(
      BlockedIndexRange(
          size_t(0), size_t(activeList.size())), // 1024 is grainsize, use only
                                                 // with simple_partitioner
      [&, beta](const BlockedIndexRange& range) -> void {
        for (auto activeIdx : boost::irange(range.begin(), range.end())) {
          auto i = activeList[activeIdx];
          double ci = static_cast<double>(txpCount[i] + priorAlphas[i]);
          muGlobal[i] = ci / effLens(i);
          txpCount[i] = 0.0;
        }
      });
   
  } else {
   tbb::parallel_for(
      BlockedIndexRange(
          size_t(0), size_t(activeList.size())), // 1024 is grainsize, use only
                                                 // with simple_partitioner
      [&, beta](const BlockedIndexRange& range) -> void {
        GeneratorType::reference gen = localGenerator.local();
        for (auto activeIdx : boost::irange(range.begin(), range.end())) {
          auto i = activeList[activeIdx];
          double ci = static_cast<double>(txpCount[i] + priorAlphas[i]);
          std::gamma_distribution<double> d(ci, 1.0 / (beta + effLens(i)));
          muGlobal[i] = d(gen);
          txpCount[i] = 0.0;
          /** DEBUG
          if (std::isnan(muGlobal[i]) or std::isinf(muGlobal[i])) {
            std::cerr << "txpCount = " << txpCount[i] << ", prior = " <<
          priorAlphas[i] << ", alpha = " << ci << ", beta = " << (1.0 / (beta +
          effLens(i))) << ", mu = " << muGlobal[i] << "\n"; std::exit(1);
          }
          **/
        }
      });
  }

  /**
   * These will store "thread local" parameters
   * for the threads doing local sampling of equivalence class counts.
   */
  class CombineableTxpCounts {
  public:
    CombineableTxpCounts(uint32_t numTxp) : txpCount(numTxp, 0) {
      // why this mess below?  see SalmonUtils.hpp : get_random_device() for more
      // details.
      #if defined(__linux) && defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128
        gen.reset(
            new pcg32_unique(pcg_extras::seed_seq_from<std::random_device>("/dev/urandom")));
      #else
        gen.reset(
          new pcg32_unique(pcg_extras::seed_seq_from<std::random_device>()));
      #endif  // defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128
    }
    std::vector<int> txpCount;
    std::unique_ptr<pcg32_unique> gen{nullptr};
  };
  tbb::combinable<CombineableTxpCounts> combineableCounts(txpCount.size());

  std::mutex writeMut;
  // resample within each equivalence class
  tbb::parallel_for(
      BlockedIndexRange(size_t(0), size_t(eqVec.size())),
      [&](const BlockedIndexRange& range) -> void {

        auto& txpCountLoc = combineableCounts.local().txpCount;
        auto& gen = *(combineableCounts.local().gen.get());
        for (auto eqid : boost::irange(range.begin(), range.end())) {
          auto& eqClass = eqVec[eqid];
          size_t offset = offsetMap[eqid];

          // get total number of reads for an equivalence class
          uint64_t classCount = eqClass.second.count;

          // for each transcript in this class
          const TranscriptGroup& tgroup = eqClass.first;
          const size_t groupSize =
              eqClass.second.weights.size(); // tgroup.txps.size();
          if (tgroup.valid) {
            const std::vector<uint32_t>& txps = tgroup.txps;
            const auto& auxs = eqClass.second.combinedWeights;
            const auto& weights = eqClass.second.weights;

            double denom = 0.0;
            // If this is a single-transcript group,
            // then it gets the full count --- otherwise,
            // sample!
            if (BOOST_LIKELY(groupSize > 1)) {
              // For each transcript in the group
              double muSum = 0.0;
              for (size_t i = 0; i < groupSize; ++i) {
                auto tid = txps[i];
                size_t gi = offset + i;
                probMap[gi] = (1000.0 * muGlobal[tid]) * weights[i];
                muSum += probMap[gi];
                denom += probMap[gi];
              }

              if (denom <= ::minEQClassWeight) {
                {
                  std::lock_guard<std::mutex> lg(writeMut);
                  std::cerr
                      << "[WARNING] eq class denom was too small : denom = "
                      << denom << ", numReads = " << classCount
                      << ". Distributing reads evenly for this class\n";
                }

                denom = 0.0;
                muSum = 0.0;
                for (size_t i = 0; i < groupSize; ++i) {
                  auto tid = txps[i];
                  size_t gi = offset + i;
                  probMap[gi] = 1.0 / effLens(tid);
                  muSum += probMap[gi];
                  denom += probMap[gi];
                }

                // If it's still too small --- divide evenly
                if (denom <= ::minEQClassWeight) {
                  for (size_t i = 0; i < groupSize; ++i) {
                    auto tid = txps[i];
                    size_t gi = offset + i;
                    probMap[gi] = 1.0;
                  }
                  denom = groupSize;
                  muSum = groupSize;
                }
              }

              if (denom > ::minEQClassWeight) {
                // Local multinomial
                std::discrete_distribution<int> dist(probMap.begin() + offset,
                                                     probMap.begin() + offset +
                                                         groupSize);
                for (size_t s = 0; s < classCount; ++s) {
                  auto ind = dist(gen);
                  ++txpCountLoc[txps[ind]];
                }
              }
            } // do nothing if group size less than 2
            else {
              auto tid = txps[0];
              txpCountLoc[tid] += static_cast<int>(classCount);
            }
          } // valid group
        }   // loop over all eq classes
      });

  auto combineCounts = [&txpCount](const CombineableTxpCounts& p) -> void {
    for (size_t i = 0; i < txpCount.size(); ++i) {
      txpCount[i] += static_cast<double>(p.txpCount[i]);
    }
  };
  combineableCounts.combine_each(combineCounts);
}

CollapsedGibbsSampler::CollapsedGibbsSampler() {}

class DistStats {
public:
  DistStats()
      : meanVal(0.0), minVal(std::numeric_limits<double>::max()), maxVal(0.0) {}
  double meanVal;
  double minVal;
  double maxVal;
};

/**
 *  Populate the prior parameters for the VBEM
 *  Note: effLens *must* be valid before calling this function.
 */
// Get rid of redundancy of this function
std::vector<double> populatePriorAlphasGibbs_(
    std::vector<Transcript>& transcripts, // transcripts
    Eigen::VectorXd& effLens,             // current effective length estimate
    double priorValue,      // the per-nucleotide prior value to use
    bool perTranscriptPrior // true if prior is per-txp, else per-nucleotide
) {
  // start out with the per-txp prior
  std::vector<double> priorAlphas(transcripts.size(), priorValue);

  // If the prior is per-nucleotide (default, then we need a potentially
  // different
  // value for each transcript based on its length).
  if (!perTranscriptPrior) {
    for (size_t i = 0; i < transcripts.size(); ++i) {
      double ml = std::max(1.0, effLens(i));
      priorAlphas[i] = priorValue * ml;
    }
  }
  return priorAlphas;
}

template <typename ExpT>
bool CollapsedGibbsSampler::sample(
    ExpT& readExp, SalmonOpts& sopt,
    std::function<bool(const std::vector<double>&)>& writeBootstrap,
    uint32_t numSamples) {

  namespace bfs = boost::filesystem;
  auto& jointLog = sopt.jointLog;
  
  // <-- deprecated in TBB --> tbb::task_scheduler_init tbbScheduler(sopt.numThreads);
  tbb::global_control c(tbb::global_control::max_allowed_parallelism, sopt.numThreads);

  std::vector<Transcript>& transcripts = readExp.transcripts();

  // Fill in the effective length vector
  Eigen::VectorXd effLens(transcripts.size());

  auto& eqBuilder = readExp.equivalenceClassBuilder();
  auto& eqVec = eqBuilder.eqVec();

  using VecT = CollapsedGibbsSampler::VecType;

  std::vector<std::vector<int>> allSamples(
      numSamples, std::vector<int>(transcripts.size(), 0));

  std::vector<double> alphasIn(transcripts.size(), 0.0);
  std::vector<double> alphasInit(transcripts.size(), 0.0);

  bool useScaledCounts = (!sopt.useQuasi and !sopt.allowOrphans);
  auto numMappedFragments = (useScaledCounts) ? readExp.upperBoundHits()
                                              : readExp.numMappedFragments();
  uint32_t numInternalRounds = sopt.thinningFactor;
  size_t numTranscripts{transcripts.size()};

  for (size_t i = 0; i < transcripts.size(); ++i) {
    auto& txp = transcripts[i];
    alphasIn[i] = txp.projectedCounts;
    alphasInit[i] = txp.projectedCounts;
    effLens(i) = txp.EffectiveLength;
  }

  bool perTranscriptPrior = (sopt.useVBOpt) ? sopt.perTranscriptPrior : true;
  // for Gibbs sampling, we really want a notion of the uncertainty, and we 
  // should steer away from sparsity more than we do when computing an MLE/MAP
  // we set the minimum prior to 1 or 1e-3 if it is per-nucleotide
  double prior = 1e-3; // prior to use if main algo was EM
  if (sopt.useVBOpt) {
    if (perTranscriptPrior) {
      prior = (sopt.vbPrior  < 1.0) ? 1.0 : sopt.vbPrior;
    } else { // per-nucleotide
      prior = (sopt.vbPrior < 1e-3) ? 1e-3 : sopt.vbPrior;
    } 
  }
  double priorValue = prior;
  std::vector<double> priorAlphas = populatePriorAlphasGibbs_(
      transcripts, effLens, priorValue, perTranscriptPrior);
  /** DEBUG
  for (size_t i = 0; i < priorAlphas.size(); ++i) {
    auto& v = priorAlphas[i];
    if (!std::isfinite(v)) {
      std::cerr << "prior for transcript " << i << " is " << v << ", eff length
  = " << effLens(i) << "\n";
    }
  }
  **/

  std::vector<bool> active(numTranscripts, false);
  size_t countMapSize{0};
  std::vector<uint32_t> offsetMap(eqVec.size(), 0);
  for (size_t i = 0; i < eqVec.size(); ++i) {
    if (eqVec[i].first.valid) {
      size_t numTranscripts = eqBuilder.getNumTranscriptsForClass(i);
      countMapSize += numTranscripts;
          // eqVec[i].second.weights.size(); 
      const auto& txpVec = eqVec[i].first.txps;
      for (size_t tidx = 0; tidx < numTranscripts; ++tidx) {
        auto t = txpVec[tidx];
        active[t] = true;
      }
      if (i < eqVec.size() - 1) {
        offsetMap[i + 1] = countMapSize;
      }
    }
  }

  std::vector<uint32_t> activeList;
  activeList.reserve(numTranscripts);
  for (size_t i = 0; i < numTranscripts; ++i) {
    if (active[i]) {
      activeList.push_back(i);
    } else {
      alphasIn[i] = 0.0;
      alphasInit[i] = 0.0;
    }
  }

  // will hold estimated counts
  std::vector<double> alphas(numTranscripts, 0.0);
  std::vector<double> mu(numTranscripts, 0.0);
  std::vector<uint64_t> countMap(countMapSize, 0);
  std::vector<double> probMap(countMapSize, 0.0);

  /*
  std::random_device rd;
  MultinomialSampler ms(rd);
  initCountMap_(eqVec, transcripts, alphasIn, priorAlphas, ms, countMap,
  probMap, effLens, allSamples[0]);
  */

  uint32_t nchains{1};
  if (numSamples >= 50) {
    nchains = 2;
  }
  if (numSamples >= 100) {
    nchains = 4;
  }
  if (numSamples >= 200) {
    nchains = 8;
  }

  std::vector<uint32_t> newChainIter{0};
  if (nchains > 1) {
    auto step = numSamples / nchains;
    for (size_t i = 1; i < nchains; ++i) {
      newChainIter.push_back(i * step);
    }
  }

  auto nextChainStart = newChainIter.begin();

  // For each sample this thread should generate
  std::unique_ptr<ez::ezETAProgressBar> pbar{nullptr};
  if (!sopt.quiet) {
    pbar.reset(new ez::ezETAProgressBar(numSamples));
    pbar->start();
  }
  //bool isFirstSample{true};
  for (size_t sampleID = 0; sampleID < numSamples; ++sampleID) {
    if (pbar) {
      ++(*pbar);
    }
    // If we should start a new chain here, then do it!
    if (nextChainStart < newChainIter.end() and sampleID == *nextChainStart) {
      alphasIn = alphasInit;
      ++nextChainStart;
    }
    /*
      if (!isFirstSample) {
          // the counts start at what they were last round.
        allSamples[sampleID] = allSamples[sampleID-1];
      }
      */

    // Thin the chain by a factor of (numInternalRounds)
    for (size_t i = 0; i < numInternalRounds; ++i) {
      sampleRoundNonCollapsedMultithreaded_(
          eqVec,      // encodes equivalence classes
          /*active,     // the set of active transcripts*/
          activeList, // the list of active transcript ids
          /*countMap,   // the count of reads in each eq coming from each eq class*/
          probMap, // the probability of reads in each eq class coming from each
                   // txp
          mu,      // transcript fractions
          effLens, // the effective transcript lengths
          priorAlphas, // the prior transcript counts
          alphasIn, // [input/output param] the (hard) fragment counts per txp
                    // from the previous iteration
          offsetMap, // where the information begins for each equivalence class
          sopt.noGammaDraw      // true if we should skip the Gamma draw, false otherwise
      );
    }

    if (sopt.dontExtrapolateCounts) {
      alphas = alphasIn;
    } else {
      double denom{0.0};
      for (size_t tn = 0; tn < numTranscripts; ++tn) {
        denom += mu[tn] * effLens[tn];
      }
      double scale = numMappedFragments / denom;
      double asum = {0.0};

      // A read cutoff for a txp to be present, adopted from Bray et al. 2016
      double minAlpha = 1e-8;
      for (size_t tn = 0; tn < numTranscripts; ++tn) {
        alphas[tn] = (mu[tn] * effLens[tn]) * scale;
        alphas[tn] = (alphas[tn] > minAlpha) ? alphas[tn] : 0.0;
        asum += alphas[tn];
      }
    }
    writeBootstrap(alphas);
    //isFirstSample = false;
  }
  return true;
}

/*
void initCountMap_(
    std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
    std::vector<Transcript>& transcriptsIn, const std::vector<double>& alphasIn,
    const std::vector<double>& priorAlphas, MultinomialSampler& msamp,
    std::vector<uint64_t>& countMap, std::vector<double>& probMap,
    Eigen::VectorXd& effLens, std::vector<int>& txpCounts) {

  size_t offset{0};
  for (auto& eqClass : eqVec) {
    uint64_t classCount = eqClass.second.count;

    // for each transcript in this class
    const TranscriptGroup& tgroup = eqClass.first;
    const size_t groupSize = tgroup.txps.size();
    if (tgroup.valid) {
      const std::vector<uint32_t>& txps = tgroup.txps;
      const auto& auxs = eqClass.second.combinedWeights;

      double denom = 0.0;
      if (BOOST_LIKELY(groupSize > 1)) {

        for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          auto aux = auxs[i];
          denom += (priorAlphas[tid] + alphasIn[tid]) * aux;
          countMap[offset + i] = 0;
        }

        if (denom > ::minEQClassWeight) {
          // Get the multinomial probabilities
          double norm = 1.0 / denom;
          for (size_t i = 0; i < groupSize; ++i) {
            auto tid = txps[i];
            auto aux = auxs[i];
            probMap[offset + i] =
                norm * ((priorAlphas[tid] + alphasIn[tid]) * aux);
          }

          // re-sample
          msamp(countMap.begin() + offset, classCount, groupSize,
                probMap.begin() + offset);
        }
      } else {
        countMap[offset] = classCount;
      }

      for (size_t i = 0; i < groupSize; ++i) {
        auto tid = txps[i];
        txpCounts[tid] += countMap[offset + i];
      }

      offset += groupSize;
    } // valid group
  }   // loop over all eq classes
}

//
 // This non-collapsed Gibbs step is largely inspired by the method first
 //introduced by
 // Turro et al. [1].  Given the current estimates `txpCount` of the read count
 //for each transcript,
 // the mean transcript fractions are sampled from a Gamma distribution
 // ~ Gam( prior[i] + txpCount[i], \Beta + effLens[i]).  Then, given these
 //transcript fractions,
 // The reads are re-assigned within each equivalence class by sampling from a
 //multinomial
 // distributed according to these means.
 //
 // [1] Haplotype and isoform specific expression estimation using multi-mapping
 //RNA-seq reads.
 // Turro E, Su S-Y, Goncalves A, Coin L, Richardson S and Lewin A. Genome
 //Biology, 2011 Feb; 12:R13.
 // doi: 10.1186/gb-2011-12-2-r13.
 //
void sampleRoundNonCollapsed_(
    std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
    std::vector<uint64_t>& countMap, std::vector<double>& probMap,
    Eigen::VectorXd& effLens, const std::vector<double>& priorAlphas,
    std::vector<int>& txpCount, MultinomialSampler& msamp) {
  std::random_device rd;
  std::mt19937 gen(rd());
  // offset for 2d to 1d count map
  size_t offset{0};

  // retain original txp count
  std::vector<int> origTxpCount = txpCount;

  // reset txpCounts to zero
  std::fill(txpCount.begin(), txpCount.end(), 0);

  // generate norm. coeff for \mu from \alpha (countMap)
  std::vector<double> muGlobal(txpCount.size(), 0.0);
  double beta = 0.1;
  double norm = 0.0;
  for (size_t i = 0; i < origTxpCount.size(); ++i) {
    std::gamma_distribution<double> d(origTxpCount[i] + priorAlphas[i],
                                      1.0 / (beta + effLens(i)));
    muGlobal[i] = d(gen);
  }

  for (auto& eqClass : eqVec) {
    // get total number of reads for an equivalence class
    uint64_t classCount = eqClass.second.count;

    // for each transcript in this class
    const TranscriptGroup& tgroup = eqClass.first;
    const size_t groupSize = tgroup.txps.size();
    if (tgroup.valid) {
      const std::vector<uint32_t>& txps = tgroup.txps;
      const auto& auxs = eqClass.second.combinedWeights;

      double denom = 0.0;
      // If this is a single-transcript group,
      // then it gets the full count --- otherwise,
      // sample!
      if (BOOST_LIKELY(groupSize > 1)) {

        std::vector<uint64_t> txpResamp(groupSize);
        std::vector<double> mu(groupSize);

        // For each transcript in the group
        double muSum = 0.0;
        for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          auto aux = auxs[i];
          // mu[i] = (origTxpCount[tid]+priorAlpha) * aux;
          mu[i] = muGlobal[tid];
          muSum += mu[i];
          denom += (priorAlphas[tid] + origTxpCount[tid]) * aux;
        }

        // calculate prob vector
        for (size_t i = 0; i < groupSize; ++i) {
          probMap[offset + i] = mu[i] / muSum;
          txpResamp[i] = 0.0;
        }

        if (denom > ::minEQClassWeight) {
          // re-sample
          msamp(txpResamp.begin(),       // count array to fill in
                classCount,              // multinomial n
                groupSize,               // multinomial k
                probMap.begin() + offset // where to find multinomial probs
                );

          for (size_t i = 0; i < groupSize; ++i) {
            auto tid = txps.at(i);
            txpCount.at(tid) += txpResamp.at(i);
            // txpCount.at(tid) -= countMap.at(offset + i);
            // countMap.at(offset + i) = txpResamp.at(i);
          }
        } // do nothing when denom less than minEQClassWeight
        else {
          std::cerr << "minEQClassWeight error";
        }
      } // do nothing if group size less than 2
      else {
        auto tid = txps.at(0);
        txpCount.at(tid) += countMap.at(offset);
      }
      offset += groupSize;
    } // valid group
  }   // loop over all eq classes
}



void sampleRound_(
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        std::vector<uint64_t>& countMap,
        std::vector<double>& probMap,
        Eigen::VectorXd& effLens,
        const std::vector<double>& priorAlphas,
        std::vector<int>& txpCount,
        MultinomialSampler& msamp) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.25, 0.75);
    size_t offset{0};
    // Choose a fraction of this class to re-sample

    // The count substracted from each transcript
    std::vector<uint64_t> txpResamp;

    for (auto& eqClass : eqVec) {
        uint64_t classCount = eqClass.second.count;
        double sampleFrac = dis(gen);

        // for each transcript in this class
        const TranscriptGroup& tgroup = eqClass.first;
        const size_t groupSize = tgroup.txps.size();
        if (tgroup.valid) {
            const std::vector<uint32_t>& txps = tgroup.txps;
            const auto& auxs = eqClass.second.combinedWeights;

            double denom = 0.0;
            // If this is a single-transcript group,
            // then it gets the full count --- otherwise,
            // sample!
            if (BOOST_LIKELY(groupSize > 1)) {

                // Subtract some fraction of the current equivalence
                // class' contribution from each transcript.
                uint64_t numResampled{0};
                if (groupSize > txpResamp.size()) {
                    txpResamp.resize(groupSize, 0);
                }

                // For each transcript in the group
                for (size_t i = 0; i < groupSize; ++i) {
                    auto tid = txps[i];
                    auto aux = auxs[i];
                    auto currCount = countMap[offset + i];
                    uint64_t currResamp = std::round(sampleFrac * currCount);
                    numResampled += currResamp;
                    txpResamp[i] = currResamp;
                    txpCount[tid] -= currResamp;
                    countMap[offset + i] -= currResamp;
                    denom += (priorAlphas[tid] + txpCount[tid]) * aux;
                }

                if (denom > ::minEQClassWeight) {
                    // Get the multinomial probabilities
                    double norm = 1.0 / denom;
                    for (size_t i = 0; i < groupSize; ++i) {
                        auto tid = txps[i];
                        auto aux = auxs[i];
                        probMap[offset + i] = norm * ((priorAlphas[tid] +
txpCount[tid]) * aux);
                    }

                    // re-sample
                    msamp(txpResamp.begin(),        // count array to fill in
                            numResampled,		// multinomial n
                            groupSize,		// multinomial k
                            probMap.begin() + offset  // where to find
multinomial probs
                         );

                    for (size_t i = 0; i < groupSize; ++i) {
                        auto tid = txps[i];
                        countMap[offset + i] += txpResamp[i];
                        txpCount[tid] += txpResamp[i];
                    }

                } else { // We didn't sample
                    // add back to txp-count!
                    for (size_t i = 0; i < groupSize; ++i) {
                        auto tid = txps[i];
                        txpCount[tid] += txpResamp[i];
                        countMap[offset + i] += txpResamp[i];
                    }
                }
            }

            offset += groupSize;
        } // valid group
    } // loop over all eq classes

}

// The original sampler!
template <typename ExpT>
bool CollapsedGibbsSampler::sampleMultipleChains(ExpT& readExp,
        SalmonOpts& sopt,
        std::function<bool(const std::vector<double>&)>& writeBootstrap,
        uint32_t numSamples) {

    namespace bfs = boost::filesystem;
    auto& jointLog = sopt.jointLog;
    tbb::task_scheduler_init tbbScheduler(sopt.numThreads);
    std::vector<Transcript>& transcripts = readExp.transcripts();

    // Fill in the effective length vector
    Eigen::VectorXd effLens(transcripts.size());

    std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec =
        readExp.equivalenceClassBuilder().eqVec();

    using VecT = CollapsedGibbsSampler::VecType;

    std::vector<std::vector<int>> allSamples(numSamples,
                                        std::vector<int>(transcripts.size(),0));

    bool perTranscriptPrior = (sopt.useVBOpt) ? sopt.perTranscriptPrior : true;
    double priorValue = (sopt.useVBOpt) ? sopt.vbPrior : 1e-8;
    std::vector<double> priorAlphas = populatePriorAlphasGibbs_(transcripts,
effLens, priorValue, perTranscriptPrior);
    std::vector<double> alphasIn(priorAlphas.size(), 0.0);

    bool useScaledCounts = (!sopt.useQuasi and !sopt.allowOrphans);
    auto numMappedFragments = (useScaledCounts) ? readExp.upperBoundHits() :
readExp.numMappedFragments();
    uint32_t numInternalRounds = sopt.thinningFactor;

    for (size_t i = 0; i < transcripts.size(); ++i) {
        auto& txp = transcripts[i];
        //txp.setMass(priorAlphas[i] + (txp.mass(false) * numMappedFragments));
        alphasIn[i] = txp.mass(false) * numMappedFragments;
        effLens(i) = txp.EffectiveLength;
    }

    tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(numSamples)),
                 [&eqVec, &transcripts, &alphasIn, &priorAlphas, &effLens,
                  &allSamples, &writeBootstrap, useScaledCounts,
numInternalRounds,
                 &jointLog, numMappedFragments]( const BlockedIndexRange& range)
-> void {


                std::random_device rd;
                MultinomialSampler ms(rd);

                size_t countMapSize{0};
                for (size_t i = 0; i < eqVec.size(); ++i) {
                    if (eqVec[i].first.valid) {
                    countMapSize += eqVec[i].first.txps.size();
                    }
                }

                size_t numTranscripts{transcripts.size()};

                // will hold estimated counts
                std::vector<int> alphas(numTranscripts, 0.0);
                std::vector<uint64_t> countMap(countMapSize, 0);
                std::vector<double> probMap(countMapSize, 0.0);

                initCountMap_(eqVec, transcripts, alphasIn, priorAlphas, ms,
countMap, probMap, effLens, allSamples[range.begin()]);

                // For each sample this thread should generate
                bool isFirstSample{true};
                for (auto sampleID : boost::irange(range.begin(), range.end()))
{
                    if (sampleID % 100 == 0) {
                        std::cerr << "gibbs sampling " << sampleID << "\n";
                    }

                    if (!isFirstSample) {
                        // the counts start at what they were last round.
                        allSamples[sampleID] = allSamples[sampleID-1];
                    }

                    // Thin the chain by a factor of (numInternalRounds)
                    for (size_t i = 0; i < numInternalRounds; ++i){
                      sampleRoundNonCollapsed_(eqVec, countMap, probMap,
effLens, priorAlphas,
                                allSamples[sampleID], ms);
                    }

                    // If we're scaling the counts, do it here.
                    if (useScaledCounts) {
                        double numMappedFrags =
static_cast<double>(numMappedFragments);
                        double alphaSum = 0.0;
                        for (auto c : allSamples[sampleID]) { alphaSum +=
static_cast<double>(c); }
                        if (alphaSum > ::minWeight) {
                            double scaleFrac = 1.0 / alphaSum;
                            // scaleFrac converts alpha to nucleotide fraction,
                            // and multiplying by numMappedFrags scales by the
total
                            // number of mapped fragments to provide an
estimated count.
                            for (size_t tn = 0; tn < numTranscripts; ++tn) {
                                alphas[tn] = static_cast<int>(
                                        std::round(
                                            numMappedFrags *
                                            (static_cast<double>(allSamples[sampleID][tn])
* scaleFrac)));
                            }
                        } else { // This shouldn't happen!
                            jointLog->error("Gibbs sampler had insufficient
number of fragments!"
                                    "Something is probably wrong; please check
that you "
                                    "have run salmon correctly and report this
to GitHub.");
                        }
                    } else { // otherwise, just copy over from the sampled
counts
                        for (size_t tn = 0; tn < numTranscripts; ++tn) {
                            alphas[tn] =
static_cast<int>(allSamples[sampleID][tn]);
                        }
                    }

                    writeBootstrap(alphas);
                    //bootstrapWriter->writeBootstrap(alphas);
                    isFirstSample = false;
                }
    });
    return true;
}
*/

using SCExpT = ReadExperiment<EquivalenceClassBuilder<SCTGValue>>;
using BulkExpT = ReadExperiment<EquivalenceClassBuilder<TGValue>>;
template <typename FragT>
using BulkAlignLibT = AlignmentLibrary<FragT, EquivalenceClassBuilder<TGValue>>;

template bool CollapsedGibbsSampler::sample<BulkExpT>(
    BulkExpT& readExp, SalmonOpts& sopt,
    std::function<bool(const std::vector<double>&)>& writeBootstrap,
    uint32_t maxIter);

template bool CollapsedGibbsSampler::sample<BulkAlignLibT<UnpairedRead>>(
    BulkAlignLibT<UnpairedRead>& readExp, SalmonOpts& sopt,
    std::function<bool(const std::vector<double>&)>& writeBootstrap,
    uint32_t maxIter);

template bool CollapsedGibbsSampler::sample<BulkAlignLibT<ReadPair>>(
    BulkAlignLibT<ReadPair>& readExp, SalmonOpts& sopt,
    std::function<bool(const std::vector<double>&)>& writeBootstrap,
    uint32_t maxIter);
/*
template
bool CollapsedGibbsSampler::sampleMultipleChains<ReadExperiment>(ReadExperiment&
readExp,
                                                   SalmonOpts& sopt,
                                                   std::function<bool(const
std::vector<double>&)>& writeBootstrap,
                                                   uint32_t maxIter);

template
bool
CollapsedGibbsSampler::sampleMultipleChains<AlignmentLibrary<UnpairedRead>>(
                                                                   AlignmentLibrary<UnpairedRead>&
readExp,
                                                                   SalmonOpts&
sopt,
                                                                   std::function<bool(const
std::vector<double>&)>& writeBootstrap,
                                                                   uint32_t
maxIter);


template
bool CollapsedGibbsSampler::sampleMultipleChains<AlignmentLibrary<ReadPair>>(
                                                               AlignmentLibrary<ReadPair>&
readExp,
                                                               SalmonOpts& sopt,
                                                               std::function<bool(const
std::vector<double>&)>& writeBootstrap,
                                                               uint32_t
maxIter);
*/

/*
    // Deprecated Gibbs output code
    auto numTranscripts = transcripts.size();
    std::vector<DistStats> ds(numTranscripts);

    // get posterior means
    tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(numTranscripts)),
                [&allSamples, &transcripts, &ds, numMappedFragments,
                 numSamples]( const BlockedIndexRange& range) -> void {

                // For each sample this thread should generate
                for (auto tid : boost::irange(range.begin(), range.end())) {
                    double meanNumReads = {0.0};
                    for (size_t i = 0; i < numSamples; ++i) {
                      auto val = allSamples[i][tid];
                      if (val < ds[tid].minVal) { ds[tid].minVal = val; }
                      if (val > ds[tid].maxVal) { ds[tid].maxVal = val; }
                      meanNumReads += (1.0 / numSamples) * val;
                    }
                    ds[tid].meanVal = meanNumReads;
                    transcripts[tid].setMass(ds[tid].meanVal);
                }
    });

    bfs::path gibbsSampleFile = sopt.outputDirectory / "samples.txt";
    sopt.jointLog->info("Writing posterior samples to {}",
   gibbsSampleFile.string());

    std::ofstream statStream(gibbsSampleFile.string());
    statStream << "# txpName\tsample_1\tsample_2\t...\tsample_n\n";

    for (size_t i = 0; i < numTranscripts; ++i) {
        statStream << transcripts[i].RefName;
        for (size_t s = 0; s < allSamples.size(); ++s) {
            statStream << '\t' << allSamples[s][i];
        }
        statStream << '\n';
    }
    statStream.close();
    sopt.jointLog->info("done writing posterior samples");

    double cutoff = priorAlpha + 1e-8;
    // Truncate tiny expression values
    double txpSumTrunc = 0.0;
    for (size_t i = 0; i < transcripts.size(); ++i) {
    // maybe use the first decile instead of the mean for the cutoff;
    // this could let too much through
        if (transcripts[i].mass(false) <= cutoff) { transcripts[i].setMass(0.0);
   }
        txpSumTrunc += transcripts[i].mass(false);
    }

    for (size_t i = 0; i < transcripts.size(); ++i) {
        // Set the mass to the normalized (after truncation)
        // relative abundance
        transcripts[i].setMass(transcripts[i].mass(false) / txpSumTrunc);
    }

*/
