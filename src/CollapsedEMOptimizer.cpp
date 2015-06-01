#include <vector>
#include <unordered_map>
#include <atomic>

#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/partitioner.h"

//#include "fastapprox.h"
#include <boost/math/special_functions/digamma.hpp>

// C++ string formatting library
#include "format.h"

#include "cuckoohash_map.hh"
#include "Eigen/Dense"

#include "CollapsedEMOptimizer.hpp"
#include "Transcript.hpp"
#include "TranscriptGroup.hpp"
#include "SalmonMath.hpp"
#include "AlignmentLibrary.hpp"
#include "ReadPair.hpp"
#include "UnpairedRead.hpp"
#include "ReadExperiment.hpp"


using BlockedIndexRange =  tbb::blocked_range<size_t>;

// intelligently chosen value adopted from
// https://github.com/pachterlab/kallisto/blob/master/src/EMAlgorithm.h#L18
constexpr double minEQClassWeight = std::numeric_limits<double>::denorm_min();
constexpr double minWeight = std::numeric_limits<double>::denorm_min();

double normalize(std::vector<tbb::atomic<double>>& vec) {
    double sum{0.0};
    for (auto& v : vec) {
        sum += v;
    }

    // too small!
    if (sum < minWeight) {
        return sum;
    }

    double invSum  = 1.0 / sum;
    for (auto& v : vec) {
        v.store(v.load() * invSum);
    }

    return sum;
}

void incLoop(tbb::atomic<double>& val, double inc) {
        double oldMass = val.load();
        double returnedMass = oldMass;
        double newMass{oldMass + inc};
        do {
            oldMass = returnedMass;
            newMass = oldMass + inc;
            returnedMass = val.compare_and_swap(newMass, oldMass);
        } while (returnedMass != oldMass);
}

void CollapsedEMOptimizer::EMUpdate_(
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        std::vector<Transcript>& transcripts,
        Eigen::VectorXd& effLens,
        const CollapsedEMOptimizer::VecType& alphaIn,
        CollapsedEMOptimizer::VecType& alphaOut) {

    assert(alphaIn.size() == alphaOut.size());
    /*
    // serial implementation
    // for each equivalence class
    for (auto& kv : eqVec) {

        // for each transcript in this class
        const TranscriptGroup& tgroup = kv.first;
        if (tgroup.valid) {
            uint64_t count = kv.second.count;
            const std::vector<uint32_t>& txps = tgroup.txps;
            const std::vector<double>& auxs = kv.second.weights;

            double denom = 0.0;
            for (size_t i = 0; i < txps.size(); ++i) {
                auto tid = txps[i];
                auto aux = auxs[i];
                // double el = effLens(tid);
                // if (el <= 0) { el = 1.0; }
                double v = alphaIn(tid) * aux;
                denom += v;
            }

            if (denom <= minEQClassWeight) {
                // tgroup.setValid(false);
            } else {
                double invDenom = 1.0 / denom;
                for (size_t i = 0; i < txps.size(); ++i) {
                    auto tid = txps[i];
                    auto aux = auxs[i];
                    double v = alphaIn(tid) * aux;
                    if (!std::isnan(v)) {
                        alphaOut(tid) += count * v * invDenom;
                    }
                }
            }
        }
    }
    */


    tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(eqVec.size())),
            [&eqVec, &alphaIn, &alphaOut](const BlockedIndexRange& range) -> void {
            for (auto eqID : boost::irange(range.begin(), range.end())) {
            auto& kv = eqVec[eqID];

            uint64_t count = kv.second.count;
            // for each transcript in this class
            const TranscriptGroup& tgroup = kv.first;
            if (tgroup.valid) {
                const std::vector<uint32_t>& txps = tgroup.txps;
                const std::vector<double>& auxs = kv.second.weights;

                double denom = 0.0;
                for (size_t i = 0; i < txps.size(); ++i) {
                auto tid = txps[i];
                auto aux = auxs[i];
                //double el = effLens(tid);
                //if (el <= 0) { el = 1.0; }
                double v = alphaIn[tid] * aux;
                denom += v;
                }

                if (denom <= ::minEQClassWeight) {
                    // tgroup.setValid(false);
                } else {

                    double invDenom = 1.0 / denom;
                    for (size_t i = 0; i < txps.size(); ++i) {
                        auto tid = txps[i];
                        auto aux = auxs[i];
                        double v = alphaIn[tid] * aux;
                        if (!std::isnan(v)) {
                            incLoop(alphaOut[tid], count * v * invDenom);
                        }
                    }
                }
            }
    }
    });

    /*
    // normalize
    //double aos = normalize(alphaOut);
    double aos = alphaOut.sum();
    alphaOut /= aos;
    */
}

void VBEMUpdate_(
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        std::vector<Transcript>& transcripts,
        Eigen::VectorXd& effLens,
        double totLen,
        const CollapsedEMOptimizer::VecType& alphaIn,
        CollapsedEMOptimizer::VecType& alphaOut,
	    CollapsedEMOptimizer::VecType& expTheta) {

    assert(alphaIn.size() == alphaOut.size());

    double alphaSum = {0.0};
    for (auto& e : alphaIn) { alphaSum += e; }
    double logNorm = boost::math::digamma(alphaSum);

    tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(transcripts.size())),
            [logNorm, totLen, &effLens, &alphaIn,
             &alphaOut, &expTheta]( const BlockedIndexRange& range) -> void {

             double prior = 0.01;
             double priorNorm = prior * totLen;

             for (auto i : boost::irange(range.begin(), range.end())) {
                if (alphaIn[i] > ::minWeight) {
                    expTheta[i] = std::exp(boost::math::digamma(alphaIn[i].load()) - logNorm);
                } else {
                    expTheta[i] = 0.0;
                }
                alphaOut[i] = (prior * effLens(i)) / priorNorm;
            }
        });

    tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(eqVec.size())),
            [&eqVec, &alphaIn,
             &alphaOut, &expTheta]( const BlockedIndexRange& range) -> void {
            for (auto eqID : boost::irange(range.begin(), range.end())) {
            auto& kv = eqVec[eqID];

            uint64_t count = kv.second.count;
            // for each transcript in this class
            const TranscriptGroup& tgroup = kv.first;
            if (tgroup.valid) {
                const std::vector<uint32_t>& txps = tgroup.txps;
                const std::vector<double>& auxs = kv.second.weights;


                double denom = 0.0;
                for (size_t i = 0; i < txps.size(); ++i) {
                    auto tid = txps[i];
                    auto aux = auxs[i];
                    //double el = effLens(tid);
                    //if (el <= 0) { el = 1.0; }
                    if (expTheta[tid] > 0.0) {
                        double v = expTheta[tid] * aux;
                        denom += v;
                   }
                }
                if (denom <= ::minEQClassWeight) {
                    // tgroup.setValid(false);
                } else {
                    double invDenom = 1.0 / denom;
                    for (size_t i = 0; i < txps.size(); ++i) {
                        auto tid = txps[i];
                        auto aux = auxs[i];
                	    if (expTheta[tid] > 0.0) {
                          double v = expTheta[tid] * aux;
                          incLoop(alphaOut[tid], count * v * invDenom);
	                    }
                    }
                }
            }
        }});
    // normalize
    //double aos = alphaOut.sum();
    //alphaOut /= aos;

}

size_t markDegenerateClasses(
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        CollapsedEMOptimizer::VecType& alphaIn,
        std::shared_ptr<spdlog::logger> jointLog,
        bool verbose=false) {

    size_t numDropped{0};
    size_t idx{0};
    for (auto& kv : eqVec) {
        uint64_t count = kv.second.count;
        // for each transcript in this class
        const TranscriptGroup& tgroup = kv.first;
        const std::vector<uint32_t>& txps = tgroup.txps;
        const std::vector<double>& auxs = kv.second.weights;

        double denom = 0.0;
        for (size_t i = 0; i < txps.size(); ++i) {
            auto tid = txps[i];
            auto aux = auxs[i];
            double v = alphaIn[tid] * aux;
            if (!std::isnan(v)) {
                denom += v;
            } else {
                std::cerr << "val is NAN; alpha( "
                          << tid << " ) = " << alphaIn[tid]
                          << ", aux = " << aux << "\n";
            }
        }
        if (denom <= minEQClassWeight) {
            fmt::MemoryWriter errstream;

            errstream << "\nDropping weighted eq class\n";
            errstream << "============================\n";

            errstream << "denom = 0, count = " << count << "\n";
            errstream << "class = { ";
            for (auto e : txps) {
                errstream << e << " ";
            }
            errstream << "}\n";
            errstream << "alphas = { ";
            for (auto e : txps) {
                errstream << alphaIn[e] << " ";
            }
            errstream << "}\n";
            errstream << "weights = { ";
            for (auto e : auxs) {
                errstream << e << " ";
            }
            errstream << "}\n";
            errstream << "============================\n\n";

            jointLog->info(errstream.str());
            ++numDropped;
            kv.first.setValid(false);
        }
    }
    return numDropped;
}


CollapsedEMOptimizer::CollapsedEMOptimizer() {}

template <typename ExpT>
bool CollapsedEMOptimizer::optimize(ExpT& readExp,
        SalmonOpts& sopt,
        double relDiffTolerance,
        uint32_t maxIter) {

    tbb::task_scheduler_init tbbScheduler(sopt.numThreads);
    std::vector<Transcript>& transcripts = readExp.transcripts();

    using VecT = CollapsedEMOptimizer::VecType;
    // With Eigen
    /*
    Eigen::VectorXd alphas(transcripts.size());
    Eigen::VectorXd alphasPrime(transcripts.size());
    alphasPrime.setZero();

    Eigen::VectorXd effLens(transcripts.size());
    Eigen::VectorXd expTheta(transcripts.size());
    */
    // With atomics
    VecType alphas(transcripts.size(), 0.0);
    VecType alphasPrime(transcripts.size(), 0.0);
    VecType expTheta(transcripts.size());
    Eigen::VectorXd effLens(transcripts.size());

    std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec =
        readExp.equivalenceClassBuilder().eqVec();

    bool useVBEM{sopt.useVBOpt};
    auto jointLog = sopt.jointLog;

    double logTotalMass = salmon::math::LOG_0;
    double totalLen{0.0};

    for (size_t i = 0; i < transcripts.size(); ++i) {
        double m = transcripts[i].mass(false);
        if (std::isnan(m)) {
            std::cerr << "FOUND NAN for txp " << i << "\n";
        }
        alphas[i] = (m == salmon::math::LOG_0) ? 0.0 : m;
        effLens(i) = std::exp(transcripts[i].getCachedEffectiveLength());
        totalLen += effLens(i);
    }

    auto numRemoved = markDegenerateClasses(eqVec, alphas, sopt.jointLog);
    sopt.jointLog->info("Marked {} weighted equivalence classes as degenerate",
            numRemoved);

    size_t itNum{0};
    double minAlpha = 1e-8;

    bool converged{false};
    double maxRelDiff = -std::numeric_limits<double>::max();
    while (itNum < maxIter and !converged) {

        if (useVBEM) {
            VBEMUpdate_(eqVec, transcripts, effLens, totalLen, alphas, alphasPrime, expTheta);
        } else {
            EMUpdate_(eqVec, transcripts, effLens, alphas, alphasPrime);
        }

        converged = true;
        maxRelDiff = -std::numeric_limits<double>::max();
        for (size_t i = 0; i < transcripts.size(); ++i) {
            if (alphas[i] > minAlpha) {
                double relDiff = std::abs(alphas[i] - alphasPrime[i]) / alphas[i];
                maxRelDiff = (relDiff > maxRelDiff) ? relDiff : maxRelDiff;
                if (relDiff > relDiffTolerance) {
                    converged = false;
                }
            }
            alphas[i] = alphasPrime[i];
            alphasPrime[i] = 0.0;
        }

        if (itNum % 100 == 0) {
            jointLog->info("iteration = {} | max rel diff. = {}",
                            itNum, maxRelDiff);
        }

        ++itNum;
    }

    jointLog->info("iteration = {} | max rel diff. = {}",
                    itNum, maxRelDiff);

    // Truncate tiny expression values
    double alphaSum = 0.0;
    for (size_t i = 0; i < alphas.size(); ++i) {
        if (alphas[i] <= minAlpha) { alphas[i] = 0.0; }
        alphaSum += alphas[i];
    }

    if (alphaSum < minWeight) {
        jointLog->error("Total alpha weight was too small! "
                        "Make sure you ran salmon correclty.");
        return false;
    }

    // Set the mass of each transcript using the
    // computed alphas.
    for (size_t i = 0; i < transcripts.size(); ++i) {
        // Set the mass to the normalized (after truncation)
        // relative abundance
        transcripts[i].setMass(alphas[i] / alphaSum);
    }

    return true;
}

template
bool CollapsedEMOptimizer::optimize<ReadExperiment>(ReadExperiment& readExp,
        SalmonOpts& sopt,
        double relDiffTolerance,
        uint32_t maxIter);

template
bool CollapsedEMOptimizer::optimize<AlignmentLibrary<UnpairedRead>>(
        AlignmentLibrary<UnpairedRead>& readExp,
        SalmonOpts& sopt,
        double relDiffTolerance,
        uint32_t maxIter);


template
bool CollapsedEMOptimizer::optimize<AlignmentLibrary<ReadPair>>(
        AlignmentLibrary<ReadPair>& readExp,
        SalmonOpts& sopt,
        double relDiffTolerance,
        uint32_t maxIter);
