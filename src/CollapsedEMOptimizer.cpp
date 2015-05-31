#include <vector>
#include <unordered_map>

#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/task_scheduler_init.h"
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

void EMUpdate_(
        tbb::task_scheduler_init& tbbInit,
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        std::vector<Transcript>& transcripts,
        Eigen::VectorXd& effLens,
        const Eigen::VectorXd& alphaIn,
        Eigen::VectorXd& alphaOut) {

    assert(alphaIn.size() == alphaOut.size());
    if (!tbbInit.is_active()) {
        tbbInit.initialize(20);
    }

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
                tgroup.setValid(false);
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

    /*
    tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(eqVec.size())),
          [&eqVec, &alphaIn, &alphaOut](const BlockedIndexRange& range) -> void {
            for (auto eqID : boost::irange(range.begin(), range.end())) {
                auto& kv = eqVec[eqID];

                uint64_t count = kv.second.count;
                // for each transcript in this class
                const TranscriptGroup& tgroup = kv.first;
                const std::vector<uint32_t>& txps = tgroup.txps;
                const std::vector<double>& auxs = kv.second.weights;

                double denom = 0.0;
                for (size_t i = 0; i < txps.size(); ++i) {
                    auto tid = txps[i];
                    auto aux = auxs[i];
                    //double el = effLens(tid);
                    //if (el <= 0) { el = 1.0; }
                    double v = alphaIn(tid) * aux;
                    denom += v;
                }

                double invDenom = 1.0 / denom;
                for (size_t i = 0; i < txps.size(); ++i) {
                    auto tid = txps[i];
                    auto aux = auxs[i];
                    double v = alphaIn(tid) * aux;
                    alphaOut(tid) += count * v * invDenom;
                }
            }
      });
    */
    // normalize
    double aos = alphaOut.sum();
    alphaOut /= aos;

}

void VBEMUpdate_(
        //cuckoohash_map<TranscriptGroup, TGValue, TranscriptGroupHasher>& eqMap,
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        std::vector<Transcript>& transcripts,
        Eigen::VectorXd& effLens,
        const Eigen::VectorXd& alphaIn,
	Eigen::VectorXd& expTheta,
        Eigen::VectorXd& alphaOut) {

    assert(alphaIn.size() == alphaOut.size());

    //Eigen::VectorXd logRho(alphaIn.size());
    expTheta.setZero();

    //std::cerr << "alphaIn = " << alphaIn << "\n";
    double alphaSum = alphaIn.sum();
    //std::cerr << "alphaSum = " << alphaSum << "\n";
    double logNorm = boost::math::digamma(alphaSum);
    //std::cerr << "logNorm = " << logNorm << "\n";

    double prior = 0.01;
    double eps = 1e-7;
    for (size_t i = 0; i < transcripts.size(); ++i) {
	    if (alphaIn(i) > prior) {
		    expTheta(i) = std::exp(boost::math::digamma(alphaIn(i)) - logNorm);
	    }
	alphaOut(i) = prior;
    }
    //std::cerr << alphaIn;
    //std::cerr << expTheta;

    // for each equivalence class
    for (auto& kv : eqVec) {
    //for (auto kv = eqMap.cbegin(); !kv.is_end(); ++kv) {

        uint64_t count = kv.second.count;
        // for each transcript in this class
        const TranscriptGroup& tgroup = kv.first;
        const std::vector<uint32_t>& txps = tgroup.txps;
        const std::vector<double>& auxs = kv.second.weights;

        double denom = 0.0;
        for (size_t i = 0; i < txps.size(); ++i) {
            auto tid = txps[i];
            auto aux = auxs[i];
            double el = effLens(tid);
            //if (el <= 0) { el = 1.0; }
	    if (expTheta(tid) > 0.0) {
            	double v = expTheta(tid) * aux;
                denom += v;
	    }
        }
        double invDenom = 1.0 / denom;
        for (size_t i = 0; i < txps.size(); ++i) {
            auto tid = txps[i];
            auto aux = auxs[i];
	    if (expTheta(tid) > 0.0) {
              double v = expTheta(tid) * aux;
              alphaOut(tid) += count * v * invDenom;
	    }
        }
    }
    // normalize
    //double aos = alphaOut.sum();
    //alphaOut /= aos;

}

size_t markDegenerateClasses(
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        Eigen::VectorXd& alphaIn,
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
            double v = alphaIn(tid) * aux;
            if (!std::isnan(v)) {
                denom += v;
            } else {
                std::cerr << "val is NAN; alpha( "
                          << tid << " ) = " << alphaIn(tid)
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
                errstream << alphaIn(e) << " ";
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


    tbb::task_scheduler_init tbbInit(tbb::task_scheduler_init::deferred);
    std::vector<Transcript>& transcripts = readExp.transcripts();
    Eigen::VectorXd alphas(transcripts.size());
    Eigen::VectorXd alphasPrime(transcripts.size());
    alphasPrime.setZero();

    Eigen::VectorXd effLens(transcripts.size());
    Eigen::VectorXd expTheta(transcripts.size());
    std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec =
        readExp.equivalenceClassBuilder().eqVec();

    bool useVBEM{false};
    auto jointLog = sopt.jointLog;

    double logTotalMass = salmon::math::LOG_0;

    for (size_t i = 0; i < transcripts.size(); ++i) {
        double m = transcripts[i].mass(false);
        if (std::isnan(m)) {
            std::cerr << "FOUND NAN for txp " << i << "\n";
        }
        alphas(i) = (m == salmon::math::LOG_0) ? 0.0 : m;
        effLens(i) = std::exp(transcripts[i].getCachedEffectiveLength());
    }

    //std::cerr << "alphas = " << alphas << "\n";

    if (!useVBEM){
        alphas /= alphas.sum();
    }

    auto numRemoved = markDegenerateClasses(eqVec, alphas, sopt.jointLog);
    sopt.jointLog->info("Marked {} weighted equivalence classes as degenerate",
            numRemoved);

    size_t itNum{0};
    double minAlpha = 1e-8;
    //double relDiffTolerance = 0.01;

    bool converged{false};
    double maxRelDiff = -std::numeric_limits<double>::max();
    while (itNum < maxIter and !converged) {

        if (useVBEM) {
            VBEMUpdate_(eqVec, transcripts, effLens, alphas, expTheta, alphasPrime);
        } else {
            EMUpdate_(tbbInit, eqVec, transcripts, effLens, alphas, alphasPrime);
        }

        converged = true;
        maxRelDiff = -std::numeric_limits<double>::max();
        for (size_t i = 0; i < transcripts.size(); ++i) {
            if (alphas(i) > minAlpha) {
                double relDiff = std::abs(alphas(i) - alphasPrime(i)) / alphas(i);
                maxRelDiff = (relDiff > maxRelDiff) ? relDiff : maxRelDiff;
                if (relDiff > relDiffTolerance) {
                    converged = false;
                    continue;
                }
            }
        }

        if (itNum % 100 == 0) {
            jointLog->info("iteration = {} | max rel diff. seen = {}",
                            itNum, maxRelDiff);
        }

        std::swap(alphas, alphasPrime);
        alphasPrime.setZero();
        ++itNum;
    }

    // If we used the VBEM, turn counts
    // into relative nucleotide fractions
    if (useVBEM) {
        alphas /= alphas.sum();
    }

    // Truncate tiny expression values
    for (size_t i = 0; i < alphas.size(); ++i) {
        if (alphas(i) <= minAlpha) { alphas(i) = 0.0; }
    }

    // Re-normalize the expression after truncation
    alphas /= alphas.sum();

    // Set the mass of each transcript using the
    // computed alphas.
    for (size_t i = 0; i < transcripts.size(); ++i) {
        transcripts[i].setMass(alphas(i));
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
