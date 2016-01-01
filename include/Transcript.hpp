#ifndef TRANSCRIPT
#define TRANSCRIPT

#include <atomic>
#include <cmath>
#include <limits>
#include "SalmonStringUtils.hpp"
#include "SalmonUtils.hpp"
#include "SalmonMath.hpp"
#include "SequenceBiasModel.hpp"
#include "FragmentLengthDistribution.hpp"
#include "tbb/atomic.h"

class Transcript {
public:
    Transcript(size_t idIn, const char* name, uint32_t len, double alpha = 0.05) :
        RefName(name), RefLength(len), EffectiveLength(-1.0), id(idIn), SAMSequence(nullptr), Sequence(nullptr),
        logPerBasePrior_(std::log(alpha)),
        priorMass_(std::log(alpha*len)),
        mass_(salmon::math::LOG_0), sharedCount_(0.0),
        avgMassBias_(salmon::math::LOG_0),
        active_(false),
	freeSeqOnDestruct(false){
            uniqueCount_.store(0);
            lastUpdate_.store(0);
            lastTimestepUpdated_.store(0);
            cachedEffectiveLength_.store(std::log(static_cast<double>(RefLength)));
        }

    ~Transcript() {
      // Free the sequence if it belongs to us
      if (freeSeqOnDestruct) { delete [] Sequence; }
      // Free the SAMSequence if it exists
      if (SAMSequence) { delete [] SAMSequence; }
    }

    Transcript(Transcript&& other) {
        id = other.id;
        //std::swap(RefName, other.RefName);
        RefName = std::move(other.RefName);
        RefLength = other.RefLength;
        EffectiveLength = other.EffectiveLength;
        SAMSequence = other.SAMSequence;
        // If this is an owned-resource, then move it
        if (other.SAMSequence) {
            other.SAMSequence = nullptr;
        }

        Sequence = other.Sequence;
        // If this is an owned-resource, then move it
        if (other.freeSeqOnDestruct) {
            freeSeqOnDestruct = true;
            other.freeSeqOnDestruct = false;
        }

        uniqueCount_.store(other.uniqueCount_);
        totalCount_.store(other.totalCount_.load());
        lastTimestepUpdated_.store(other.lastTimestepUpdated_.load());
        sharedCount_.store(other.sharedCount_.load());
        mass_.store(other.mass_.load());
        lastUpdate_.store(other.lastUpdate_.load());
        cachedEffectiveLength_.store(other.cachedEffectiveLength_.load());
        lengthClassIndex_ = other.lengthClassIndex_;
        logPerBasePrior_ = other.logPerBasePrior_;
        priorMass_ = other.priorMass_;
        avgMassBias_.store(other.avgMassBias_.load());
        hasAnchorFragment_.store(other.hasAnchorFragment_.load());
        active_ = other.active_;
    }

    Transcript& operator=(Transcript&& other) {
        id = other.id;
        //std::swap(RefName, other.RefName);
        RefName = std::move(other.RefName);
        RefLength = other.RefLength;
        EffectiveLength = other.EffectiveLength;
        SAMSequence = other.SAMSequence;
        // If this is an owned-resource, then move it
        if (other.SAMSequence) {
            other.SAMSequence = nullptr;
        }

        Sequence = other.Sequence;
        // If this is an owned-resource, then move it
        if (other.freeSeqOnDestruct) {
            freeSeqOnDestruct = true;
            other.freeSeqOnDestruct = false;
        }

        uniqueCount_.store(other.uniqueCount_);
        totalCount_.store(other.totalCount_.load());
        lastTimestepUpdated_.store(other.lastTimestepUpdated_.load());
        sharedCount_.store(other.sharedCount_.load());
        mass_.store(other.mass_.load());
        lastUpdate_.store(other.lastUpdate_.load());
        cachedEffectiveLength_.store(other.cachedEffectiveLength_.load());
        lengthClassIndex_ = other.lengthClassIndex_;
        logPerBasePrior_ = other.logPerBasePrior_;
        priorMass_ = other.priorMass_;
        avgMassBias_.store(other.avgMassBias_.load());
        hasAnchorFragment_.store(other.hasAnchorFragment_.load());
        active_ = other.active_;
        return *this;
    }


    inline double sharedCount() { return sharedCount_.load(); }
    inline size_t uniqueCount() { return uniqueCount_.load(); }
    inline size_t totalCount() { return totalCount_.load(); }

    inline void addUniqueCount(size_t newCount) { uniqueCount_ += newCount; }
    inline void addTotalCount(size_t newCount) { totalCount_ += newCount; }

    inline double uniqueUpdateFraction() {
        double ambigCount = static_cast<double>(totalCount_ - uniqueCount_);
        return uniqueCount_ / ambigCount;
    }

    inline char charBaseAt(size_t idx,
                              salmon::stringtools::strand dir = salmon::stringtools::strand::forward) {
        return salmon::stringtools::samCodeToChar[baseAt(idx, dir)];
    }

    inline uint8_t baseAt(size_t idx,
                          salmon::stringtools::strand dir = salmon::stringtools::strand::forward) {
        using salmon::stringtools::strand;
        using salmon::stringtools::encodedRevComp;
        size_t byte = idx >> 1;
        size_t nibble = idx & 0x1;

        switch(dir) {
        case strand::forward:
            if (nibble) {
                return SAMSequence[byte] & 0x0F;
            } else {
                return ((SAMSequence[byte] & 0xF0) >> 4) & 0x0F;
            }
            break;
        case strand::reverse:
            if (nibble) {
                return encodedRevComp[SAMSequence[byte] & 0x0F];
            } else {
                return encodedRevComp[((SAMSequence[byte] & 0xF0) >> 4) & 0x0F];
            }
            break;
        }

        return std::numeric_limits<uint8_t>::max();
    }

    inline void setSharedCount(double sc) {
        sharedCount_.store(sc);
    }

    inline void addSharedCount(double sc) {
	    salmon::utils::incLoop(sharedCount_, sc);
    }

    inline void setLastTimestepUpdated(uint64_t currentTimestep) {
        uint64_t oldTimestep = lastTimestepUpdated_;
        if (currentTimestep > oldTimestep) {
            lastTimestepUpdated_ = currentTimestep;
        }
    }

    inline void addBias(double bias) {
	salmon::utils::incLoopLog(avgMassBias_, bias);
    }

    inline void addMass(double mass) {
	salmon::utils::incLoopLog(mass_, mass);
    }

    inline void setMass(double mass) {
        mass_.store(mass);
    }

    inline double mass(bool withPrior=true) {
        return (withPrior) ? salmon::math::logAdd(priorMass_, mass_.load()) : mass_.load();
    }

    void setActive() { active_ = true; }
    bool getActive() { return active_; }

    inline double bias() {
        return (totalCount_.load() > 0) ?
                    avgMassBias_ - std::log(totalCount_.load()) :
                    salmon::math::LOG_1;
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
      *  NOTE: Adopted from "est_effective_length" at (https://github.com/adarob/eXpress/blob/master/src/targets.cpp)
      *  originally written by Adam Roberts.
      *
      *
      */
    double computeLogEffectiveLength(
            std::vector<double>& logPMF,
            double logFLDMean,
            size_t minVal,
            size_t maxVal) {

        double effectiveLength = salmon::math::LOG_0;
        double refLen = static_cast<double>(RefLength);
        double logRefLength = std::log(refLen);

        if (logRefLength <= logFLDMean) {
            effectiveLength = logRefLength;
        } else {
            uint32_t mval = maxVal;
            size_t clen = minVal;
            size_t maxLen = std::min(RefLength, mval);
            while (clen <= maxLen) {
                size_t i = clen - minVal;
                effectiveLength = salmon::math::logAdd(
                        effectiveLength,
                        logPMF[i] + std::log(refLen - clen + 1));
                ++clen;
            }
        }
        if (std::exp(effectiveLength) <= 1.0) {
            effectiveLength = salmon::math::LOG_1;
        }

        return effectiveLength;
    }

    /**
     * Return the cached value for the log of the effective length.
     */
    double getCachedLogEffectiveLength() {
        return cachedEffectiveLength_.load();
    }

    void updateEffectiveLength(
            std::vector<double>& logPMF,
            double logFLDMean,
            size_t minVal,
            size_t maxVal) {
        double cel = computeLogEffectiveLength(logPMF, logFLDMean, minVal, maxVal);
        cachedEffectiveLength_.store(cel);
    }

    /**
     * If we should update the effective length, then do it and cache the result.
     * Otherwise, return the cached result.
     */
    /*
    double getLogEffectiveLength(const FragmentLengthDistribution& fragLengthDist,
                                 size_t currObs, size_t burnInObs, bool forceUpdate=false) {
        if (forceUpdate or
            (lastUpdate_ == 0) or
            (currObs - lastUpdate_ >= 250000) or
            (lastUpdate_ < burnInObs and currObs > burnInObs)) {
            // compute new number
            lastUpdate_.store(currObs);
            double cel = computeLogEffectiveLength(fragLengthDist);
            cachedEffectiveLength_.store(cel);
            //priorMass_ = cel + logPerBasePrior_;
            return cachedEffectiveLength_.load();
        } else {
            // return cached number
            return cachedEffectiveLength_.load();
        }
    }
    */

    double perBasePrior() { return std::exp(logPerBasePrior_); }
    inline size_t lastTimestepUpdated() { return lastTimestepUpdated_.load(); }

    void lengthClassIndex(uint32_t ind) { lengthClassIndex_ = ind; }
    uint32_t lengthClassIndex() { return lengthClassIndex_; }

    void setAnchorFragment() {
        hasAnchorFragment_.store(true);
    }

    bool hasAnchorFragment() {
        return hasAnchorFragment_.load();
    }

    std::string RefName;
    uint32_t RefLength;
    double EffectiveLength;
    uint32_t id;

    double uniqueCounts{0.0};
    double totalCounts{0.0};
    double projectedCounts{0.0};
    double sharedCounts{0.0};

    uint8_t* SAMSequence;
    const char* Sequence;
    bool freeSeqOnDestruct;

private:
    std::atomic<size_t> uniqueCount_;
    std::atomic<size_t> totalCount_;
    // The most recent timestep at which this transcript's mass was updated.
    std::atomic<size_t> lastTimestepUpdated_;
    double priorMass_;
    tbb::atomic<double> mass_;
    tbb::atomic<double> sharedCount_;
    tbb::atomic<double> cachedEffectiveLength_;
    tbb::atomic<size_t> lastUpdate_;
    tbb::atomic<double> avgMassBias_;
    uint32_t lengthClassIndex_;
    double logPerBasePrior_;
    // In a paired-end protocol, a transcript has
    // an "anchor" fragment if it has a proper
    // pair of reads mapping to it.
    std::atomic<bool> hasAnchorFragment_{false};
    bool active_;
};

#endif //TRANSCRIPT
