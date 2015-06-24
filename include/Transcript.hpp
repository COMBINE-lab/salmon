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
        RefName(name), RefLength(len), id(idIn), Sequence(nullptr),
        logPerBasePrior_(std::log(alpha)),
        priorMass_(std::log(alpha*len)),
        mass_(salmon::math::LOG_0), sharedCount_(0.0),
        avgMassBias_(salmon::math::LOG_0) {
            uniqueCount_.store(0);
            lastUpdate_.store(0);
            lastTimestepUpdated_.store(0);
            cachedEffectiveLength_.store(std::log(static_cast<double>(RefLength)));
        }

    Transcript(Transcript&& other) {
        id = other.id;
        //std::swap(RefName, other.RefName);
        RefName = std::move(other.RefName);
        RefLength = other.RefLength;
        Sequence = other.Sequence;
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
    }

    Transcript& operator=(Transcript&& other) {
        id = other.id;
        //std::swap(RefName, other.RefName);
        RefName = std::move(other.RefName);
        RefLength = other.RefLength;
        Sequence = other.Sequence;
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
                return Sequence[byte] & 0x0F;
            } else {
                return ((Sequence[byte] & 0xF0) >> 4) & 0x0F;
            }
            break;
        case strand::reverse:
            if (nibble) {
                return encodedRevComp[Sequence[byte] & 0x0F];
            } else {
                return encodedRevComp[((Sequence[byte] & 0xF0) >> 4) & 0x0F];
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
    double updateEffectiveLength(const FragmentLengthDistribution& fragLengthDist) {

        double effectiveLength = salmon::math::LOG_0;
        double refLen = static_cast<double>(RefLength);
        double logLength = std::log(refLen);

        if (logLength < fragLengthDist.mean()) {
            effectiveLength = logLength;
        } else {
            uint32_t mval = fragLengthDist.maxVal();
            for (size_t l = fragLengthDist.minVal(); l <= std::min(RefLength, mval); ++l) {
                effectiveLength = salmon::math::logAdd(
                        effectiveLength,
                        fragLengthDist.pmf(l) + std::log(refLen - l + 1));
            }
        }

        return effectiveLength;
    }

    double getCachedEffectiveLength() {
        return cachedEffectiveLength_.load();
    }

    double getEffectiveLength(const FragmentLengthDistribution& fragLengthDist,
                              size_t currObs,
                              size_t burnInObs) {
        if (lastUpdate_ == 0 or
            (currObs - lastUpdate_ >= 250000) or
            (lastUpdate_ < burnInObs and currObs > burnInObs)) {
            // compute new number
            double cel = updateEffectiveLength(fragLengthDist);
            cachedEffectiveLength_.store(cel);
            lastUpdate_.store(currObs);
            //priorMass_ = cel + logPerBasePrior_;
            return cachedEffectiveLength_.load();
        } else {
            // return cached number
            return cachedEffectiveLength_.load();
        }
    }

    double perBasePrior() { return std::exp(logPerBasePrior_); }
    inline size_t lastTimestepUpdated() { return lastTimestepUpdated_.load(); }

    void lengthClassIndex(uint32_t ind) { lengthClassIndex_ = ind; }
    uint32_t lengthClassIndex() { return lengthClassIndex_; }

    std::string RefName;
    uint32_t RefLength;
    uint32_t id;

    double uniqueCounts{0.0};
    double totalCounts{0.0};
    double projectedCounts{0.0};
    double sharedCounts{0.0};

    uint8_t* Sequence;

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
};

#endif //TRANSCRIPT
