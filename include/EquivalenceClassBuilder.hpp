#ifndef EQUIVALENCE_CLASS_BUILDER_HPP
#define EQUIVALENCE_CLASS_BUILDER_HPP

#include <unordered_map>
#include <vector>
#include <thread>
#include <memory>

// Logger includes
#include "spdlog/spdlog.h"

#include "cuckoohash_map.hh"
#include "concurrentqueue.h"
#include "TranscriptGroup.hpp"


struct TGValue {
    TGValue(const TGValue& o) {
        weights = o.weights;
        count.store(o.count.load());
    }

    TGValue(std::vector<double>& weightIn, uint64_t countIn) :
        weights(weightIn) { count.store(countIn); }

    // const is a lie
    void normalizeAux() const {
        double sumOfAux = salmon::math::LOG_0;
        for (size_t i = 0; i < weights.size(); ++i) {
            sumOfAux = salmon::math::logAdd(sumOfAux, weights[i]);
        }
        for (size_t i = 0; i < weights.size(); ++i) {
            weights[i] = std::exp(weights[i] - sumOfAux);
        }
    }

    // forget synchronizing this for the time being
    mutable std::vector<double> weights;
    std::atomic<uint64_t> count{0};
};

class EquivalenceClassBuilder {
    public:
        EquivalenceClassBuilder(std::shared_ptr<spdlog::logger> loggerIn) :
		logger_(loggerIn) {
            countMap_.reserve(1000000);
        }

        ~EquivalenceClassBuilder() {}

        void start() { active_ = true; }

        bool finish() {
            active_ = false;
            for (auto kv = countMap_.begin(); !kv.is_end(); ++kv) {
                kv->second.normalizeAux();
            }
            countVec_ = countMap_.snapshot_table();
    	    logger_->info("Computed {} weighted equivalence classes "
			  "for further processing", countVec_.size());
            return true;
        }

        inline void addGroup(TranscriptGroup&& g,
                             std::vector<double>& weights) {

            auto upfn = [&weights](TGValue& g) -> TGValue& {

                g.count++;

                for (size_t i = 0; i < g.weights.size(); ++i) {
                    g.weights[i] =
                        salmon::math::logAdd(g.weights[i], weights[i]);
                }

                return g;
            };
            TGValue v(weights, 1);
            countMap_.upsert(g, upfn, v);
        }

	    cuckoohash_map<TranscriptGroup, TGValue, TranscriptGroupHasher>& eqMap() {
            return countMap_;
        }

        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec() {
            return countVec_;
        }

    private:
        std::atomic<bool> active_;
	    cuckoohash_map<TranscriptGroup, TGValue, TranscriptGroupHasher> countMap_;
        std::vector<std::pair<const TranscriptGroup, TGValue>> countVec_;
    	std::shared_ptr<spdlog::logger> logger_;
};

#endif // EQUIVALENCE_CLASS_BUILDER_HPP

