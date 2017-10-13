#ifndef __TRANSCRIPT_CLUSTER_HPP__
#define __TRANSCRIPT_CLUSTER_HPP__

#include <atomic>
#include <iostream>
#include <list>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include "SalmonMath.hpp"

class ClusterForest;

class TranscriptCluster {
  friend class ClusterForest;

public:
  TranscriptCluster()
      : members_(std::list<size_t>()), count_(0), logMass_(salmon::math::LOG_0),
        active_(true) {}
  TranscriptCluster(size_t initialMember)
      : members_(std::list<size_t>(1, initialMember)), count_(0),
        logMass_(salmon::math::LOG_0), active_(true) {}

  // void incrementCount(size_t num) { count_ += num; }
  void incrementCount(double num) { count_ += num; }
  void addMass(double logNewMass) {
    logMass_ = salmon::math::logAdd(logMass_, logNewMass);
  }
  void merge(TranscriptCluster& other) {
    members_.splice(members_.begin(), other.members_);
    logMass_ = salmon::math::logAdd(logMass_, other.logMass_);
    count_ += other.count_;
  }

  std::list<size_t>& members() { return members_; }
  // size_t numHits() { return count_.load(); }
  double numHits() { return count_; }
  bool isActive() { return active_; }
  void deactivate() { active_ = false; }
  double logMass() { return logMass_; }

  // Adapted from https://github.com/adarob/eXpress/blob/master/src/targets.cpp
  void projectToPolytope(std::vector<Transcript>& allTranscripts) {
    using salmon::math::approxEqual;

    // The transcript belonging to this cluster
    double clusterCounts{static_cast<double>(count_)};
    std::vector<Transcript*> transcripts;
    for (auto tid : members_) {
      transcripts.push_back(&allTranscripts[tid]);
    }
    // The cluster size
    size_t clusterSize = transcripts.size();
    size_t round{0};
    boost::dynamic_bitset<> polytopeBound(clusterSize);
    while (true) {
      double unboundCounts{0.0};
      double boundCounts{0.0};
      for (size_t i = 0; i < clusterSize; ++i) {
        Transcript* transcript = transcripts[i];
        if (transcript->projectedCounts > transcript->totalCounts) {
          transcript->projectedCounts = transcript->totalCounts;
          polytopeBound[i] = true;
        } else if (transcript->projectedCounts < transcript->uniqueCounts) {
          transcript->projectedCounts = transcript->uniqueCounts;
          polytopeBound[i] = true;
        }

        if (polytopeBound[i]) {
          boundCounts += transcript->projectedCounts;
        } else {
          unboundCounts += transcript->projectedCounts;
        }
      }

      if (approxEqual(unboundCounts + boundCounts, clusterCounts)) {
        return;
      }

      if (unboundCounts == 0) {
        polytopeBound = boost::dynamic_bitset<>(clusterSize);
        unboundCounts = boundCounts;
        boundCounts = 0;
      }

      double normalizer = (clusterCounts - boundCounts) / unboundCounts;
      for (size_t i = 0; i < clusterSize; ++i) {
        if (!polytopeBound[i]) {
          Transcript* transcript = transcripts[i];
          transcript->projectedCounts *= normalizer;
        }
      }
      ++round;
      if (round > 5000) {
        return;
      }

    } // end while
  }

private:
  std::list<size_t> members_;
  // std::atomic<size_t> count_;
  double count_;
  double logMass_;
  bool active_;
};

#endif // __TRANSCRIPT_CLUSTER_HPP__
