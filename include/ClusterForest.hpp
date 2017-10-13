#ifndef __CLUSTER_FOREST_HPP__
#define __CLUSTER_FOREST_HPP__

#include <boost/pending/disjoint_sets.hpp>

#include "SalmonSpinLock.hpp"
#include "Transcript.hpp"
#include "TranscriptCluster.hpp"

#include <mutex>
#include <unordered_set>
#include <vector>

/** A forest of transcript clusters */
class ClusterForest {
public:
  ClusterForest(size_t numTranscripts, std::vector<Transcript>& refs)
      : rank_(std::vector<size_t>(numTranscripts, 0)),
        parent_(std::vector<size_t>(numTranscripts, 0)),
        disjointSets_(&rank_[0], &parent_[0]),
        clusters_(std::vector<TranscriptCluster>(numTranscripts)) {
    // Initially make a unique set for each transcript
    for (size_t tnum = 0; tnum < numTranscripts; ++tnum) {
      disjointSets_.make_set(tnum);
      clusters_[tnum].members_.push_front(tnum);
      clusters_[tnum].addMass(refs[tnum].mass());
    }
  }

  template <typename FragT>
  void mergeClusters(typename std::vector<FragT>::iterator start,
                     typename std::vector<FragT>::iterator finish) {
// Use a lock_guard to ensure this is a locked (and exception-safe) operation
#if defined __APPLE__
    spin_lock::scoped_lock sl(clusterMutex_);
#else
    std::lock_guard<std::mutex> lock(clusterMutex_);
#endif
    size_t firstCluster, otherCluster;
    auto firstTranscriptID = start->transcriptID();
    ++start;

    for (auto it = start; it != finish; ++it) {
      firstCluster = disjointSets_.find_set(firstTranscriptID);
      otherCluster = disjointSets_.find_set(it->transcriptID());
      if (otherCluster != firstCluster) {
        disjointSets_.link(firstCluster, otherCluster);
        auto parentClust = disjointSets_.find_set(it->transcriptID());
        auto childClust =
            (parentClust == firstCluster) ? otherCluster : firstCluster;
        if (parentClust == firstCluster or parentClust == otherCluster) {
          clusters_[parentClust].merge(clusters_[childClust]);
          clusters_[childClust].deactivate();
        } else {
          std::cerr << "DANGER\n";
        }
      }
    }
  }

  template <typename FragT>
  void mergeClusters(typename std::vector<FragT*>::iterator start,
                     typename std::vector<FragT*>::iterator finish) {
// Use a lock_guard to ensure this is a locked (and exception-safe) operation
#if defined __APPLE__
    spin_lock::scoped_lock sl(clusterMutex_);
#else
    std::lock_guard<std::mutex> lock(clusterMutex_);
#endif

    auto firstTranscriptID = (*start)->transcriptID();
    decltype(firstTranscriptID) firstCluster, otherCluster;
    ++start;

    for (auto it = start; it != finish; ++it) {
      firstCluster = disjointSets_.find_set(firstTranscriptID);
      otherCluster = disjointSets_.find_set((*it)->transcriptID());
      if (otherCluster != firstCluster) {
        disjointSets_.link(firstCluster, otherCluster);
        auto parentClust = disjointSets_.find_set((*it)->transcriptID());
        auto childClust =
            (parentClust == firstCluster) ? otherCluster : firstCluster;
        if (parentClust == firstCluster or parentClust == otherCluster) {
          clusters_[parentClust].merge(clusters_[childClust]);
          clusters_[childClust].deactivate();
        } else {
          std::cerr << "DANGER\n";
        }
      }
    }
  }

  /*
  void mergeClusters(AlignmentBatch<ReadPair>::iterator start,
  AlignmentBatch<ReadPair>::iterator finish) {
      // Use a lock_guard to ensure this is a locked (and exception-safe)
  operation std::lock_guard<std::mutex> lock(clusterMutex_); size_t
  firstCluster, otherCluster; auto firstTranscriptID = start->read1->core.tid;
      ++start;

      for (auto it = start; it != finish; ++it) {
          firstCluster = disjointSets_.find_set(firstTranscriptID);
          otherCluster = disjointSets_.find_set(it->read1->core.tid);
          if (otherCluster != firstCluster)  {
              disjointSets_.link(firstCluster, otherCluster);
              auto parentClust = disjointSets_.find_set(it->read1->core.tid);
              auto childClust = (parentClust == firstCluster)  ? otherCluster :
  firstCluster; if (parentClust == firstCluster or parentClust == otherCluster)
  { clusters_[parentClust].merge(clusters_[childClust]);
                  clusters_[childClust].deactivate();
              } else { std::cerr << "DANGER\n"; }
          }
      }
  }
*/
  void updateCluster(size_t memberTranscript, size_t newCount,
                     double logNewMass, bool updateCount) {
// Use a lock_guard to ensure this is a locked (and exception-safe) operation
#if defined __APPLE__
    spin_lock::scoped_lock sl(clusterMutex_);
#else
    std::lock_guard<std::mutex> lock(clusterMutex_);
#endif
    auto clusterID = disjointSets_.find_set(memberTranscript);
    auto& cluster = clusters_[clusterID];
    if (updateCount) {
      cluster.incrementCount(newCount);
    }
    cluster.addMass(logNewMass);
  }

  std::vector<TranscriptCluster*> getClusters() {
    std::vector<TranscriptCluster*> clusters;
    std::unordered_set<size_t> observedReps;
    for (size_t i = 0; i < clusters_.size(); ++i) {
      auto rep = disjointSets_.find_set(i);
      if (observedReps.find(rep) == observedReps.end()) {
        if (!clusters_[rep].isActive()) {
          std::cerr << "returning a non-active cluster!\n";
          std::exit(1);
        }
        clusters.push_back(&clusters_[rep]);
        observedReps.insert(rep);
      }
    }
    return clusters;
  }

private:
  std::vector<size_t> rank_;
  std::vector<size_t> parent_;
  boost::disjoint_sets<size_t*, size_t*> disjointSets_;
  std::vector<TranscriptCluster> clusters_;
#if defined __APPLE__
  spin_lock clusterMutex_;
#else
  std::mutex clusterMutex_;
#endif
};

#endif // __CLUSTER_FOREST_HPP__
