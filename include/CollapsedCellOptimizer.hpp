#ifndef COLLAPSED_CELL_OPTIMIZER_HPP
#define COLLAPSED_CELL_OPTIMIZER_HPP

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <thread>

#include "tbb/atomic.h"
#include "tbb/task_scheduler_init.h"

#include <boost/filesystem.hpp>

#include "ReadExperiment.hpp"
#include "SalmonOpts.hpp"
#include "GZipWriter.hpp"
#include "CollapsedEMOptimizer.hpp"
#include "EquivalenceClassBuilder.hpp"

#include "AlevinOpts.hpp"
#include "SingleCellProtocols.hpp"
#include "WhiteList.hpp"
#include "Dedup.hpp"

#include "cuckoohash_map.hh"
#include "Eigen/Dense"
#include "concurrentqueue.h"

using JqueueT = moodycamel::ConcurrentQueue<uint32_t>;
using eqMapT = cuckoohash_map<TranscriptGroup, SCTGValue, TranscriptGroupHasher>;
using tgrouplabelt = std::vector<uint32_t>;
using tgroupweightvec = std::vector<double>;
namespace bfs = boost::filesystem;
using SCExpT = ReadExperiment<EquivalenceClassBuilder<SCTGValue>>;

bool runPerCellEM(std::vector<std::vector<uint32_t>>& txpGroups,
                  std::vector<std::vector<double>>& txpGroupCombinedWeights,
                  std::vector<uint64_t>& txpGroupCounts,
                  const std::vector<Transcript>& transcripts,
                  uint64_t totalNumFrags, double uniformTxpWeight,
                  std::shared_ptr<spdlog::logger>& jointlog,
                  bfs::path& outDirPath,
                  std::unordered_set<uint32_t>& activeTxps);

void optimizeCell(SCExpT& experiment,
                  std::vector<std::string>& trueBarcodes,
                  std::atomic<uint32_t>& barcode,
                  size_t totalCells,
                  eqMapT& eqMap,
                  std::deque<TranscriptGroup>& orderedTgroup,
                  std::shared_ptr<spdlog::logger>& jointlog,
                  bfs::path& outDir, std::vector<uint32_t>& umiCount,
                  tbb::atomic<uint32_t>& skippedCBcount,
                  bool verbose, GZipWriter& gzw, size_t umiLength, bool noEM,
                  spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                  std::vector<std::vector<double>>& countMatrix);

class CollapsedCellOptimizer {
public:
  using VecType = std::vector<tbb::atomic<double>>;
  using SerialVecType = std::vector<double>;
  CollapsedCellOptimizer();

  template <class ProtocolT>
  bool optimize(SCExpT& readExp,
                AlevinOpts<ProtocolT>& aopt,
                GZipWriter& gzw,
                std::vector<std::string>& trueBarcodes,
                std::vector<uint32_t>& umiCount,
                CFreqMapT& freqCounter);
};

using VecT = CollapsedCellOptimizer::SerialVecType;

#endif // COLLAPSED_CELL_OPTIMIZER_HPP
