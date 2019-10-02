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
#include "EquivalenceClassBuilder.hpp"

#include "AlevinOpts.hpp"
#include "SingleCellProtocols.hpp"
#include "WhiteList.hpp"
#include "DedupUMI.hpp"
#include "TranscriptGroup.hpp"

#include "cuckoohash_map.hh"
#include "Eigen/Dense"
#include "concurrentqueue.h"

#include <boost/math/special_functions/digamma.hpp>

namespace bfs = boost::filesystem;
using JqueueT = moodycamel::ConcurrentQueue<uint32_t>;
using eqMapT = cuckoohash_map<TranscriptGroup, SCTGValue, TranscriptGroupHasher>;
using tgrouplabelt = std::vector<uint32_t>;
using tgroupweightvec = std::vector<double>;
using SCExpT = ReadExperiment<EquivalenceClassBuilder<SCTGValue>>;
using EqMapT = cuckoohash_map<TranscriptGroup, SCTGValue, TranscriptGroupHasher>;

constexpr double digammaMin = 1e-10;

struct CellState {
  bool inActive;
};

class CollapsedCellOptimizer {
public:
  using VecType = std::vector<tbb::atomic<double>>;
  using SerialVecType = std::vector<double>;
  CollapsedCellOptimizer();

  template <class ProtocolT>
  bool optimize(EqMapT& freqMap,
                spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                AlevinOpts<ProtocolT>& aopt,
                GZipWriter& gzw,
                std::vector<std::string>& trueBarcodes,
                std::vector<uint32_t>& umiCount,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode);
};

bool runPerCellEM(double& totalNumFrags, size_t numGenes,
                  CollapsedCellOptimizer::SerialVecType& alphas,
                  const CollapsedCellOptimizer::SerialVecType& priorAlphas,
                  std::vector<SalmonEqClass>& salmonEqclasses,
                  std::shared_ptr<spdlog::logger>& jointlog,
                  bool initUniform, bool useVBEM);

void optimizeCell(std::vector<std::string>& trueBarcodes,
                  const std::vector<std::vector<double>>& priorAlphas,
                  std::atomic<uint32_t>& barcode,
                  size_t totalCells, uint32_t umiEditDistance, eqMapT& eqMap,
                  std::deque<std::pair<TranscriptGroup, uint32_t>>& orderedTgroup,
                  std::shared_ptr<spdlog::logger>& jointlog,
                  std::vector<uint32_t>& umiCount,
                  std::vector<CellState>& skippedCB,
                  bool verbose, GZipWriter& gzw, bool noEM, bool useVBEM,
                  bool quiet, tbb::atomic<double>& totalDedupCounts,
                  tbb::atomic<uint32_t>& totalExpGeneCounts,
                  spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                  uint32_t numGenes, uint32_t umiLength, uint32_t numBootstraps,
                  bool naiveEqclass, bool dumpUmiGraph, bool useAllBootstraps,
                  bool initUniform, CFreqMapT& freqCounter,
                  spp::sparse_hash_set<uint32_t>& mRnaGenes,
                  spp::sparse_hash_set<uint32_t>& rRnaGenes,
                  std::atomic<uint64_t>& totalUniEdgesCounts,
                  std::atomic<uint64_t>& totalBiEdgesCounts);

using VecT = CollapsedCellOptimizer::SerialVecType;

#endif // COLLAPSED_CELL_OPTIMIZER_HPP
