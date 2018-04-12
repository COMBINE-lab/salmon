#ifndef DEDUP_HPP
#define DEDUP_HPP

#include "ReadExperiment.hpp"
#include "SingleCellProtocols.hpp"
#include "edlib.h"
#include "SalmonOpts.hpp"
#include "CollapsedEMOptimizer.hpp"
#include "GZipWriter.hpp"
#include "AlevinKmer.hpp"
#include "RapMapUtils.hpp"
#include "AlevinUtils.hpp"

namespace apt = alevin::protocols;
using UGroupT = spp::sparse_hash_map<uint64_t, uint32_t>;

uint32_t edLibCollapse(const spp::sparse_hash_set<uint64_t>& umiList,
                       std::vector<uint64_t>& vList,
                       const size_t& length,
                       const UGroupT& ugroup,
                       AlignerEngine& ae,
                       std::vector<std::string>& umiSeqs);

uint32_t neighborCollapse(const spp::sparse_hash_set<uint64_t>& umiList,
                          std::vector<uint64_t>& vList,
                          const size_t& length,
                          const UGroupT& ugroup,
                          std::vector<std::string>& umiSeqs);

uint32_t dedupReads(
                    const size_t umiLength, std::shared_ptr<spdlog::logger>& jointlog,
                    const UGroupT& ugroup,
                    spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint64_t>>& umiBiasList,
                    const TranscriptGroup& tgroup);

#endif // DEDUP_HPP
