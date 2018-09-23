#ifndef DEDUP_UMI_HPP
#define DEDUP_UMI_HPP

#include "edlib.h"
#include "AlevinUtils.hpp"
#include <boost/graph/adjacency_list.hpp>

using UGroupT = spp::sparse_hash_map<uint64_t, uint32_t>;
using TGroupT = std::vector<uint32_t>;
bool dedupClasses(std::vector<TGroupT> txpGroups,
                  std::vector<UGroupT> umiGroups,
                  spp::sparse_hash_map<uint32_t, uint32_t> &txpToGeneMap);

#endif // DEDUP_HPP
