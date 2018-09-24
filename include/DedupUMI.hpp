#ifndef DEDUP_UMI_HPP
#define DEDUP_UMI_HPP

#include "Graph.hpp"
#include "AlevinUtils.hpp"

struct SalmonEqClass {
  std::vector<uint32_t> labels;
  uint32_t count;
};

using UGroupT = spp::sparse_hash_map<uint64_t, uint32_t>;
using TGroupT = std::vector<uint32_t>;

bool dedupClasses(std::vector<double>& geneAlphas,
                  uint64_t& totalUMICount,
                  std::vector<TGroupT>& txpGroups,
                  std::vector<UGroupT>& umiGroups,
                  std::vector<SalmonEqClass>& salmonEqclasses,
                  spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap);

#endif // DEDUP_UMI_HPP
