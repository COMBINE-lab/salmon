#ifndef DEDUP_UMI_HPP
#define DEDUP_UMI_HPP

#include "Graph.hpp"
#include "GZipWriter.hpp"
#include "AlevinUtils.hpp"

struct SalmonEqClass {
  std::vector<uint32_t> labels;
  uint32_t count;
};

using UGroupT = spp::sparse_hash_map<uint64_t, uint32_t>;
using TGroupT = std::vector<uint32_t>;

// takem from https://stackoverflow.com/questions/896155/tr1unordered-set-union-and-intersection/896440
template <typename InIt1, typename InIt2, typename OutIt>
OutIt unordered_set_intersection(InIt1 b1, InIt1 e1, InIt2 b2, InIt2 e2, OutIt out) {
  while (!(b1 == e1)) {
    if (!(std::find(b2, e2, *b1) == e2)) {
      *out = *b1;
      ++out;
    }

    ++b1;
  }

  return out;
}

bool dedupClasses(std::vector<double>& geneAlphas,
                  double& totalUMICount,
                  std::vector<TGroupT>& txpGroups,
                  std::vector<UGroupT>& umiGroups,
                  std::vector<SalmonEqClass>& salmonEqclasses,
                  uint32_t umiLength,
                  spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                  std::vector<uint8_t>& tiers,
                  GZipWriter& gzw, uint32_t umiEditDistance,
                  bool dumpUmiGraph, std::string& trueBarcodeStr,
                  std::vector<spp::sparse_hash_map<uint16_t, uint16_t>>& arboEqClassCount,
                  bool dumpArborescences,
                  std::atomic<uint64_t>& totalUniEdgesCounts,
                  std::atomic<uint64_t>& totalBiEdgesCounts);

#endif // DEDUP_UMI_HPP
