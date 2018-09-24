#ifndef DEDUP_UMI_HPP
#define DEDUP_UMI_HPP

#include "edlib.h"
#include "AlevinUtils.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/property_map/property_map.hpp>

enum EdgeType {
  NoEdge,
  BiDirected,
  XToY,
  YToX,
};

struct SalmonEqClass {
  std::vector<uint32_t> labels;
  uint32_t count;
};
struct VertexType {
  uint32_t eqclassId;
  uint32_t umiId;
};

using UGroupT = spp::sparse_hash_map<uint64_t, uint32_t>;
using TGroupT = std::vector<uint32_t>;

typedef boost::adjacency_list<boost::listS,
                              boost::vecS,
                              boost::directedS,
                              boost::property<boost::vertex_name_t,
                                              VertexType>> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor VertexT;

// taken from https://stackoverflow.com/questions/10405030/c-unordered-map-fail-when-used-with-a-vector-as-key
template <typename Container>
struct container_hash {
  std::size_t operator()(Container const& c) const {
    return boost::hash_range(c.begin(), c.end());
  }
};

bool dedupClasses(std::vector<double>& geneAlphas,
                  uint64_t& totalUMICount,
                  std::vector<TGroupT>& txpGroups,
                  std::vector<UGroupT>& umiGroups,
                  std::vector<SalmonEqClass>& salmonEqclasses,
                  spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap);

#endif // DEDUP_HPP
