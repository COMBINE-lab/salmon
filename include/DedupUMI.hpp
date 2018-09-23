#ifndef DEDUP_UMI_HPP
#define DEDUP_UMI_HPP

#include "edlib.h"
#include "AlevinUtils.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

enum EdgeType {
  NoEdge,
  BiDirected,
  XToY,
  YToX,
};

struct VertexType {
  std::pair<uint32_t, uint32_t> vertex;
};

using UGroupT = spp::sparse_hash_map<uint64_t, uint32_t>;
using TGroupT = std::vector<uint32_t>;

typedef boost::adjacency_list<boost::listS,
                              boost::vecS,
                              boost::directedS,
                              std::pair<uint32_t, uint32_t>,
                              EdgeType> DirectedGraph;
typedef boost::graph_traits<DirectedGraph>::vertex_descriptor VertexT;
typedef boost::graph_traits<DirectedGraph>::edge_descriptor EdgeT;

bool dedupClasses(std::vector<TGroupT> txpGroups,
                  std::vector<UGroupT> umiGroups,
                  spp::sparse_hash_map<uint32_t, uint32_t> &txpToGeneMap);

#endif // DEDUP_HPP
