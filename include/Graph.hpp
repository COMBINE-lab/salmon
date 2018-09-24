#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <boost/functional/hash.hpp>
#include "AlevinUtils.hpp"
#include "edlib.h"

namespace alevin {
  namespace graph {

    using VertexT = std::pair<uint32_t, uint32_t>;

    enum EdgeType {
      NoEdge,
      BiDirected,
      XToY,
      YToX,
    };

    struct Graph {
      std::vector<VertexT> vertexNames;
      spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t>> edges;

      void add_edge(uint32_t source, uint32_t sink) {
        edges[source].insert(sink);
      }
    };

    uint32_t getVertexIndex(spp::sparse_hash_map<VertexT, uint32_t, boost::hash<VertexT>>& vertMap,
                            VertexT& node);
    EdgeType hasEdge(std::pair<std::string, uint32_t> &x,
                     std::pair<std::string, uint32_t> &y,
                     AlignerEngine& ae);
  }
}

#endif // GRAPH_HPP
