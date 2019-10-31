#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <boost/graph/connected_components.hpp>
#include <boost/graph/adjacency_list.hpp>
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

      size_t num_vertices() {
        return vertexNames.size();
      }

      size_t num_edges() {
        size_t i = 0;
        for(auto& it: edges) {
          i += it.second.size();
        }
        return i;
      }

      uint32_t getEqclassId(uint32_t vertex) {
        return vertexNames[vertex].first;
      }

      uint32_t getUId(uint32_t vertex) {
        return vertexNames[vertex].second;
      }

      const spp::sparse_hash_set<uint32_t>& getNeighbors(uint32_t vertex) {
        return edges[vertex];
      }

      uint32_t connected_components(std::vector<uint32_t>& component){
        // resize the component based on the number of vertices
        component.resize( vertexNames.size() );

        typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS > AdjList;
        AdjList adjList ( vertexNames.size() );

        // iterating over edges and filling the graph
        for (auto& it: edges) {
          uint32_t source = it.first;
          for(uint32_t target: it.second) {
            boost::add_edge(source, target, adjList);
          }
        }

        return boost::connected_components(adjList, component.data());
      }
    };

    uint32_t getVertexIndex(spp::sparse_hash_map<VertexT, uint32_t, boost::hash<VertexT>>& vertMap,
                            VertexT& node);

    EdgeType hasEdge(std::pair<uint64_t, uint32_t> &x,
                     std::pair<uint64_t, uint32_t> &y,
                     uint32_t umiEditDistance,
                     uint32_t umiLength);
  }
}

#endif // GRAPH_HPP
