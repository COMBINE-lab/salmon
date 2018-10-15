#include "Graph.hpp"

namespace alevin {
  namespace graph {
    uint32_t getVertexIndex(spp::sparse_hash_map<VertexT, uint32_t, boost::hash<VertexT>>& vertMap,
                            VertexT& node){
      uint32_t index;

      if ( vertMap.contains(node) ) {
        index = vertMap[node];
      }
      else{
        index = vertMap.size();
        vertMap[node] = index;
      }

      return index;
    }

    bool is_one_edit(std::string& first,
                     std::string& second) {
      size_t seqLen = first.size();
      uint32_t distance = 0;
      for(size_t i=0; i<seqLen; i++) {
        if (first[i] != second[i]) {
          distance += 1;
          if (distance > 1) {
            return false;
          }
        }
      }

      return true;
    }

    EdgeType hasEdge(std::pair<std::string, uint32_t> &x,
                     std::pair<std::string, uint32_t> &y) {
      if ( x.first.compare(y.first) == 0 ) {
        return EdgeType::BiDirected;
      }
      if ( x.second > (2*y.second-1) ) {
        if ( is_one_edit(x.first, y.first) ){
          return EdgeType::XToY;
        }
        else {
          return EdgeType::NoEdge;
        }
      }
      else if (y.second > (2*x.second-1) ) {
        if ( is_one_edit(x.first, y.first) ){
          return EdgeType::YToX;
        }
        else {
          return EdgeType::NoEdge;
        }
      }
      else{
        return EdgeType::NoEdge;
      }
    }
  }
}
