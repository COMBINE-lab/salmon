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

    EdgeType hasEdge(std::pair<std::string, uint32_t> &x,
                     std::pair<std::string, uint32_t> &y,
                     AlignerEngine& ae) {
      if ( x.first.compare(y.first) == 0 ) {
        return EdgeType::BiDirected;
      }
      if ( x.second > (2*y.second-1) ) {
        ae(x.first.c_str(), x.first.size(),
           y.first.c_str(), y.first.size(),
           edlibNewAlignConfig(1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));

        if ( ae.result().editDistance < 2 ){
          return EdgeType::XToY;
        }
        else {
          return EdgeType::NoEdge;
        }
      }
      else if (y.second > (2*x.second-1) ) {
        ae(x.first.c_str(), x.first.size(),
           y.first.c_str(), y.first.size(),
           edlibNewAlignConfig(1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));

        if ( ae.result().editDistance < 2 ){
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
