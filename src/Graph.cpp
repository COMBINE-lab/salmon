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

    // returns true if distance is 1
    // returns false if distance is > 1
    bool oneHamming(uint64_t k1, uint64_t k2) {
      auto k1XORk2 = k1 ^ k2;
      int pcnt = __builtin_popcountl(k1XORk2);
      if (pcnt == 1) { return true; }
      if (pcnt > 2) { return false; }
      int lzero = __builtin_clzl(k1XORk2);
      // if it's odd, then substract 1, otherwise keep the shift the same
      int lshift = (lzero & 1) ? (lzero - 1) : lzero;
      return (((k1XORk2 << lshift) & 0x3FFFFFFFFFFFFFFF) > 0) ? false : true;
    }

    EdgeType hasEdge(std::pair<uint64_t, uint32_t> &x,
                     std::pair<uint64_t, uint32_t> &y) {
      if ( x.first == y.first ) {
        return EdgeType::BiDirected;
      }
      if ( x.second > (2*y.second-1) ) {
        if ( oneHamming(x.first, y.first) ){
          return EdgeType::XToY;
        }
        else {
          return EdgeType::NoEdge;
        }
      }
      else if (y.second > (2*x.second-1) ) {
        if ( oneHamming(x.first, y.first) ){
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
