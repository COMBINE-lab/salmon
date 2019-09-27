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

    // NOTE : I like our clever tricks here, but we should
    // make note of how widely distributable the compiled code
    // will be (i.e. which instructions will be used for popcountl
    // and clzl, and which processors will they work on).
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

    bool kHamming(uint64_t k1, uint64_t k2, size_t maxDist, uint32_t umiLength) {
      size_t dist {0};
      for (size_t i=0; i<2*umiLength; i+=2) {
        auto pcnt1 = ( k1 << i) & 0xC000000000000000;
        auto pcnt2 = ( k2 << i) & 0xC000000000000000;
        if (pcnt1 != pcnt2) { dist++; }
        if (dist > maxDist) { return false; }
      }
      return true;
    }

    EdgeType hasEdge(std::pair<uint64_t, uint32_t> &x,
                     std::pair<uint64_t, uint32_t> &y,
                     uint32_t umiEditDistance,
                     uint32_t umiLength) {
      if (x.first == y.first) { return EdgeType::BiDirected; }

      bool isCollapsable;
      if (umiEditDistance == 1) {
        isCollapsable = oneHamming(x.first, y.first);
      } else {
        isCollapsable = kHamming(x.first, y.first, umiEditDistance, umiLength);
      }

      if ((x.second > (2*y.second -1)) and isCollapsable) {
        return EdgeType::XToY ;
      } else if ((y.second > (2*x.second-1)) and isCollapsable) {
        return EdgeType::YToX ;
      } else if (isCollapsable) {
        return EdgeType::BiDirected;
      }
      return EdgeType::NoEdge;
    }
  }
}
