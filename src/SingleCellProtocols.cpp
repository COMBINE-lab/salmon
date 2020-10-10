#include "SingleCellProtocols.hpp"


namespace alevin{
  namespace protocols {
      
  std::ostream& operator<<(std::ostream& os, const TagGeometry& tg) {

      os << "Tag geometry description :: [\n";
      if (tg.uses_r1()) {
        os << "\tRead Num : 1 \n\t[";
        for (size_t i = 0; i < tg.substr_locs1.size(); ++i) {
          os << "  [ p:" << tg.substr_locs1[i].first
             << ", l:" << tg.substr_locs1[i].second << ")";
        }
        os << "\t]\n";
      }
      if (tg.uses_r1()) {
        os << "\tRead Num : 2 \n\t[";
        for (size_t i = 0; i < tg.substr_locs2.size(); ++i) {
          os << "  [ p:" << tg.substr_locs2[i].first
             << ", l:" << tg.substr_locs2[i].second << ")";
        }
        os << "\t]\n";
      }
      os << "\tTotal Length : " << tg.length1 + tg.length2;
      os << "]\n";
      return os;
    }

  }// protocols
}//alevin