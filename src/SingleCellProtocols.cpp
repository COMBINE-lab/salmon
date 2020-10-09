#include "SingleCellProtocols.hpp"


namespace alevin{
  namespace protocols {
      
  std::ostream& operator<<(std::ostream& os, const TagGeometry& tg) {
      os << "Tag geometry description :: [\n";
      os << "\tRead Num : " << static_cast<int32_t>(tg.read_num) << "\n\t[";
      for (size_t i = 0; i < tg.substr_locs.size(); ++i) {
        os << "  [ p:" << tg.substr_locs[i].first << ", l:" << tg.substr_locs[i].second << ")";
      }
      os << "\t]\n";
      os << "\tTotal Length : " << tg.length;
      os << "]\n";
      return os;
    }

  }// protocols
}//alevin