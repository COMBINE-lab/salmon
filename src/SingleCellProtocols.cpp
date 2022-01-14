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

  // store regex for reads 1 and 2
  extern std::string CustomGeo::reg[2];
  // store positions of matches for bc and umi
  extern std::vector<int> CustomGeo::b[2], CustomGeo::u[2];
  // bioRead stores the read number for biological read and bioPat stores match pattern number on regex
  extern unsigned int CustomGeo::bioRead, CustomGeo::bioPat; // biological read would be contigous and on only 1 of the read
  extern uint32_t CustomGeo::minBcLen, CustomGeo::maxBcLen;
  extern uint32_t CustomGeo::minUmiLen, CustomGeo::maxUmiLen;
  extern uint32_t CustomGeo::barcodeLength, CustomGeo::umiLength;
  // bool alevin::protocols::CustomGeo::bioReadFound;
  extern bool CustomGeo::bioReadFound;
  extern boost::regex CustomGeo::rgx[2];
  }// protocols
}//alevin