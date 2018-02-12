#ifndef __WHITE_LIST_HPP__
#define __WHITE_LIST_HPP__


#include <unordered_map>
#include <vector>
#include <string>
#include <unordered_set>

#include "SingleCellProtocols.hpp"
#include "BarcodeGroup.hpp"

namespace alevin {
  namespace whitelist {
    template <typename ProtocolT>
    bool performWhitelisting(AlevinOpts<ProtocolT>& aopt,
                             std::vector<uint32_t>& umiCount,
                             std::vector<std::string>& trueBarcodes,
                             CFreqMapT& freqCounter);
  }
}

#endif // __WHITE_LIST_HPP__
