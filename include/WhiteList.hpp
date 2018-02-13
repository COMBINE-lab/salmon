#ifndef __WHITE_LIST_HPP__
#define __WHITE_LIST_HPP__


#include <unordered_map>
#include <vector>
#include <string>
#include <unordered_set>

#include "SingleCellProtocols.hpp"
#include "BarcodeGroup.hpp"
#include "RapMapUtils.hpp"

namespace alevin {
  namespace whitelist {
    template <typename ProtocolT>
    bool performWhitelisting(AlevinOpts<ProtocolT>& aopt,
                             std::vector<uint32_t>& umiCount,
                             std::vector<std::string>& trueBarcodes,
                             CFreqMapT& freqCounter, size_t numGenes,
                             std::vector<std::vector<double>>& countMatrix,
                             spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap);
  }
}

#endif // __WHITE_LIST_HPP__
