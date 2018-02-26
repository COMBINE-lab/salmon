#ifndef __WHITE_LIST_HPP__
#define __WHITE_LIST_HPP__


#include <unordered_map>
#include <vector>
#include <string>
#include <unordered_set>
#include <cmath>
#include <cassert>
#include <fstream>
#include <numeric>

#include "SingleCellProtocols.hpp"
#include "BarcodeGroup.hpp"
#include "RapMapUtils.hpp"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include <boost/range/irange.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

namespace alevin {
  namespace whitelist {
    template <typename ProtocolT>
    bool performWhitelisting(AlevinOpts<ProtocolT>& aopt,
                             std::vector<uint32_t>& umiCount,
                             std::vector<std::string>& trueBarcodes,
                             CFreqMapT& freqCounter,
                             spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                             spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                             size_t numLowConfidentBarcode);
  }
}

#endif // __WHITE_LIST_HPP__
