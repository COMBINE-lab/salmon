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
//#include "RapMapUtils.hpp"
#include "SingleCellProtocols.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include <boost/range/irange.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

namespace alevin {
  namespace whitelist {

    using BlockedIndexRange = tbb::blocked_range<size_t>;
    using DoubleMatrixT = std::vector<std::vector<double>> ;
    using DoubleVectorT = std::vector<double> ;

    uint32_t populate_count_matrix(boost::filesystem::path& outDir,
                                   size_t numElem,
                                   DoubleMatrixT& countMatrix);

    template <typename ProtocolT>
    bool performWhitelisting(AlevinOpts<ProtocolT>& aopt,
                             std::vector<std::string>& trueBarcodes,
                             bool useRibo, bool useMito,
                             size_t numLowConfidentBarcode);
  }
}

#endif // __WHITE_LIST_HPP__
