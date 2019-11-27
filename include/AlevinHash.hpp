#ifndef __ALEVIN_HASH_HPP__
#define __ALEVIN_HASH_HPP__

#include <unordered_set>
#include <unordered_map>

#include "spdlog/spdlog.h"

#include "AlevinOpts.hpp"
#include "AlevinUtils.hpp"
#include "SingleCellProtocols.hpp"
#include "GZipWriter.hpp"
#include "TranscriptGroup.hpp"
#include "CollapsedCellOptimizer.hpp"
#include "BarcodeGroup.hpp"
#include "SalmonIndex.hpp"

namespace apt = alevin::protocols;
namespace bfs = boost::filesystem;

template <typename ProtocolT>
void alevinOptimize( std::vector<std::string>& trueBarcodesVec,
                     spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                     spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                     EqMapT& fullEqMap,
                     AlevinOpts<ProtocolT>& aopt,
                     GZipWriter& gzw,
                     CFreqMapT& freqCounter,
                     size_t numLowConfidentBarcode);

template <typename ProtocolT>
int salmonHashQuantify(AlevinOpts<ProtocolT>& aopt,
                       bfs::path& outputDirectory,
                       CFreqMapT& freqCounter);
#endif // __ALEVIN_HASH_HPP__
