#ifndef __ALEVIN_UTILS_HPP__
#define __ALEVIN_UTILS_HPP__

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <exception>
#include <functional>
#include <iterator>
#include <map>
#include <mutex>
#include <queue>
#include <random>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <limits>
#include <string>

#include "spdlog/spdlog.h"

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/variant.hpp>
#include <boost/any.hpp>

#include "cereal/archives/json.hpp"
#include "metro/metrohash.h"
#include "nonstd/optional.hpp"

#include "AlevinOpts.hpp"
#include "SingleCellProtocols.hpp"
#include "BarcodeGroup.hpp"
#include "SalmonDefaults.hpp"

#include "SalmonConfig.hpp"
#include "SalmonUtils.hpp"
// #include "IndexHeader.hpp"

#include "pufferfish/sparsepp/spp.h"

namespace alevin{
  namespace utils{

    namespace apt = alevin::protocols;
    namespace bfs = boost::filesystem;

    constexpr uint32_t uint32_max = std::numeric_limits<uint32_t>::max();

    void getIndelNeighbors(
                           const std::string& barcodeSeq,
                           std::unordered_set<uint32_t>& neighbors);

    void findNeighbors(size_t length,
                       const std::string& barcode,
                       std::unordered_set<std::string>& neighbors);

    //template <typename ProtocolT>
    bool sequenceCheck(const std::string& barcode,
                       //AlevinOpts<ProtocolT>& aopt,
                       //std::mutex& iomutex,
                       Sequence seq = Sequence::BARCODE);

    bool recoverBarcode(std::string& sequence);

    void readWhitelist(bfs::path& filePath,
                       TrueBcsT& trueBarcodes);

    template <typename ProtocolT>
    bool processAlevinOpts(AlevinOpts<ProtocolT>& aopt,
                           SalmonOpts& sopt, bool noTgMap,
                           boost::program_options::variables_map& vm);

    template <typename ProtocolT>
    bool extractUMI(std::string& read,
                    std::string& read2,
                    ProtocolT& pt,
                    std::string& umi);

    template <typename ProtocolT>
    std::string* getReadSequence(ProtocolT& pt,
                         std::string& seq,
                         std::string& seq2,
                         std::string& subseq);

    template <typename ProtocolT>
    bool extractBarcode(std::string& read, std::string& read2, ProtocolT& pt, std::string& bc);

    template <typename OrderedOptionsT>
    bool writeCmdInfo(boost::filesystem::path cmdInfoPath,
                      OrderedOptionsT& orderedOptions) {
      std::ofstream os(cmdInfoPath.string());
      cereal::JSONOutputArchive oa(os);
      oa(cereal::make_nvp("salmon_version:", std::string(salmon::version)));
      for (auto& opt : orderedOptions.options) {
        if (opt.value.size() == 1) {
          oa(cereal::make_nvp(opt.string_key, opt.value.front()));
        } else {
          oa(cereal::make_nvp(opt.string_key, opt.value));
        }
      }
      return true;
    }

    void getTxpToGeneMap(spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                         spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                         const std::string& t2gFile, const std::string& refNamesFile,
                         const std::string& refLengthFile,
                         const std::string& headerFile,
                         std::shared_ptr<spdlog::logger>& jointLog,
                         bool noTgMap);

    bool checkSetCoverage(std::vector<std::vector<uint32_t>>& tgroup,
                          std::vector<uint32_t> txps);

    void combinationUtil(std::vector<uint32_t>& arr, int n, int r,
                         int index, std::vector<uint32_t> data,
                         int i, std::vector<std::vector<uint32_t>>& comb);
    bool hasOneGene(const std::vector<uint32_t>& txps, uint32_t& geneId,
                    spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                    const size_t numGenes);
  }
}
#endif // __ALEVIN_UTILS_HPP__
