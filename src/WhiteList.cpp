#include "WhiteList.hpp"

namespace alevin {
  namespace whitelist {
    template <typename ProtocolT>
    bool performWhitelisting(AlevinOpts<ProtocolT>& aopt,
                             std::vector<uint32_t>& umiCount,
                             std::vector<std::string>& trueBarcodes,
                             CFreqMapT& freqCounter){
      // freqCounter has frequency of reads for all detected Barcodes
      // umiCount has frequency of UMI per-cell after knee selection
      // Have to read the sparse matrix file from dish for deduplicated
      // counts
      return false;
    }
    template bool performWhitelisting(AlevinOpts<alevin::protocols::DropSeq>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::InDrop>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Chromium>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter);
    template bool performWhitelisting(AlevinOpts<alevin::protocols::Custom>& aopt,
                                      std::vector<uint32_t>& umiCount,
                                      std::vector<std::string>& trueBarcodes,
                                      CFreqMapT& freqCounter);
  }
}
