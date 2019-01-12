#include <unordered_set>

#include "AlevinOpts.hpp"
#include "SingleCellProtocols.hpp"
#include "BarcodeGroup.hpp"
#include "GZipWriter.hpp"
#include "CollapsedCellOptimizer.hpp"

namespace apt = alevin::protocols;
namespace bfs = boost::filesystem;

template <typename ProtocolT>
void alevinOptimize( std::vector<std::string>& trueBarcodesVec,
                     std::vector<Transcript>& transcripts,
                     EqMapT& fullEqMap,
                     AlevinOpts<ProtocolT>& aopt,
                     GZipWriter& gzw,
                     CFreqMapT& freqCounter,
                     size_t numLowConfidentBarcode);

template <typename ProtocolT>
int salmonHashQuantify(AlevinOpts<ProtocolT>& aopt,
                       CFreqMapT& freqCounter) {
  aopt.jointLog->info("Reading BFH");
  aopt.jointLog->flush();

  bfs::path eqFilePath = aopt.bfhFile;
  std::ifstream equivFile(eqFilePath.string());

  //  auto& transcripts = experiment.transcripts();
  //  const auto& eqVec =
  //    experiment.equivalenceClassBuilder().eqMap().lock_table();

  size_t num_txps, num_bcs, num_eqclasses;

  // Number of transcripts
  equivFile >> num_txps;

  // Number of barcodes
  equivFile >> num_bcs;

  // Number of equivalence classes
  equivFile >> num_eqclasses;

  std::vector<std::string> txp_names (num_txps);
  for (size_t i=0; i<num_txps; i++) {
    equivFile >> txp_names[i] ;
  }

  size_t bc_length {aopt.protocol.barcodeLength};
  std::vector<std::string> bc_names (num_bcs);
  for (size_t i=0; i<num_bcs; i++) {
    equivFile >> bc_names[i] ;
    if (bc_names[i].size() != bc_length) {
      aopt.jointLog->error("CB {} has wrong length", bc_names[i]);
      aopt.jointLog->flush();
      exit(1);
    }
  }

  //alevin::types::AlevinUMIKmer umiObj;

  //for (size_t i=0; i<num_eqclasses; i++) {
  //  uint64_t count = eq.second.count;
  //  // for each transcript in this class
  //  const TranscriptGroup& tgroup = eq.first;
  //  const std::vector<uint32_t>& txps = tgroup.txps;

  //  // group size
  //  equivFile << txps.size() << '\t';
  //  // each group member
  //  for (auto tid : txps) { equivFile << tid << '\t'; }
  //  const auto& bgroup = eq.second.barcodeGroup;
  //  equivFile << count << "\t" << bgroup.size();
  //  for (auto  bcIt : bgroup){
  //    auto bc = bcIt.first;
  //    auto ugroup = bcIt.second;
  //    equivFile << "\t" << bc << "\t" << ugroup.size();
  //    for (auto umiIt : ugroup){
  //      auto umi = umiIt.first;
  //      umiObj.word__(0) = umi;
  //      auto count = umiIt.second;

  //      std::string s = umiObj.toStr();
  //      std::reverse(s.begin(), s.end());
  //      equivFile << "\t" << s << "\t" << count;
  //    }
  //  }
  //  equivFile << "\n";
  //}

  equivFile.close();
  return 0;
}


template
int salmonHashQuantify(AlevinOpts<apt::Chromium>& aopt,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::ChromiumV3>& aopt,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::Gemcode>& aopt,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::DropSeq>& aopt,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::InDrop>& aopt,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::CELSeq>& aopt,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::CELSeq2>& aopt,
                       CFreqMapT& freqCounter);
template
int salmonHashQuantify(AlevinOpts<apt::Custom>& aopt,
                       CFreqMapT& freqCounter);
