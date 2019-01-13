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

  size_t numTxps, numBcs, numEqclasses;

  // Number of transcripts
  equivFile >> numTxps;

  // Number of barcodes
  equivFile >> numBcs;

  // Number of equivalence classes
  equivFile >> numEqclasses;

  std::vector<std::string> txpNames (numTxps);
  for (size_t i=0; i<numTxps; i++) {
    equivFile >> txpNames[i] ;
  }

  size_t bcLength {aopt.protocol.barcodeLength};
  std::vector<std::string> bcNames (numBcs);
  for (size_t i=0; i<numBcs; i++) {
    equivFile >> bcNames[i] ;
    if (bcNames[i].size() != bcLength) {
      aopt.jointLog->error("CB {} has wrong length", bcNames[i]);
      aopt.jointLog->flush();
      exit(1);
    }
  }

  alevin::types::AlevinUMIKmer umiObj;
  //printing on screen progress
  const char RESET_COLOR[] = "\x1b[0m";
  char green[] = "\x1b[30m";
  green[3] = '0' + static_cast<char>(fmt::GREEN);
  char red[] = "\x1b[30m";
  red[3] = '0' + static_cast<char>(fmt::RED);
  std::cerr<<std::endl;

  for (size_t i=0; i<numEqclasses; i++) {
    uint64_t count;
    size_t labelSize ;
    equivFile >> labelSize;

    std::vector<uint32_t> txps(labelSize);
    for (auto& tid : txps) { equivFile >> tid; }

    size_t bgroupSize;
    equivFile >> count >> bgroupSize;
    SparseBarcodeMapType bgroup;

    for (size_t j=0; j<bgroupSize; j++){
      uint32_t bc;
      size_t ugroupSize;

      equivFile >> bc >> ugroupSize;
      auto& ugroup = bgroup[bc];

      for (size_t k=0; k<ugroupSize; k++){
        std::string umiSeq;
        uint64_t umiIndex;
        uint32_t umiCount;
        equivFile >> umiSeq >> umiCount;

        bool isUmiIdxOk = umiObj.fromChars(umiSeq);
        if(isUmiIdxOk){
          umiIndex = umiObj.word(0);
        } else {
          aopt.jointLog->error("Umi Alevin Object conversion error");
          aopt.jointLog->flush();
          exit(1);
        }

        ugroup[umiIndex] = umiCount;
      }// end-ugroup for
    }//end-bgroup for

    double completionFrac = i*100.0/numEqclasses;
    uint32_t percentCompletion {static_cast<uint32_t>(completionFrac)};
    if ( percentCompletion % 10 == 0 || percentCompletion > 95) {
      fmt::print(stderr, "\r{}Done Reading : {}{}%{}",
                 green, red, percentCompletion, RESET_COLOR);
    }
  }
  std::cerr<<std::endl;

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
