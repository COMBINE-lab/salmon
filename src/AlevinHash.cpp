#include <unordered_set>

#include <AlevinOpts.hpp>
#include "SingleCellProtocols.hpp"
#include "BarcodeGroup.hpp"

namespace apt = alevin::protocols;

template <typename ProtocolT>
int salmonHashQuantify(AlevinOpts<ProtocolT>& aopt,
                       CFreqMapT& freqCounter) {
  aopt.jointLog->info("Reading Hash");
  aopt.jointLog->flush();
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
