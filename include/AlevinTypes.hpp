#ifndef __ALEVIN_TYPES_HPP__
#define __ALEVIN_TYPES_HPP__

#include "pufferfish/Kmer.hpp"

namespace alevin {
namespace types {

  using AlevinUMIKmer = combinelib::kmers::Kmer<32,2>;
  // Can only hold the barcode if it's <= 32 nucleotides
  using AlevinCellBarcodeKmer = combinelib::kmers::Kmer<32,5>;

}
}

#endif//__ALEVIN_TYPES_HPP__

