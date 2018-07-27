#ifndef __ALEVIN_TYPES_HPP__
#define __ALEVIN_TYPES_HPP__

#include "rapmap/Kmer.hpp"

namespace alevin {
namespace types {

  using AlevinUMIKmer = combinelib::kmers::Kmer<32,2>;

}
}

#endif//__ALEVIN_TYPES_HPP__

