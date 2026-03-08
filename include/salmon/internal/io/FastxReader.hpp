#pragma once

#define fastx_parser salmon_fqfeeder
#include <FastxParser.hpp>
#undef fastx_parser
#undef __FASTX_PARSER__

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace salmon::io::fastx {

using ReadSeq = salmon_fqfeeder::ReadSeq;
using ReadPair = salmon_fqfeeder::ReadPair;

template <typename T> using FastxParser = salmon_fqfeeder::FastxParser<T>;

struct ReadPairView {
  klibpp::KSeq& first;
  klibpp::KSeq& second;
};

using CompatReadSeq = klibpp::KSeq;
using CompatReadPair = ReadPairView;

inline ReadPairView toReadPairView(ReadPair& read) {
  return ReadPairView{read.first(), read.second()};
}

inline klibpp::KSeq& toReadSeqView(ReadSeq& read) { return read.first(); }

} // namespace salmon::io::fastx
