#pragma once

#define fastx_parser salmon_fqfeeder
#include "include/FastxParser.hpp"
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
using ParserConfig = salmon_fqfeeder::ParserConfig;

struct CompatReadSeq {
  std::string name;
  std::string comment;
  std::string seq;
  std::string qual;
};

struct CompatReadPair {
  klibpp::KSeq first;
  klibpp::KSeq second;
};

inline CompatReadSeq toCompatReadSeq(const ReadSeq& read) {
  const auto& k = read.first();
  return CompatReadSeq{k.name, k.comment, k.seq, k.qual};
}

inline CompatReadPair toCompatReadPair(const ReadPair& read) {
  return CompatReadPair{read.first(), read.second()};
}

} // namespace salmon::io::fastx
