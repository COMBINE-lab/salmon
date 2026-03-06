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

template <typename T> using ReadGroup = salmon_fqfeeder::ReadGroup<T>;

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

template <typename T> class FastxReader {
public:
  FastxReader(std::vector<std::string> files, uint32_t numConsumers,
              uint32_t numParsers = 1, uint32_t chunkSize = 1000)
      : config_{numConsumers, numParsers, chunkSize, false},
        parser_(config_, std::move(files)) {}

  FastxReader(std::vector<std::string> files, std::vector<std::string> files2,
              uint32_t numConsumers, uint32_t numParsers = 1,
              uint32_t chunkSize = 1000)
      : config_{numConsumers, numParsers, chunkSize, false},
        parser_(config_, std::move(files), std::move(files2)) {}

  bool start() { return parser_.start(); }
  bool stop() { return parser_.stop(); }

  ReadGroup<T> getReadGroup() { return parser_.getReadGroup(); }
  bool refill(ReadGroup<T>& rg) { return parser_.refill(rg); }
  void finishedWithGroup(ReadGroup<T>& rg) { parser_.finishedWithGroup(rg); }

private:
  salmon_fqfeeder::ParserConfig config_;
  salmon_fqfeeder::FastxParser<T> parser_;
};

} // namespace salmon::io::fastx
