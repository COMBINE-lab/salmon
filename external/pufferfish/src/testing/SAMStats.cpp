#include "FastxParser.hpp"
#include "sparsepp/spp.h"
#include <iostream>
#include <fstream>
#include <vector>

// https://stackoverflow.com/questions/236129/the-most-elegant-way-to-iterate-the-words-of-a-string
template<typename T>
std::vector<T>
split(const T & str, const T & delimiters) {
  std::vector<T> v;
  typename T::size_type start = 0;
  auto pos = str.find_first_of(delimiters, start);
  while(pos != T::npos) {
    if(pos != start) // ignore empty tokens
      v.emplace_back(str, start, pos - start);
    start = pos + 1;
    pos = str.find_first_of(delimiters, start);
  }
  if(start < str.length()) // ignore trailing delimiter
    v.emplace_back(str, start, str.length() - start); // add what's left of the string
  return v;
}

int main(int argc, char* argv[]) {

  std::ifstream sfile(argv[1]);
  spp::sparse_hash_map<std::string, bool> mappedReads;

  std::string tab = "\t";
  std::string slash= "/";
  std::string semicolon = ";";
  size_t i = 0;
  std::string prevName = "";
  size_t correctNum{0};
  bool haveCorrectMapping = false;
  for (std::string line; std::getline(sfile, line);) {
    // skip header lines
    if (line[0] != '@')  {
      std::vector<std::string> toks = split(line, tab);
      std::vector<std::string> ntok = split(toks[0], slash);
      std::vector<std::string> trueTxps = split(ntok[1], semicolon);
      auto& readName = ntok[0];
      auto& trueTxp = trueTxps[0];
      auto& mappedTxp = toks[2];
      if (readName != prevName) {
        prevName = readName;
        haveCorrectMapping = false;
      }
      if (!haveCorrectMapping) {
        bool foundTruth = (trueTxp == mappedTxp);
        if (foundTruth) {
          mappedReads[readName] = true;
          haveCorrectMapping = true;
          ++correctNum;
        } else {
          mappedReads[readName] = false;
        }
      }
      ++i;
      if (i % 1000000 == 0) {
        std::cerr << "processed " << mappedReads.size() << " reads; " << correctNum << " have been correctly mapped\n";
      }
    }
  }


  std::cerr << "\n\n total mapped = " << mappedReads.size() << "\n\n";
  size_t correct{0};
  for (auto& kv : mappedReads) {
    correct += kv.second ? 1 : 0;
  }
  std::cerr << "\n\n correctly mapped = " << correct << "\n\n";
  return 0;
}
