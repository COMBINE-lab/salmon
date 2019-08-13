#include "pufferfish/ProgOpts.hpp"
#include "SalmonStringUtils.hpp"
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

std::vector<uint8_t> salmon::stringtools::encodeSequenceInSAM(const char* src, size_t len) {
  //uint8_t* target = new uint8_t[static_cast<size_t>(ceil(len / 2.0))]();
  std::vector<uint8_t> target(static_cast<size_t>(ceil(len / 2.0)), 0);
  for (size_t i = 0; i < len; ++i) {
    size_t byte = i >> 1;
    size_t nibble = i & 0x1;
    if (nibble) {
      target[byte] |= charToSamEncode[static_cast<uint8_t>(src[i])];
    } else {
      target[byte] |= (charToSamEncode[static_cast<uint8_t>(src[i])] << 4);
    }
  }
  return target;
}
