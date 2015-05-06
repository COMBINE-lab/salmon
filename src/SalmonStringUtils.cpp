#include "SalmonStringUtils.hpp"
#include <iostream>
#include <cstdint>
#include <cmath>

uint8_t* salmon::stringtools::encodeSequenceInSAM(const char* src, size_t len) {
    uint8_t* target = new uint8_t[static_cast<size_t>(ceil(len / 2.0))]();
    for(size_t i = 0; i < len; ++i) {
        size_t byte = i >> 1;
        size_t nibble = i & 0x1;
        if (nibble) {
            target[byte] |= charToSamEncode[src[i]];
        } else {
            target[byte] |= (charToSamEncode[src[i]] << 4);
        }
    }
    return target;
}

