// metrohash64.h
//
// The MIT License (MIT)
//
// Copyright (c) 2015 J. Andrew Rogers
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#ifndef METROHASH_METROHASH_64_H
#define METROHASH_METROHASH_64_H

#include <stdint.h>

class MetroHash64
{
public:
    static const uint32_t bits = 64;

    // Constructor initializes the same as Initialize()
    MetroHash64(const uint64_t seed=0);
    
    // Initializes internal state for new hash with optional seed
    void Initialize(const uint64_t seed=0);
    
    // Update the hash state with a string of bytes. If the length
    // is sufficiently long, the implementation switches to a bulk
    // hashing algorithm directly on the argument buffer for speed.
    void Update(const uint8_t * buffer, const uint64_t length);
    
    // Constructs the final hash and writes it to the argument buffer.
    // After a hash is finalized, this instance must be Initialized()-ed
    // again or the behavior of Update() and Finalize() is undefined.
    void Finalize(uint8_t * const hash);
    
    // A non-incremental function implementation. This can be significantly
    // faster than the incremental implementation for some usage patterns.
    static void Hash(const uint8_t * buffer, const uint64_t length, uint8_t * const hash, const uint64_t seed=0);

    // Does implementation correctly execute test vectors?
    static bool ImplementationVerified();

    // test vectors -- Hash(test_string, seed=0) => test_seed_0
    static const char * test_string;
    static const uint8_t test_seed_0[8];
    static const uint8_t test_seed_1[8];

private:
    static const uint64_t k0 = 0xD6D018F5;
    static const uint64_t k1 = 0xA2AA033B;
    static const uint64_t k2 = 0x62992FC1;
    static const uint64_t k3 = 0x30BC5B29;
    
    struct { uint64_t v[4]; } state;
    struct { uint8_t b[32]; } input;
    uint64_t bytes;
    uint64_t vseed;
};


// Legacy 64-bit hash functions -- do not use
void metrohash64_1(const uint8_t * key, uint64_t len, uint32_t seed, uint8_t * out);
void metrohash64_2(const uint8_t * key, uint64_t len, uint32_t seed, uint8_t * out);


#endif // #ifndef METROHASH_METROHASH_64_H
