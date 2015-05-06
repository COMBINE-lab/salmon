/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __JELLYFISH_MER_DNA_BLOOM_COUNTER_HPP_
#define __JELLYFISH_MER_DNA_BLOOM_COUNTER_HPP_

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/bloom_counter2.hpp>
#include <jellyfish/bloom_filter.hpp>
#include <jellyfish/rectangular_binary_matrix.hpp>
#include <jellyfish/misc.hpp>

namespace jellyfish {
template<>
struct hash_pair<mer_dna> {
  RectangularBinaryMatrix m1, m2;

  hash_pair() : m1(8 * sizeof(uint64_t), mer_dna::k() * 2), m2(8 * sizeof(uint64_t), mer_dna::k() * 2) {
    m1.randomize(random_bits);
    m2.randomize(random_bits);
  }

  hash_pair(RectangularBinaryMatrix&& m1_, RectangularBinaryMatrix&& m2_) : m1(m1_), m2(m2_) { }

  void operator()(const mer_dna& k, uint64_t* hashes) const {
    hashes[0] = m1.times(k);
    hashes[1] = m2.times(k);
  }
};

typedef bloom_counter2<mer_dna> mer_dna_bloom_counter;
typedef bloom_counter2_file<mer_dna> mer_dna_bloom_counter_file;
typedef bloom_filter<mer_dna> mer_dna_bloom_filter;
typedef bloom_filter_file<mer_dna> mer_dna_bloom_filter_file;
}

#endif /* __JELLYFISH_MER_DNA_BLOOM_COUNTER_HPP_ */
