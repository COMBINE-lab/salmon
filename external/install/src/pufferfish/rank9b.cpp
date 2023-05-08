/*		 
 * Sux: Succinct data structures
 *
 * Copyright (C) 2007-2013 Sebastiano Vigna 
 *
 *  This library is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This library is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
 *  for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <cassert>
#include <cstring>
#include "rank9b.hpp"

rank9b::rank9b() {}

rank9b::rank9b( const uint64_t * const bits, const uint64_t num_bits ) {
	this->bits = bits;
	num_words = ( num_bits + 63 ) / 64;
	num_counts = ( ( num_bits + 64 * 8 - 1 ) / ( 64 * 8 ) ) * 2;
	
	// Init rank structure
	counts = new uint64_t[ num_counts + 1 ];
	memset( counts, 0, ( num_counts + 1 ) * sizeof *counts );

	uint64_t c = 0;
	uint64_t pos = 0;
	for( uint64_t i = 0; i < num_words; i += 8, pos += 2 ) {
		counts[ pos ] = c;
		c += __builtin_popcountll( bits[ i ] );
		for( int j = 1;  j < 8; j++ ) {
			counts[ pos + 1 ] |= ( c - counts[ pos ] ) << (63 - (9 * j));
			if ( i + j < num_words ) c += __builtin_popcountll( bits[ i + j ] );
		}
	}

	counts[ num_counts ] = c;

	assert( c <= num_bits );
}


rank9b& rank9b::operator=(rank9b&& other) {
  num_words = other.num_words;
  num_counts = other.num_counts;
  inventory_size = other.inventory_size;
  ones_per_inventory = other.ones_per_inventory;
  log2_ones_per_inventory = other.log2_ones_per_inventory;
  num_ones = other.num_ones;
  bits = other.bits;
  counts = other.counts;
  other.counts = nullptr;
  inventory = other.inventory;
  other.inventory = nullptr;
  return *this;
}

rank9b::rank9b(rank9b&& other) {
  num_words = other.num_words;
  num_counts = other.num_counts;
  inventory_size = other.inventory_size;
  ones_per_inventory = other.ones_per_inventory;
  log2_ones_per_inventory = other.log2_ones_per_inventory;
  num_ones = other.num_ones;
  bits = other.bits;
  counts = other.counts;
  other.counts = nullptr;
  inventory = other.inventory;
  other.inventory = nullptr;
}

rank9b::~rank9b() {
  if (counts != nullptr) {
	delete [] counts;
  }
}



uint64_t rank9b::rank( const uint64_t k ) {
	const uint64_t word = k / 64;
	const uint64_t block = word / 4 & ~1;
	const int offset = word % 8;
	return counts[ block ] + ( counts[ block + 1 ] >> ( 63 - offset * 9 ) & 0x1FF ) + __builtin_popcountll( bits[ word ] & ( ( 1ULL << k % 64 ) - 1 ) );
}

uint64_t rank9b::bit_count() {
	return num_counts * 64;
}

void rank9b::print_counts() {}
