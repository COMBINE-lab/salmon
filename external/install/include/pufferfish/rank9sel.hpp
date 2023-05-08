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

#ifndef rank9sel_h
#define rank9sel_h
#include <stdint.h>

#include "macros.hpp"
#include "select.hpp"
#include "compact_vector/compact_vector.hpp"

class rank9sel {
private:
	//const compact::vector<uint64_t, 1>* bits;
	uint64_t* bits;
	uint64_t *counts, *inventory, *subinventory;
	uint64_t num_words, num_counts, inventory_size, ones_per_inventory, log2_ones_per_inventory, num_ones;

public:
  rank9sel ();
	rank9sel( compact::vector<uint64_t, 1>* compact_bits, uint64_t num_bits );
  rank9sel( rank9sel&& other);
  rank9sel( rank9sel& other) = delete;

  rank9sel& operator=(rank9sel&& other);
  rank9sel& operator=(rank9sel& other) = delete;

	~rank9sel();
	void arm_hack() const;
	uint64_t rank( const uint64_t pos ) const;
	uint64_t select( const uint64_t rank ) const;
	uint64_t get_word(const uint64_t index);
	// Just for analysis purposes
	void print_counts();
	uint64_t bit_count();
};

#endif
