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

#ifndef __JELLYFISH_JELLYFISH_HPP__
#define __JELLYFISH_JELLYFISH_HPP__

#include <stdint.h>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/text_dumper.hpp>
#include <jellyfish/binary_dumper.hpp>

typedef jellyfish::cooperative::hash_counter<jellyfish::mer_dna> mer_hash;
typedef mer_hash::array mer_array;
typedef jellyfish::text_dumper<mer_array> text_dumper;
typedef jellyfish::text_reader<jellyfish::mer_dna, uint64_t> text_reader;
typedef jellyfish::binary_dumper<mer_array> binary_dumper;
typedef jellyfish::binary_reader<jellyfish::mer_dna, uint64_t> binary_reader;
typedef jellyfish::binary_query_base<jellyfish::mer_dna, uint64_t> binary_query;
typedef jellyfish::binary_writer<jellyfish::mer_dna, uint64_t> binary_writer;
typedef jellyfish::text_writer<jellyfish::mer_dna, uint64_t> text_writer;


#endif /* __JELLYFISH_JELLYFISH_HPP__ */
