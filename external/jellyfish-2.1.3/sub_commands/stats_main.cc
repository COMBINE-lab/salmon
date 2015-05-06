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

#include <config.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <iostream>
#include <fstream>

#include <jellyfish/err.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/fstream_default.hpp>
#include <jellyfish/jellyfish.hpp>
#include <sub_commands/stats_main_cmdline.hpp>

template<typename reader_type>
void compute_stats(reader_type& reader, uint64_t low, uint64_t high,
                   uint64_t& uniq, uint64_t& distinct, uint64_t& total,
                   uint64_t& max) {
  uniq = distinct = total = max = 0;

  while(reader.next()) {
    if(reader.val() < low || reader.val() > high) continue;
    uniq  += reader.val() == 1;
    total += reader.val();
    max    = std::max(max, reader.val());
    ++distinct;
  }
}


int stats_main(int argc, char *argv[])
{
  stats_main_cmdline args(argc, argv);

  std::ifstream is(args.db_arg);
  if(!is.good())
    die << "Failed to open input file '" << args.db_arg << "'" << jellyfish::err::no;
  jellyfish::file_header header;
  header.read(is);
  jellyfish::mer_dna::k(header.key_len() / 2);

  ofstream_default out(args.output_given ? args.output_arg : 0, std::cout);
  if(!out.good())
    die << "Error opening output file '" << args.output_arg << "'" << jellyfish::err::no;

  if(!args.upper_count_given)
    args.upper_count_arg = std::numeric_limits<uint64_t>::max();
  uint64_t uniq = 0, distinct = 0, total = 0, max = 0;
  if(!header.format().compare(binary_dumper::format)) {
    binary_reader reader(is, &header);
    compute_stats(reader, args.lower_count_arg, args.upper_count_arg, uniq, distinct, total, max);
  } else if(!header.format().compare(text_dumper::format)) {
    text_reader reader(is, &header);
    compute_stats(reader, args.lower_count_arg, args.upper_count_arg, uniq, distinct, total, max);
  } else {
    die << "Unknown format '" << header.format() << "'";
  }

  out << "Unique:    " << uniq << "\n"
      << "Distinct:  " << distinct << "\n"
      << "Total:     " << total << "\n"
      << "Max_count: " << max << "\n";
  out.close();

  return 0;
}
