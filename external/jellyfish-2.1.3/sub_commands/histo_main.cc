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
#include <vector>

#include <jellyfish/err.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/fstream_default.hpp>
#include <jellyfish/jellyfish.hpp>
#include <sub_commands/histo_main_cmdline.hpp>

template<typename reader_type>
void compute_histo(reader_type& reader, const uint64_t base, const uint64_t ceil,
                   uint64_t* histo, const uint64_t nb_buckets, const uint64_t inc) {
  while(reader.next()) {
    if(reader.val() < base)
      ++histo[0];
    else if(reader.val() > ceil)
      ++histo[nb_buckets - 1];
    else
      ++histo[(reader.val() - base) / inc];
  }
}


int histo_main(int argc, char *argv[])
{
  histo_main_cmdline args(argc, argv);

  std::ifstream is(args.db_arg);
  if(!is.good())
    die << "Failed to open input file '" << args.db_arg << "'" << jellyfish::err::no;
  jellyfish::file_header header;
  header.read(is);
  jellyfish::mer_dna::k(header.key_len() / 2);

  if(args.high_arg < args.low_arg)
    histo_main_cmdline::error("High count value must be >= to low count value");
  ofstream_default out(args.output_given ? args.output_arg : 0, std::cout);
  if(!out.good())
    die << "Error opening output file '" << args.output_arg << "'" << jellyfish::err::no;

  const uint64_t base = args.increment_arg >= args.low_arg ? 0 : args.low_arg - args.increment_arg;
  const uint64_t ceil = args.high_arg + args.increment_arg;
  const uint64_t inc  = args.increment_arg;

  const uint64_t nb_buckets  = (ceil + inc - base) / inc;
  uint64_t*      histo       = new uint64_t[nb_buckets];
  memset(histo, '\0', sizeof(uint64_t) * nb_buckets);

  if(!header.format().compare(binary_dumper::format)) {
    binary_reader reader(is, &header);
    compute_histo(reader, base, ceil, histo, nb_buckets, inc);
  } else if(!header.format().compare(text_dumper::format)) {
    text_reader reader(is, &header);
    compute_histo(reader, base, ceil, histo, nb_buckets, inc);
  } else {
    die << "Unknown format '" << header.format() << "'";
  }

  for(uint64_t i = 0, col = base; i < nb_buckets; ++i, col += inc)
    if(histo[i] > 0 || args.full_flag)
      out << col << " " << histo[i] << "\n";

  delete [] histo;
  out.close();

  return 0;
}
