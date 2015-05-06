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
#include <limits>

#include <jellyfish/err.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/fstream_default.hpp>
#include <jellyfish/jellyfish.hpp>
#include <sub_commands/dump_main_cmdline.hpp>

static dump_main_cmdline args; // Command line switches and arguments

template<typename iterator>
void dump(iterator& it, std::ostream &out,
          uint64_t lower_count, uint64_t upper_count) {
  if(args.column_flag) {
    char spacer = args.tab_flag ? '\t' : ' ';
    while(it.next()) {
      if(it.val() < lower_count || it.val() > upper_count)
        continue;
      out << it.key() << spacer << it.val() << "\n";
    }
  } else {
    while(it.next()) {
      if(it.val() < lower_count || it.val() > upper_count)
        continue;
      out << ">" << it.val() << "\n" << it.key() << "\n";
    }
  }
}

int dump_main(int argc, char *argv[])
{
  args.parse(argc, argv);
  std::ios::sync_with_stdio(false); // No sync with stdio -> faster

  ofstream_default out(args.output_given ? args.output_arg : 0, std::cout);
  if(!out.good())
    die << "Error opening output file '" << args.output_arg << "'";

  std::ifstream is(args.db_arg);
  if(!is.good())
    die << "Failed to open input file '" << args.db_arg << "'" << jellyfish::err::no;
  jellyfish::file_header header;
  header.read(is);
  jellyfish::mer_dna::k(header.key_len() / 2);

  if(!args.lower_count_given)
    args.lower_count_arg = 0;
  if(!args.upper_count_given)
    args.upper_count_arg = std::numeric_limits<uint64_t>::max();

  if(!header.format().compare(binary_dumper::format)) {
    binary_reader reader(is, &header);
    dump(reader, out, args.lower_count_arg, args.upper_count_arg);
  } else if(!header.format().compare(text_dumper::format)) {
    text_reader reader(is, &header);
    dump(reader, out, args.lower_count_arg, args.upper_count_arg);
  } else {
    die << "Unknown format '" << header.format() << "'";
  }

  out.close();

  return 0;
}
