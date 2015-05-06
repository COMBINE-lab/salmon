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

#include <vector>

#include <jellyfish/err.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/file_header.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <jellyfish/fstream_default.hpp>
#include <jellyfish/jellyfish.hpp>
#include <sub_commands/query_main_cmdline.hpp>

using jellyfish::mer_dna;
using jellyfish::mer_dna_bloom_counter;
typedef std::vector<const char*> file_vector;
typedef jellyfish::mer_overlap_sequence_parser<jellyfish::stream_manager<file_vector::iterator> > sequence_parser;
typedef jellyfish::mer_iterator<sequence_parser, mer_dna> mer_iterator;

static query_main_cmdline args;

// mer_dna_bloom_counter query_load_bloom_filter(const char* path) {
//   return res;
// }

template<typename PathIterator, typename Database>
void query_from_sequence(PathIterator file_begin, PathIterator file_end, const Database& db,
                         std::ostream& out, bool canonical) {
  jellyfish::stream_manager<PathIterator> streams(file_begin, file_end);
  sequence_parser parser(mer_dna::k(), 1, 3, 4096, streams);
  for(mer_iterator mers(parser, canonical); mers; ++mers)
    out << *mers << " " << db.check(*mers) << "\n";
}

template<typename Database>
void query_from_cmdline(std::vector<const char*> mers, const Database& db, std::ostream& out,
                        bool canonical) {
  mer_dna m;
  for(auto it = mers.cbegin(); it != mers.cend(); ++it) {
    try {
      m = *it;
      if(canonical)
        m.canonicalize();
      out << m << " " << db.check(m) << "\n";
    } catch(std::length_error e) {
      std::cerr << "Invalid mer '" << *it << "'\n";
    }
  }
}

template<typename Database>
void query_from_stdin(const Database& db, std::ostream& out, bool canonical) {
  std::string buffer;
  mer_dna     m;

  while(getline(std::cin, buffer)) {
    try {
      m = buffer;
      if(canonical)
        m.canonicalize();
      out << db.check(m) << std::endl;  // a flush is need for interactive use
    } catch(std::length_error e) {
      std::cerr << "Invalid mer '" << buffer << "'" << std::endl;
    }
  }
}

int query_main(int argc, char *argv[])
{
  args.parse(argc, argv);

  ofstream_default out(args.output_given ? args.output_arg : 0, std::cout);
  if(!out.good())
    die << "Error opening output file '" << args.output_arg << "'";

  std::ifstream in(args.file_arg, std::ios::in|std::ios::binary);
  jellyfish::file_header header(in);
  if(!in.good())
    die << "Failed to parse header of file '" << args.file_arg << "'";
  mer_dna::k(header.key_len() / 2);
  if(header.format() == "bloomcounter") {
    jellyfish::hash_pair<mer_dna> fns(header.matrix(1), header.matrix(2));
    mer_dna_bloom_counter filter(header.size(), header.nb_hashes(), in, fns);
    if(!in.good())
      die << "Bloom filter file is truncated";
    in.close();
    query_from_sequence(args.sequence_arg.begin(), args.sequence_arg.end(), filter, out, header.canonical());
    query_from_cmdline(args.mers_arg, filter, out, header.canonical());
    if(args.interactive_flag)  query_from_stdin(filter, out, header.canonical());
  } else if(header.format() == binary_dumper::format) {
    jellyfish::mapped_file binary_map(args.file_arg);
    binary_query bq(binary_map.base() + header.offset(), header.key_len(), header.counter_len(), header.matrix(),
                               header.size() - 1, binary_map.length() - header.offset());
    query_from_sequence(args.sequence_arg.begin(), args.sequence_arg.end(), bq, out, header.canonical());
    query_from_cmdline(args.mers_arg, bq, out, header.canonical());
    if(args.interactive_flag)  query_from_stdin(bq, out, header.canonical());
  } else {
    die << "Unsupported format '" << header.format() << "'. Must be a bloom counter or binary list.";
  }

  return 0;
}
