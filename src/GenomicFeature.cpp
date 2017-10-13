/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Salmon.

    Salmon is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Salmon is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Salmon.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/

#include "GenomicFeature.hpp"
#include <atomic>
#include <boost/thread/thread.hpp>
#include <fstream>
#include <thread>

std::ostream& operator<<(std::ostream& out, const TranscriptGeneID& ids) {
  out << "transcript_id => " << ids.transcript_id << ", gene_id => "
      << ids.gene_id;
  return out;
}

template <typename StaticAttributes>
std::ostream& operator<<(std::ostream& out,
                         const GenomicFeature<StaticAttributes>& gf) {
  out << "Feature: ["
      << "id :" << gf.seqID << ", "
      << "source : " << gf.source << ", "
      << "type : " << gf.type << ", "
      << "start : " << gf.start << ", "
      << "end : " << gf.end << ", "
      << "score : " << gf.score << ", "
      << "strand : " << gf.strand << ", "
      << "phase : " << gf.phase << " ** ";
  out << gf.sattr << ", ";
  /*for ( auto& kv : gf.attributes ) {
      out << kv.first << " => " << kv.second << ", ";
   }
   */
  return out;
}

template <typename StaticAttributes>
std::istream& operator>>(std::istream& in,
                         GenomicFeature<StaticAttributes>& gf) {
  std::string line;
  std::getline(in, gf.seqID, '\t');
  std::getline(in, gf.source, '\t');
  std::getline(in, gf.type, '\t');

  std::getline(in, line, '\t');
  gf.start = atoi(line.c_str());

  std::getline(in, line, '\t');
  gf.end = atoi(line.c_str());

  std::getline(in, line, '\t');
  gf.score = atoi(line.c_str());

  std::getline(in, line, '\t');
  gf.strand = line[0];

  std::getline(in, line, '\t');
  gf.phase = line[0];

  std::getline(in, line);
  using tokenizer = boost::tokenizer<boost::char_separator<char>>;
  boost::char_separator<char> sep(";");
  tokenizer tokens(line, sep);

  for (auto tokIt : tokens) {
    // Currently, we'll handle the following separators
    // '\s+'
    // '\s*=\s*'
    tokIt = tokIt.substr(tokIt.find_first_not_of(' '));
    auto kvsepStart = tokIt.find('=');

    // If we reached the end of the key, value token, then the string must have
    // been separated by some set of spaces, and NO '='.  If this is the case,
    // find the 'spaces' so that we can split on it.
    if (kvsepStart == tokIt.npos) {
      kvsepStart = tokIt.find(' ');
    }

    auto key = tokIt.substr(0, kvsepStart);
    key = key.substr(0, key.find(' '));

    auto kvsepStop =
        1 + kvsepStart + tokIt.substr(kvsepStart + 1).find_first_not_of(' ');
    auto val =
        (tokIt[kvsepStop] == '"')
            ? tokIt.substr(kvsepStop + 1, (tokIt.length() - (kvsepStop + 2)))
            : tokIt.substr(kvsepStop, (tokIt.length() - (kvsepStop + 1)));
    gf.sattr.parseAttribute(key, val);
  }

  return in;
}
