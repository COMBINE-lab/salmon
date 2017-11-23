/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Sailfish.

    Sailfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sailfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sailfish.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/

#ifndef __READPRODUCER_HPP__
#define __READPRODUCER_HPP__

#include "StreamingSequenceParser.hpp"
#include "jellyfish/compacted_hash.hpp"
#include "jellyfish/dna_codes.hpp"
#include "jellyfish/mapped_file.hpp"
#include "jellyfish/mer_counting.hpp"
#include "jellyfish/misc.hpp"
#include "jellyfish/parse_dna.hpp"
#include "jellyfish/parse_read.hpp"
#include "jellyfish/sequence_parser.hpp"

template <typename Parser> class ReadProducer {};

template <> class ReadProducer<jellyfish::parse_read> {
public:
  ReadProducer(jellyfish::parse_read& parser)
      : s_(new ReadSeq), stream_(parser.new_thread()) {}
  ~ReadProducer() { delete s_; }

  bool nextRead(ReadSeq*& s) {
    if ((read_ = stream_.next_read())) {
      // we iterate over the entire read
      const char* start = read_->seq_s;
      const char* const end = read_->seq_e;
      uint32_t readLen = std::distance(start, end);

      s_->seq = const_cast<char*>(read_->seq_s);
      s_->len = readLen;
      s_->name = const_cast<char*>(read_->header);
      s_->nlen = read_->hlen;
      s = s_;
      return true;
    } else {
      return false;
    }
  }

  void finishedWithRead(ReadSeq*& s) { s = nullptr; }

private:
  ReadSeq* s_;
  jellyfish::parse_read::read_t* read_;  //{parser.new_thread()};
  jellyfish::parse_read::thread stream_; //{parser.new_thread()};
};

template <> class ReadProducer<StreamingReadParser> {
public:
  ReadProducer(StreamingReadParser& parser) : parser_(parser) {}

  bool nextRead(ReadSeq*& s) {
    if (parser_.nextRead(s)) {
      return true;
    } else {
      s = nullptr;
      return false;
    }
  }

  void finishedWithRead(ReadSeq*& s) { parser_.finishedWithRead(s); }

private:
  StreamingReadParser& parser_;
};

#endif // __READPRODUCER_HPP__
