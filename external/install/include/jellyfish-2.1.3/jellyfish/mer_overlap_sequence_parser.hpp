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

#ifndef __JELLYFISH_MER_OVELAP_SEQUENCE_PARSER_H_
#define __JELLYFISH_MER_OVELAP_SEQUENCE_PARSER_H_

#include <stdint.h>

#include <memory>

#include <jellyfish/err.hpp>
#include <jellyfish/cooperative_pool2.hpp>
#include <jellyfish/cpp_array.hpp>

namespace jellyfish {

struct sequence_ptr {
  char* start;
  char* end;
};

template<typename StreamIterator>
class mer_overlap_sequence_parser : public jellyfish::cooperative_pool2<mer_overlap_sequence_parser<StreamIterator>, sequence_ptr> {
  typedef jellyfish::cooperative_pool2<mer_overlap_sequence_parser<StreamIterator>, sequence_ptr> super;
  enum file_type { DONE_TYPE, FASTA_TYPE, FASTQ_TYPE };
  typedef std::unique_ptr<std::istream> stream_type;

  struct stream_status {
    char*       seam;
    size_t      seq_len;
    bool        have_seam;
    file_type   type;
    stream_type stream;

    stream_status() : seam(0), seq_len(0), have_seam(false), type(DONE_TYPE) { }
  };

  uint16_t                       mer_len_;
  size_t                         buf_size_;
  char*                          buffer;
  char*                          seam_buffer;
  locks::pthread::mutex          streams_mutex;
  char*                          data;
  cpp_array<stream_status>       streams_;
  StreamIterator&                streams_iterator_;

public:
  /// Max_producers is the maximum number of concurrent threads than
  /// can produce data simultaneously. Size is the number of buffer to
  /// keep around. It should be larger than the number of thread
  /// expected to read from this class. buf_size is the size of each
  /// buffer. A StreamIterator is expected to have a next() method,
  /// which is thread safe, and which returns (move) a
  /// std::unique<std::istream> object.
  mer_overlap_sequence_parser(uint16_t mer_len, uint32_t max_producers, uint32_t size, size_t buf_size,
                              StreamIterator& streams) :
    super(max_producers, size),
    mer_len_(mer_len),
    buf_size_(buf_size),
    buffer(new char[size * buf_size]),
    seam_buffer(new char[max_producers * (mer_len - 1)]),
    streams_(max_producers),
    streams_iterator_(streams)
  {
    for(sequence_ptr* it = super::element_begin(); it != super::element_end(); ++it)
      it->start = it->end = buffer + (it - super::element_begin()) * buf_size;
    for(uint32_t i = 0; i < max_producers; ++i) {
      streams_.init(i);
      streams_[i].seam = seam_buffer + i * (mer_len - 1);
      open_next_file(streams_[i]);
    }
  }

  ~mer_overlap_sequence_parser() {
    delete [] buffer;
    delete [] seam_buffer;
  }

  //  file_type get_type() const { return type; }

  inline bool produce(uint32_t i, sequence_ptr& buff) {
    stream_status& st = streams_[i];

    switch(st.type) {
    case FASTA_TYPE:
      read_fasta(st, buff);
      break;
    case FASTQ_TYPE:
      read_fastq(st, buff);
      break;
    case DONE_TYPE:
      return true;
    }

    if(st.stream->good())
      return false;

    // Reach the end of file, close current and try to open the next one
    st.have_seam = false;
    open_next_file(st);
    return false;
  }

protected:
  bool open_next_file(stream_status& st) {
    // The stream must be released, with .reset(), before calling
    // .next() on the streams_iterator_, to ensure that the
    // streams_iterator_ noticed that we closed that stream before
    // requesting a new one.
    st.stream.reset();
    st.stream = streams_iterator_.next();
    if(!st.stream) {
      st.type = DONE_TYPE;
      return false;
    }

    switch(st.stream->peek()) {
    case EOF: return open_next_file(st);
    case '>':
      st.type = FASTA_TYPE;
      ignore_line(*st.stream); // Pass header
      break;
    case '@':
      st.type = FASTQ_TYPE;
      ignore_line(*st.stream); // Pass header
      break;
    default:
      eraise(std::runtime_error) << "Unsupported format"; // Better error management
    }
    return true;
  }

  void read_fasta(stream_status& st, sequence_ptr& buff) {
    size_t read = 0;
    if(st.have_seam) {
      memcpy(buff.start, st.seam, mer_len_ - 1);
      read = mer_len_ - 1;
    }

    // Here, the current stream is assumed to always point to some
    // sequence (or EOF). Never at header.
    while(st.stream->good() && read < buf_size_ - mer_len_ - 1) {
      read += read_sequence(*st.stream, read, buff.start, '>');
      if(st.stream->peek() == '>') {
        *(buff.start + read++) = 'N'; // Add N between reads
        ignore_line(*st.stream); // Skip to next sequence (skip headers, quals, ...)
      }
    }
    buff.end = buff.start + read;

    st.have_seam = read >= (size_t)(mer_len_ - 1);
    if(st.have_seam)
      memcpy(st.seam, buff.end - mer_len_ + 1, mer_len_ - 1);
  }

  void read_fastq(stream_status& st, sequence_ptr& buff) {
    size_t read = 0;
    if(st.have_seam) {
      memcpy(buff.start, st.seam, mer_len_ - 1);
      read = mer_len_ - 1;
    }

    // Here, the st.stream is assumed to always point to some
    // sequence (or EOF). Never at header.
    while(st.stream->good() && read < buf_size_ - mer_len_ - 1) {
      size_t nread  = read_sequence(*st.stream, read, buff.start, '+');
      read         += nread;
      st.seq_len   += nread;
      if(st.stream->peek() == '+') {
        skip_quals(*st.stream, st.seq_len);
        if(st.stream->good()) {
          *(buff.start + read++) = 'N'; // Add N between reads
          ignore_line(*st.stream); // Skip sequence header
        }
        st.seq_len = 0;
      }
    }
    buff.end = buff.start + read;

    st.have_seam = read >= (size_t)(mer_len_ - 1);
    if(st.have_seam)
      memcpy(st.seam, buff.end - mer_len_ + 1, mer_len_ - 1);
  }

  size_t read_sequence(std::istream& is, const size_t read, char* const start, const char stop) {
    size_t nread = read;
    while(is && nread < buf_size_ - 1 && is.peek() != stop) {
      // Skip new lines -> get below does like them
      skip_newlines(is);
      is.get(start + nread, buf_size_ - nread);
      nread += is.gcount();
      skip_newlines(is);
    }
    return nread - read;
  }

  inline void ignore_line(std::istream& is) {
    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  inline void skip_newlines(std::istream& is) {
    while(is.peek() == '\n')
      is.get();
  }

  // Skip quals header and qual values (read_len) of them.
  void skip_quals(std::istream& is, size_t read_len) {
    ignore_line(is);
    size_t quals = 0;
    while(is.good() && quals < read_len) {
      skip_newlines(is);
      is.ignore(read_len - quals + 1, '\n');
      quals += is.gcount();
      if(is)
        ++read_len;
    }
    skip_newlines(is);
    if(quals == read_len && (is.peek() == '@' || is.peek() == EOF))
      return;

    eraise(std::runtime_error) << "Invalid fastq sequence";
  }

  char peek(std::istream& is) { return is.peek(); }
};
}

#endif /* __JELLYFISH_MER_OVELAP_SEQUENCE_PARSER_H_ */
