#ifndef __JELLYFISH_WHOLE_SEQUENCE_PARSER_HPP__
#define __JELLYFISH_WHOLE_SEQUENCE_PARSER_HPP__

#include <string>
#include <memory>

#include <jellyfish/err.hpp>
#include <jellyfish/cooperative_pool2.hpp>
#include <jellyfish/cpp_array.hpp>

namespace jellyfish {
struct header_sequence_qual {
  std::string header;
  std::string seq;
  std::string qual;
};
struct sequence_list {
  size_t nb_filled;
  std::vector<header_sequence_qual> data;
};

template<typename StreamIterator>
class whole_sequence_parser : public jellyfish::cooperative_pool2<whole_sequence_parser<StreamIterator>, sequence_list> {
  typedef jellyfish::cooperative_pool2<whole_sequence_parser<StreamIterator>, sequence_list> super;
  typedef std::unique_ptr<std::istream> stream_type;
  enum file_type { DONE_TYPE, FASTA_TYPE, FASTQ_TYPE };

  struct stream_status {
    file_type   type;
    std::string buffer;
    stream_type stream;
    stream_status() : type(DONE_TYPE) { }
  };
  cpp_array<stream_status> streams_;
  StreamIterator&          streams_iterator_;

public:
  /// Size is the number of buffers to keep around. It should be
  /// larger than the number of thread expected to read from this
  /// class. nb_sequences is the number of sequences to read into a
  /// buffer. 'begin' and 'end' are iterators to a range of istream.
  whole_sequence_parser(uint32_t size, uint32_t nb_sequences,
                        uint32_t max_producers, StreamIterator& streams) :
    super(max_producers, size),
    streams_(max_producers),
    streams_iterator_(streams)
  {
    for(auto it = super::element_begin(); it != super::element_end(); ++it) {
      it->nb_filled = 0;
      it->data.resize(nb_sequences);
    }
    for(uint32_t i = 0; i < max_producers; ++i) {
      streams_.init(i);
      open_next_file(streams_[i]);
    }
  }

  inline bool produce(uint32_t i, sequence_list& buff) {
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
    open_next_file(st);
    return false;
  }

protected:
  void open_next_file(stream_status& st) {
    st.stream.reset();
    st.stream = streams_iterator_.next();
    if(!st.stream) {
      st.type = DONE_TYPE;
      return;
    }

    // Update the type of the current file and move past first header
    // to beginning of sequence.
    switch(st.stream->peek()) {
    case EOF: return open_next_file(st);
    case '>':
      st.type = FASTA_TYPE;
      break;
    case '@':
      st.type = FASTQ_TYPE;
      break;
    default:
      eraise(std::runtime_error) << "Unsupported format"; // Better error management
    }
  }

  void read_fasta(stream_status& st, sequence_list& buff) {
    size_t&      nb_filled = buff.nb_filled;
    const size_t data_size = buff.data.size();

    for(nb_filled = 0; nb_filled < data_size && st.stream->peek() != EOF; ++nb_filled) {
      header_sequence_qual& fill_buff = buff.data[nb_filled];
      st.stream->get(); // Skip '>'
      std::getline(*st.stream, fill_buff.header);
      fill_buff.seq.clear();
      while(st.stream->peek() != '>' && st.stream->peek() != EOF) {
        std::getline(*st.stream, st.buffer); // Wish there was an easy way to combine the
        fill_buff.seq.append(st.buffer);             // two lines avoiding copying
      }
    }
  }

  void read_fastq(stream_status& st, sequence_list& buff) {
    size_t&      nb_filled = buff.nb_filled;
    const size_t data_size = buff.data.size();

    for(nb_filled = 0; nb_filled < data_size && st.stream->peek() != EOF; ++nb_filled) {
      header_sequence_qual& fill_buff = buff.data[nb_filled];
      st.stream->get(); // Skip '@'
      std::getline(*st.stream, fill_buff.header);
      fill_buff.seq.clear();
      while(st.stream->peek() != '+' && st.stream->peek() != EOF) {
        std::getline(*st.stream, st.buffer); // Wish there was an easy way to combine the
        fill_buff.seq.append(st.buffer);             // two lines avoiding copying
      }
      if(!st.stream->good())
        eraise(std::runtime_error) << "Truncated fastq file";
      st.stream->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      fill_buff.qual.clear();
      while(fill_buff.qual.size() < fill_buff.seq.size() && st.stream->good()) {
        std::getline(*st.stream, st.buffer);
        fill_buff.qual.append(st.buffer);
      }
      if(fill_buff.qual.size() != fill_buff.seq.size())
        eraise(std::runtime_error) << "Invalid fastq file: wrong number of quals";
      if(st.stream->peek() != EOF && st.stream->peek() != '@')
        eraise(std::runtime_error) << "Invalid fastq file: header missing";
    }
  }
};
} // namespace jellyfish

#endif /* __JELLYFISH_WHOLE_SEQUENCE_PARSER_HPP__ */
