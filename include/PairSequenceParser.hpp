#ifndef __PAIR_SEQUENCE_PARSER_HPP__
#define __PAIR_SEQUENCE_PARSER_HPP__

#include <string>
#include <memory>
#include <utility>
#include <vector>
#include <thread>
#include <mutex>
#include <fstream>

#include <jellyfish/err.hpp>
#include <jellyfish/cooperative_pool2.hpp>
#include <jellyfish/cpp_array.hpp>

struct header_sequence_qual {
  std::string header;
  std::string seq;
  std::string qual;
};
struct sequence_list {
  size_t nb_filled;
  std::vector<std::pair<header_sequence_qual, header_sequence_qual> > data;
};

template<typename PathIterator>
class pair_sequence_parser : public jellyfish::cooperative_pool2<pair_sequence_parser<PathIterator>, sequence_list> {
  typedef jellyfish::cooperative_pool2<pair_sequence_parser<PathIterator>, sequence_list> super;
  typedef std::unique_ptr<std::istream> stream_type;
  enum file_type { DONE_TYPE, FASTA_TYPE, FASTQ_TYPE, ERROR_TYPE };

  struct stream_status {
    file_type   type;
    std::string buffer;
    stream_type stream1;
    stream_type stream2;

    stream_status() : type(DONE_TYPE) { }
  };
  jellyfish::cpp_array<stream_status> streams_;
  PathIterator                        path_begin_, path_end_;
  std::mutex                          path_mutex_;

public:
  /// Size is the number of buffers to keep around. It should be
  /// larger than the number of thread expected to read from this
  /// class. nb_sequences is the number of sequences to read into a
  /// buffer. 'begin' and 'end' are iterators to a range of istream.
  pair_sequence_parser(uint32_t size, uint32_t nb_sequences,
                       uint32_t max_producers,
                       PathIterator path_begin, PathIterator path_end) :
    super(max_producers, size),
    streams_(max_producers),
    path_begin_(path_begin), path_end_(path_end)
  {
    for(auto it = super::element_begin(); it != super::element_end(); ++it) {
      it->nb_filled = 0;
      it->data.resize(nb_sequences);
    }
    for(uint32_t i = 0; i < max_producers; ++i) {
      streams_.init(i);
      open_next_files(streams_[i]);
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
    case ERROR_TYPE:
      return true;
    }

    if(st.stream1->good() && st.stream2->good())
      return false;

    // Reach the end of file, close current and try to open the next one
    open_next_files(st);
    return false;
  }

protected:
  file_type peek_file_type(std::istream& is) {
    switch(is.peek()) {
    case EOF: return DONE_TYPE;
    case '>': return FASTA_TYPE;
    case '@': return FASTQ_TYPE;
    default: return ERROR_TYPE;
    }
  }

  void open_next_files(stream_status& st) {
    st.stream1.reset();
    st.stream2.reset();
    const char *p1 = 0, *p2 = 0;
    {
      std::lock_guard<std::mutex> lck(path_mutex_);
      if(path_begin_ < path_end_) {
        p1 = *path_begin_;
        ++path_begin_;
      }
      if(path_begin_ < path_end_) {
        p2 = *path_begin_;
        ++path_begin_;
      }
    }

    if(!p1 || !p2) {
      st.type = DONE_TYPE;
      return;
    }
    st.stream1.reset(new std::ifstream(p1));
    st.stream2.reset(new std::ifstream(p2));
    if(!*st.stream1 || !*st.stream2) {
      st.type = DONE_TYPE;
      return;
    }

    // Update the type of the current file and move past first header
    // to beginning of sequence.
    file_type type1 = peek_file_type(*st.stream1);
    file_type type2 = peek_file_type(*st.stream2);
    if(type1 == DONE_TYPE || type2 == DONE_TYPE)
      return open_next_files(st);
    if(type1 != type2)
      throw std::runtime_error("Paired files are of different format");
    if(type1 == ERROR_TYPE || type2 == ERROR_TYPE)
      throw std::runtime_error("Unsupported format");
    st.type = type1;
  }

  void read_fasta_one_sequence(std::istream& is, std::string& tmp, header_sequence_qual& hsq) {
    is.get(); // Skip '>'
    std::getline(is, hsq.header);
    hsq.seq.clear();
    while(is.peek() != '>' && is.peek() != EOF) {
      std::getline(is, tmp); // Wish there was an easy way to combine the
      hsq.seq.append(tmp);             // two lines avoiding copying
    }
  }

  void read_fasta(stream_status& st, sequence_list& buff) {
    size_t&      nb_filled = buff.nb_filled;
    const size_t data_size = buff.data.size();

    for(nb_filled = 0; nb_filled < data_size && st.stream1->peek() != EOF && st.stream2->peek() != EOF; ++nb_filled) {
      read_fasta_one_sequence(*st.stream1, st.buffer, buff.data[nb_filled].first);
      read_fasta_one_sequence(*st.stream2, st.buffer, buff.data[nb_filled].second);
    }
  }

  void read_fastq_one_sequence(std::istream& is, std::string& tmp, header_sequence_qual& hsq) {
    is.get(); // Skip '@'
    std::getline(is, hsq.header);
    hsq.seq.clear();
    while(is.peek() != '+' && is.peek() != EOF) {
      std::getline(is, tmp); // Wish there was an easy way to combine the
      hsq.seq.append(tmp);             // two lines avoiding copying
    }
    if(!is.good())
      throw std::runtime_error("Truncated fastq file");
    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    hsq.qual.clear();
    while(hsq.qual.size() < hsq.seq.size() && is.good()) {
      std::getline(is, tmp);
      hsq.qual.append(tmp);
    }
    if(hsq.qual.size() != hsq.seq.size())
      throw std::runtime_error("Invalid fastq file: wrong number of quals");
    if(is.peek() != EOF && is.peek() != '@')
      throw std::runtime_error("Invalid fastq file: header missing");

  }

  void read_fastq(stream_status& st, sequence_list& buff) {
    size_t&      nb_filled = buff.nb_filled;
    const size_t data_size = buff.data.size();

    for(nb_filled = 0; nb_filled < data_size && st.stream1->peek() != EOF && st.stream2->peek() != EOF; ++nb_filled) {
      read_fastq_one_sequence(*st.stream1, st.buffer, buff.data[nb_filled].first);
      read_fastq_one_sequence(*st.stream2, st.buffer, buff.data[nb_filled].second);
    }
  }
};

#endif /* __PAIR_SEQUENCE_PARSER_HPP__ */
