#include "salmon/internal/util/InputFileUtils.hpp"

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <spdlog/fmt/fmt.h>

#include "salmon/internal/io/AlignmentIO.hpp"
#include "salmon/internal/util/FmtCompat.hpp"

namespace salmon::utils {

bool peekBAMIsPaired(const boost::filesystem::path& file) {
  namespace bfs = boost::filesystem;
  std::string readMode = "r";

  if (bfs::is_regular_file(file)) {
    if (bfs::is_empty(file)) {
      fmt::MemoryWriter errstr;
      errstr << "file [" << file.string()
             << "] appears to be empty "
                "(i.e. it has size 0).  This is likely an error. "
                "Please re-run salmon with a corrected input file.\n\n";
      throw std::invalid_argument(errstr.str());
      return false;
    }
  }
  if (file.extension() == ".bam") {
    readMode = "rb";
  }

  auto* fp = scram_open(file.c_str(), readMode.c_str());

  // If we couldn't open the file, then report this and exit.
  if (fp == NULL) {
    fmt::MemoryWriter errstr;
    errstr << "ERROR: Failed to open file " << file.string() << ", exiting!\n";
    throw std::invalid_argument(errstr.str());
    return false;
  }

  bam_seq_t* read = nullptr;
  read = salmon::io::bam_init();

  bool didRead = (scram_get_seq(fp, &read) >= 0);
  bool isPaired{false};

  if (didRead) {
    isPaired = bam_flag(read) & BAM_FPAIRED;
  } else {
    fmt::MemoryWriter errstr;
    errstr << "ERROR: Failed to read alignment from " << file.string()
           << ", exiting!\n";
    salmon::io::bam_destroy(read);
    throw std::invalid_argument(errstr.str());
    return false;
  }

  scram_close(fp);
  salmon::io::bam_destroy(read);
  return isPaired;
}

size_t numberOfReadsInFastaFile(const std::string& fname) {
  constexpr size_t bufferSize = 16184;
  char buffer[bufferSize];
  std::ifstream ifile(fname, std::ifstream::in);
  ifile.rdbuf()->pubsetbuf(buffer, bufferSize);

  size_t numReads = 0;
  std::string s;
  while (ifile >> s) {
    if (s.front() == '>') {
      ++numReads;
    }
  }

  ifile.close();

  return numReads;
}

bool readKmerOrder(const std::string& fname, std::vector<uint64_t>& kmers) {
  std::ifstream mlist(fname, std::ios::in | std::ios::binary);
  // Get the number of kmers from file
  size_t numKmers{0};
  mlist.read(reinterpret_cast<char*>(&numKmers), sizeof(size_t));

  // Resize the array that will hold the sorted kmers
  kmers.resize(numKmers, 0);
  mlist.read(reinterpret_cast<char*>(&kmers[0]),
             sizeof(uint64_t) * kmers.size());

  mlist.close();

  return true;
}

} // namespace salmon::utils
