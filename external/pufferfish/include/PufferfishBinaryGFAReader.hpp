#ifndef OUR_GFA_READER_H
#define OUR_GFA_READER_H

#include "CLI/Timer.hpp"
#include "Util.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/vector.hpp"
//#include "sdsl/int_vector.hpp"
#include "sparsepp/spp.h"
#include "spdlog/spdlog.h"
#include "compact_vector/compact_vector.hpp"
#include "rank9sel.hpp"

#include "string_view.hpp"

#include "zstr/zstr.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

namespace pufferfish {

class BinaryGFAReader {
private:
  std::string filename_;
  std::unique_ptr<zstr::ifstream> file;
  size_t k;
  struct Contig {
    std::string seq;
    std::string id;
  };


  pufferfish::util::PackedContigInfoVec contigid2seq;
  
  // path maps each transcript_id to a pair of <contig_id, orientation>
  // orientation : +/true main, -/false reverse
  spp::sparse_hash_map<uint64_t, std::vector<std::pair<uint64_t, bool>>> path;
  spp::sparse_hash_map<uint64_t, uint32_t> refIDs;

  // spp::sparse_hash_map<uint64_t, std::string> refMap;
  std::vector<std::string> refMap;
  std::vector<uint32_t> refLengths;
  uint64_t maxRefLength{0};

  compact::vector<uint64_t, 2> seqVec_;
  compact::vector<uint64_t, 1> rankVec_;

  //edge table
  //ATGC|ATGC = 8 bits
  compact::vector<uint64_t, 8> edgeVec_;
  //predecessor,stores the same
  //transcript in reverse order
  //improve walkability
  //sdsl::int_vector<8> edgeVec2_;

  size_t fillContigInfoMap_();

  // Avoiding un-necessary stream creation + replacing strings with string view
  // is a bit > than a 2x win!
  // implementation from : https://marcoarena.wordpress.com/tag/string_view/
  std::vector<stx::string_view> split(stx::string_view str, char delims);

  bool buildEdgeVec_{false};
  bool buildEqClses_{false};
  std::shared_ptr<spdlog::logger> logger_{nullptr};
  std::unique_ptr<compact::vector<uint64_t>> cpos_offsets{nullptr};

public:
//  spp::sparse_hash_map<uint64_t, std::vector<pufferfish::util::Position>>
  std::vector<pufferfish::util::Position> contig2pos;

  BinaryGFAReader(const char* gfaFileName, size_t input_k,
            bool buildEqClses, bool buildEdgeVEc,
            std::shared_ptr<spdlog::logger> logger);

  /*void encodeSeq(sdsl::int_vector<2>& seqVec, size_t offset,
                 stx::string_view str);
  */ 
  void encodeSeq(compact::vector<uint64_t, 2>& seqVec, size_t offset,
                 stx::string_view str);

  //spp::sparse_hash_map<uint64_t, pufferfish::util::PackedContigInfo>& getContigNameMap();
  pufferfish::util::PackedContigInfoVec& getContigNameMap();

  std::vector<std::string>& getRefIDs();
  std::vector<uint32_t>& getRefLengths();
  compact::vector<uint64_t, 2>& getContigSeqVec();
  compact::vector<uint64_t, 1>& getRankVec();
  compact::vector<uint64_t, 8>& getEdgeVec();

  void parseFile();
  void mapContig2Pos();
  void clearContigTable();
  void serializeContigTable(const std::string& odir,
          const std::vector<std::pair<std::string, uint16_t>>& shortRefsNameLen,
          const std::vector<uint32_t>& refIdExtensions);

    void deserializeContigTable();
  // void writeFile(std::string fileName);
};

} // end namespace pufferfish
#endif
