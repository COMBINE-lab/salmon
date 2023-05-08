#ifndef OUR_GFA_READER_H
#define OUR_GFA_READER_H

#include "FatPufferGraph.hpp"
#include "Util.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/vector.hpp"
//#include "sdsl/int_vector.hpp"
#include "sparsepp/spp.h"
#include "spdlog/spdlog.h"
#include "compact_vector/compact_vector.hpp"

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

class GFAReader {
private:
  std::string filename_;
  std::unique_ptr<zstr::ifstream> file;
  size_t k;
  struct Contig {
    std::string seq;
    std::string id;
  };

  spp::sparse_hash_map<uint64_t, pufferfish::util::PackedContigInfo>
      contigid2seq; // map of contig_id to # of letters in contig (contig
                    // length)
  spp::sparse_hash_map<std::string, std::string> seq2contigid;
  // path maps each transcript_id to a pair of <contig_id, orientation>
  // orientation : +/true main, -/false reverse
  spp::sparse_hash_map<uint64_t, std::vector<std::pair<uint64_t, bool>>> path;
  spp::sparse_hash_map<uint64_t, uint32_t> refIDs;

  // spp::sparse_hash_map<uint64_t, std::string> refMap;
  std::vector<std::string> refMap;
  std::vector<uint32_t> refLengths;
  // maps each contig to a list of positions in different transcripts
  std::vector<std::pair<uint64_t, bool>> explode(const stx::string_view str,
                                                 const char& ch);

  std::map<std::pair<std::string, bool>, bool, pufferfish::util::cmpByPair> pathStart;
  std::map<std::pair<std::string, bool>, bool, pufferfish::util::cmpByPair> pathEnd;

  compact::vector<uint64_t, 2> seqVec_;

  //edge table
  //ATGC|ATGC = 8 bits
  compact::vector<uint64_t, 8> edgeVec_;
  //predecessor,stores the same
  //transcript in reverse order
  //improve walkability
  //sdsl::int_vector<8> edgeVec2_;

  std::vector<std::pair<std::string, std::string>> newSegments;
  pufg::Graph semiCG;

  size_t fillContigInfoMap_();
  bool is_number(const std::string& s);

  // Avoiding un-necessary stream creation + replacing strings with string view
  // is a bit > than a 2x win!
  // implementation from : https://marcoarena.wordpress.com/tag/string_view/
  std::vector<stx::string_view> split(stx::string_view str, char delims);

  bool buildEdgeVec_{false};
  std::shared_ptr<spdlog::logger> logger_{nullptr};

public:
  spp::sparse_hash_map<uint64_t, std::vector<pufferfish::util::Position>> contig2pos;

  GFAReader(const char* gfaFileName, size_t input_k,
            bool buildEdgeVEc, std::shared_ptr<spdlog::logger> logger);

  /*void encodeSeq(sdsl::int_vector<2>& seqVec, size_t offset,
                 stx::string_view str);
  */
 
  void encodeSeq(compact::vector<uint64_t, 2>& seqVec, size_t offset,
                 stx::string_view str);
  // spp::sparse_hash_map<uint64_t, std::string>& getContigNameMap();
  spp::sparse_hash_map<uint64_t, pufferfish::util::PackedContigInfo>& getContigNameMap();

  spp::sparse_hash_map<std::string, std::string>& getContigIDMap();
  // spp::sparse_hash_map<uint32_t, std::string>& getRefIDs();
  std::vector<std::string>& getRefIDs();
  std::vector<uint32_t>& getRefLengths();
  std::map<std::pair<std::string, bool>, bool, pufferfish::util::cmpByPair>& getPathStart();
  std::map<std::pair<std::string, bool>, bool, pufferfish::util::cmpByPair>& getPathEnd();
  std::vector<std::pair<std::string, std::string>>& getNewSegments();
  compact::vector<uint64_t, 2>& getContigSeqVec();
  compact::vector<uint64_t, 8>& getEdgeVec();
  compact::vector<uint64_t, 8>& getEdgeVec2();
  
  // spp::sparse_hash_map<std::string, std::vector< std::pair<std::string, bool>
  // > >& getPaths() {return path;}
  // spp::sparse_hash_map<std::string, std::vector< std::pair<std::string, bool>
  // > >& getPaths() {return path;}
  pufg::Graph& getSemiCG();

  void parseFile();
  void mapContig2Pos();
  void clearContigTable();
  void serializeContigTable(const std::string& odir);
  void deserializeContigTable();
  // void writeFile(std::string fileName);
};

} // end namespace pufferfish
#endif
