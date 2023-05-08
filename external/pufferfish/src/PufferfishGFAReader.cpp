#include "PufferfishGFAReader.hpp"
#include "CanonicalKmer.hpp"
#include "cereal/archives/binary.hpp"
#include "xxhash.h"
#include "Kmer.hpp"
#include "string_view.hpp"
#include <chrono>
#include <algorithm>
#include <string>
#include <bitset>


namespace kmers = combinelib::kmers;

namespace pufferfish {

enum class Direction : bool { PREPEND, APPEND } ;


namespace gfa_reader {
  namespace detail {
#define R -1
#define I -2
#define O -3
#define A 3
#define C 0
#define G 1
#define T 2
    static constexpr int shift_table[256] = {
      O, O, O, O, O, O, O, O, O, O, I, O, O, O, O, O, O, O, O, O, O, O, O, O,
      O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, R, O, O,
      O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, A, R, C, R, O, O, G,
      R, O, O, R, O, R, R, O, O, O, R, R, T, O, R, R, R, R, O, O, O, O, O, O,
      O, A, R, C, R, O, O, G, R, O, O, R, O, R, R, O, O, O, R, R, T, O, R, R,
      R, R, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
      O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
      O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
      O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
      O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
      O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O};
#undef R
#undef I
#undef O
#undef A
#undef C
#undef G
#undef T
  }
}

uint8_t encodeEdge(char c, Direction dir){
  //std::map<char,uint8_t> shift_table = {{'A',3}, {'T',2}, {'G',1}, {'C',0}};
  uint8_t val{1} ;
  if(dir == Direction::APPEND)
    return (val << gfa_reader::detail::shift_table[static_cast<uint8_t>(c)]) ;
  else
    return (val << (gfa_reader::detail::shift_table[static_cast<uint8_t>(c)]+4));
}



std::vector<std::pair<uint64_t, bool>>
GFAReader::explode(const stx::string_view str, const char& ch) {
  std::string next;
  std::vector<std::pair<uint64_t, bool>> result;
  // For each character in the string
  for (auto it = str.begin(); it != str.end(); it++) {
    // If we've hit the terminal character
    if (*it == '+' or *it == '-') {
      bool orientation = true;
      // If we have some characters accumulated
      // Add them to the result vector
      if (!next.empty()) {
        if (*it == '-') {
          orientation = false;
        }
        result.emplace_back(std::stoll(next), orientation);
        next.clear();
      }
    } else if (*it != ch) {
      // Accumulate the next character into the sequence
      next += *it;
    }
  }
  if (!next.empty())
    result.emplace_back(std::stoll(next),
                        true); // this case shouldn't even happen
  return result;
}

bool GFAReader::is_number(const std::string& s) {
  return !s.empty() && std::find_if(s.begin(), s.end(), [](char c) {
                         return !std::isdigit(c);
                       }) == s.end();
}

// Avoiding un-necessary stream creation + replacing strings with string view
// is a bit > than a 2x win!
// implementation from : https://marcoarena.wordpress.com/tag/string_view/
std::vector<stx::string_view> GFAReader::split(stx::string_view str,
                                               char delims) {
  std::vector<stx::string_view> ret;

  stx::string_view::size_type start = 0;
  auto pos = str.find_first_of(delims, start);
  while (pos != stx::string_view::npos) {
    if (pos != start) {
      ret.push_back(str.substr(start, pos - start));
    }
    start = pos + 1;
    pos = str.find_first_of(delims, start);
  }
  if (start < str.length()) {
    ret.push_back(str.substr(start, str.length() - start));
  }
  return ret;
}

GFAReader::GFAReader(const char* gfaFileName, size_t input_k, bool buildEdgeVec, std::shared_ptr<spdlog::logger> logger) {
  logger_ = logger;
  filename_ = std::string(gfaFileName);
  logger_->info("Reading GFA file {}", gfaFileName);
  file.reset(new zstr::ifstream(gfaFileName));
  k = input_k;
  buildEdgeVec_ = buildEdgeVec;
}

/*
void GFAReader::scanContigLengths() {
  std::string ln;
  std::string tag, id, value;
  size_t contig_cnt{0};
  size_t contig_len{0};
  size_t ref_cnt{0};
  while(std::getline(*file, ln)) {
    char firstC = ln[0];
    if (firstC != 'S') continue;
    stx::string_view lnview(ln);
    std::vector<stx::string_view> splited = split(lnview, '\t');
    if (is_number(id)) {
      contig_len += splited[1].length();
      tag = splited[0].to_string();
      id = splited[1].to_string();
      value = splited[2].to_string();
      contigid2seq[id] = value;
      contig_cnt++;
    }
  }
}
*/

size_t GFAReader::fillContigInfoMap_() {
  std::string ln;
  //std::string tag, id, value;
  size_t contig_ctr{0};
  //size_t contig_len{0};
  while (std::getline(*file, ln)) {
    char firstC = ln[0];
    if (firstC != 'S')
      continue;
    stx::string_view lnview(ln);
    std::vector<stx::string_view> splited = split(lnview, '\t');
    try {
      auto nid = std::stoull(splited[1].to_string());
      (void)nid;
      auto clen = splited[2].length();
      contigid2seq[nid] = {contig_ctr, 0, static_cast<uint32_t>(clen)};
      ++contig_ctr;
      // contigid2seq[id] = value;
      // contig_cnt++;
    } catch (std::exception& e) {
      // not numeric contig
    }
  }
  size_t total_len{0};
  uint64_t uId{0} ;
  for (auto& kv : contigid2seq) {
    kv.second.fileOrder = uId ;
    kv.second.offset = total_len;
    total_len += kv.second.length;
    ++uId ;
  }
  return total_len;
}

void GFAReader::encodeSeq(compact::vector<uint64_t, 2>& seqVec, size_t offset,
                          stx::string_view str) {
  for (size_t i = 0; i < str.length(); ++i) {
    auto c = kmers::codeForChar(str[i]);
    seqVec[offset + i] = c;
  }
}

compact::vector<uint64_t, 2>& GFAReader::getContigSeqVec() { return seqVec_; }
compact::vector<uint64_t, 8>& GFAReader::getEdgeVec() { return edgeVec_; }


void GFAReader::parseFile() {
  auto startTime = std::chrono::system_clock::now();
  size_t total_len = fillContigInfoMap_();
  auto now = std::chrono::system_clock::now();
  logger_->info("filling the contig map tool {:n} seconds.", std::chrono::duration_cast<std::chrono::seconds>(now - startTime).count());
  file.reset(new zstr::ifstream(filename_));
  logger_->info("total contig length = {:n} ", total_len);
  logger_->info("packing contigs into contig vector");
  seqVec_.resize(total_len);

  std::string ln;
  std::string tag, id, value;
  size_t contig_cnt{0};
  size_t ref_cnt{0};

  k = k + 1 ;
  CanonicalKmer::k(k) ;

  // start and end kmer-hash over the contigs
  // might get deprecated later
  uint64_t maxnid{0} ;
  while (std::getline(*file, ln)) {
    char firstC = ln[0];
    if (firstC != 'S' and firstC != 'P')
      continue;
    stx::string_view lnview(ln);
    std::vector<stx::string_view> splited = split(lnview, '\t');
    tag = splited[0].to_string();
    id = splited[1].to_string();
    // value = splited[2].to_string();
    // A segment line
    if (tag == "S") {
      try {
        uint64_t nid = std::stoll(id);
        if(nid > maxnid)
          maxnid = nid ;
        encodeSeq(seqVec_, contigid2seq[nid].offset, splited[2]);
        // contigid2seq[nid] = value;
      } catch (std::exception& e) {
        // not a numeric contig id
      }
      contig_cnt++;
    }


    // A path line
    if (tag == "P") {
      auto pvalue = splited[2];
      std::vector<std::pair<uint64_t, bool>> contigVec = explode(pvalue, ',');
      //go over the contigVec
      //uint64_t kn{0} ;
      //uint64_t knn{0} ;

      // parse value and add all conitgs to contigVec

      path[ref_cnt] = contigVec;

      uint32_t refLength{0};
      bool firstContig{true};
      for (auto& ctig : contigVec) {
        int32_t l = contigid2seq[ctig.first].length - (firstContig ? 0 : (k-1));
        refLength += l;
        firstContig = false;
      }

      refLengths.push_back(refLength);
      refMap.push_back(id);
      ref_cnt++;
      // refMap[ref_cnt] = id;
      // refIDs[id] = ref_cnt++;
    }
  }

  //Initialize edgeVec_
  //bad way, have to re-think
  if (buildEdgeVec_) {
    edgeVec_.resize(contig_cnt);
    for(uint64_t i=0; i<contig_cnt; i++) edgeVec_[i]=0;

    for(auto const& ent: path){
      const std::vector<std::pair<uint64_t, bool>>& contigs = ent.second;

      for(size_t i = 0 ; i < contigs.size() - 1 ; i++){
        auto cid = contigs[i].first ;
        bool ore = contigs[i].second ;
        size_t forder = contigid2seq[cid].fileOrder ;
        auto nextcid = contigs[i+1].first ;
        bool nextore = contigs[i+1].second ;

        bool nextForder = contigid2seq[nextcid].fileOrder ;
        // a+,b+ end kmer of a , start kmer of b
        // a+,b- end kmer of a , rc(end kmer of b)
        // a-,b+ rc(start kmer of a) , start kmer of b
        // a-,b- rc(start kmer of a) , rc(end kmer of b)
        //  1. `a+,(*)` we need to append a nucl to the `end-kmer`
        //  2. `a-,(*)` we need to prepend rc(nucl) to the `start-kmer`
        //  3. `(*),a+` we need to prepend a nucl to `start-kmer`
        //  4. `(*),a-` we need to append a rc(nucl) to `end-kmer`

        CanonicalKmer lastKmerInContig;
        CanonicalKmer firstKmerInNextContig;
        Direction contigDirection;
        Direction nextContigDirection;
        // If a is in the forward orientation, the last k-mer comes from the end, otherwise it is the reverse complement of the first k-mer
        if (ore) {
          lastKmerInContig.fromNum(seqVec_.get_int(2 * (contigid2seq[cid].offset + contigid2seq[cid].length - k), 2 * k));
          contigDirection = Direction::APPEND;
        } else {
          lastKmerInContig.fromNum(seqVec_.get_int(2 * contigid2seq[cid].offset, 2*k));
          lastKmerInContig.swap();
          contigDirection = Direction::PREPEND;
        }

        // If a is in the forward orientation, the first k-mer comes from the beginning, otherwise it is the reverse complement of the last k-mer
        if (nextore) {
          firstKmerInNextContig.fromNum(seqVec_.get_int(2 * contigid2seq[nextcid].offset, 2*k));
          nextContigDirection = Direction::PREPEND;
        } else {
          firstKmerInNextContig.fromNum(seqVec_.get_int(2 * (contigid2seq[nextcid].offset + contigid2seq[nextcid].length - k), 2 * k));
          firstKmerInNextContig.swap();
          nextContigDirection = Direction::APPEND;
        }

        // The character to append / prepend to contig to get to next contig
        const char contigChar = firstKmerInNextContig.to_str()[k-1];
        // The character to prepend / append to next contig to get to contig
        const char nextContigChar = lastKmerInContig.to_str()[0];

        edgeVec_[forder] = edgeVec_[forder] | encodeEdge(contigChar, contigDirection);
        edgeVec_[nextForder] = edgeVec_[nextForder] | encodeEdge(nextContigChar, nextContigDirection);
      }
    }
  }
  k = k - 1;
  logger_->info("Total # of Contigs : {:n}", contig_cnt);
  logger_->info("Total # of numerical Contigs : {:n}", contigid2seq.size());
}

// spp::sparse_hash_map<uint64_t, std::string>& GFAReader::getContigNameMap() {
spp::sparse_hash_map<uint64_t, pufferfish::util::PackedContigInfo>&
GFAReader::getContigNameMap() {
  return contigid2seq;
}
spp::sparse_hash_map<std::string, std::string>& GFAReader::getContigIDMap() {
  return seq2contigid;
}

/*
spp::sparse_hash_map<uint32_t, std::string>& GFAReader::getRefIDs() {
  return refMap;
}
*/
std::vector<std::string>& GFAReader::getRefIDs() { return refMap; }

std::vector<uint32_t>& GFAReader::getRefLengths() { return refLengths;}

std::map<std::pair<std::string, bool>, bool, pufferfish::util::cmpByPair>&
GFAReader::getPathStart() {
  return pathStart;
}
std::map<std::pair<std::string, bool>, bool, pufferfish::util::cmpByPair>&
GFAReader::getPathEnd() {
  return pathEnd;
}

std::vector<std::pair<std::string, std::string>>& GFAReader::getNewSegments() {
  return newSegments;
}

pufg::Graph& GFAReader::getSemiCG() { return semiCG; }

/*
void GFAReader::writeFile(std::string fileName){
    std::ofstream gfa_file(fileName) ;
    for(auto& cseq : contigid2seq){
        gfa_file << "S" << "\t" << cseq.first <<"\t" << cseq.second << "\n" ;
    }
    for(auto& p : path){
        auto tid = p.first ;
        gfa_file << "P" << "\t" << tid << "\t"  ;
        auto vec = p.second ;
        for(size_t i = 0 ; i < vec.size()-1 ; i++){
            gfa_file << vec[i].first << ((vec[i].second)?"+":"-") << "," ;
        }
        gfa_file << vec[vec.size()-1].first <<
((vec[vec.size()-1].second)?"+":"-") << "\t*\n";

    }
}
*/

void GFAReader::mapContig2Pos() {
  uint64_t pos = 0;
  uint64_t accumPos;
  uint64_t currContigLength = 0;
  uint64_t total_output_lines = 0;
  for (auto const& ent : path) {
    const uint64_t& tr = ent.first;
    const std::vector<std::pair<uint64_t, bool>>& contigs = ent.second;
    accumPos = 0;
    for (size_t i = 0; i < contigs.size(); i++) {
      if (contig2pos.find(contigs[i].first) == contig2pos.end()) {
        contig2pos[contigs[i].first] = {};
        total_output_lines += 1;
      }
      if (contigid2seq.find(contigs[i].first) == contigid2seq.end()) {
        logger_->info("{}", contigs[i].first);
      }
      pos = accumPos;
      currContigLength = contigid2seq[contigs[i].first].length;
      accumPos += currContigLength - k;
      (contig2pos[contigs[i].first])
          .push_back(pufferfish::util::Position(tr, pos, contigs[i].second));
    }
  }
  logger_->info("\nTotal # of segments we have position for : {:n}", total_output_lines);
}

void GFAReader::clearContigTable() {
  refMap.clear();
  refLengths.clear();
  contig2pos.clear();
}

// Note : We assume that odir is the name of a valid (i.e., existing) directory.
void GFAReader::serializeContigTable(const std::string& odir) {
  std::string ofile = odir + "/ctable.bin";
  std::string eqfile = odir + "/eqtable.bin";
  std::string rlfile = odir + "/reflengths.bin";
  std::ofstream ct(ofile);
  std::ofstream et(eqfile);
  std::ofstream rl(rlfile);
  cereal::BinaryOutputArchive ar(ct);
  cereal::BinaryOutputArchive eqAr(et);
  cereal::BinaryOutputArchive rlAr(rl);
  {
    // Write out the reference lengths
    rlAr(refLengths);

    // We want to iterate over the contigs in precisely the
    // order they appear in the contig array (i.e., the iterator
    // order of contigid2seq).
    std::vector<std::string> refNames;
    refNames.reserve(refMap.size());
    for (size_t i = 0; i < refMap.size(); ++i) {
      refNames.push_back(refMap[i]);
    }
    ar(refNames);

    class VecHasher {
    public:
      size_t operator()(const std::vector<uint32_t>& vec) const {
        return XXH64(const_cast<std::vector<uint32_t>&>(vec).data(),
                     vec.size() * sizeof(decltype(vec.front())), 0);
      }
    };

    // Write out contig offsets and lengths
    /*
    {
      std::vector<pufferfish::util::ContigPosInfo> cpi;
      cpi.reserve(contigid2seq.size());
      for (auto& kv : contigid2seq) {
        cpi.push_back({kv.second.offset, kv.second.length});
      }
      ar(cpi);
    }
    */

    spp::sparse_hash_map<std::vector<uint32_t>, uint32_t, VecHasher> eqMap;
    std::vector<uint32_t> eqIDs;
    //std::vector<std::vector<pufferfish::util::Position>> cpos;

    // Compute sizes to reserve
    size_t contigVecSize{0};
    size_t contigOffsetSize{1};
    for (auto& kv : contigid2seq) {
      contigOffsetSize++;
      contigVecSize += contig2pos[kv.first].size();
    }

    logger_->info("total contig vec entries {:n}", contigVecSize);
    std::vector<pufferfish::util::Position> cpos;
    cpos.reserve(contigVecSize);

    // We need the +1 here because we store the last entry that is 1 greater than the last offset
    // so we must be able to represent of number of size contigVecSize+1, not contigVecSize.
    size_t w = std::ceil(std::log2(contigVecSize+1));
    logger_->info("bits per offset entry {:n}", w);
    compact::vector<uint64_t> cpos_offsets(w, contigOffsetSize);

    //std::vector<uint64_t> cpos_offsets;
    //cpos_offsets.reserve(contigOffsetSize);
    //cpos_offsets.push_back(0);
    size_t idx{0};
    cpos_offsets[0] = 0;
    for (auto& kv : contigid2seq) {
      ++idx;
      //cpos.push_back(contig2pos[kv.first]);
      auto& b = contig2pos[kv.first];
      cpos_offsets[idx] = cpos_offsets[idx-1] + b.size();
      //cpos_offsets.push_back(cpos_offsets.back() + b.size());
      cpos.insert(cpos.end(), std::make_move_iterator(b.begin()), std::make_move_iterator(b.end()));
      std::vector<uint32_t> tlist;
      for (auto& p : contig2pos[kv.first]) {
        tlist.push_back(p.transcript_id());
      }
      std::sort(tlist.begin(), tlist.end());
      tlist.erase(std::unique(tlist.begin(), tlist.end()), tlist.end());
      size_t eqID = eqMap.size();
      if (eqMap.contains(tlist)) {
        eqID = eqMap[tlist];
      } else {
        eqMap[tlist] = eqID;
      }
      eqIDs.push_back(eqID);
      // ar(contig2pos[kv.first]);
      /*
      auto& ent = contig2pos[kv.first];
      ct << kv.first;
      for (auto const & ent2 : ent) {
        ct << '\t' << ent2.transcript_id << ',' << ent2.pos << ',' <<
      ent2.orien;
      }
      ct << '\n';
      */
    }
    logger_->info("there were {:n}  equivalence classes", eqMap.size());
    eqAr(eqIDs);
    eqIDs.clear();
    eqIDs.shrink_to_fit();
    std::vector<std::vector<uint32_t>> eqLabels;
    eqLabels.reserve(eqMap.size());
    for (auto& kv : eqMap) {
      eqLabels.push_back(kv.first);
    }
    std::sort(eqLabels.begin(), eqLabels.end(),
              [&](const std::vector<uint32_t>& l1,
                  const std::vector<uint32_t>& l2) -> bool {
                return eqMap[l1] < eqMap[l2];
              });
    eqAr(eqLabels);
    ar(cpos);
    //ar(cpos_offsets);
    {
      std::string fname = odir + pufferfish::util::CONTIG_OFFSETS;
      std::ofstream bfile(fname, std::ios::binary);
      cpos_offsets.serialize(bfile);
      bfile.close();
    }

  }
  /*
    ct << refIDs.size() << '\n';
    for (auto const & ent : refIDs) {
        ct << ent.first << '\t' << ent.second << '\n';
    }
    for (auto const & ent : contig2pos) {
        ct << ent.first;
        for (auto const & ent2 : ent.second) {
            ct << '\t' << ent2.transcript_id << ',' << ent2.pos << ',' <<
    ent2.orien;
        }
        ct << '\n';
    }
  */
}

void GFAReader::deserializeContigTable() {
  // TODO read the file in the same order as you've written it down.
}
/*
int main(int argc, char* argv[]) {
    //std::vector<std::string> gfa_file = {argv[1]};
    //std::vector<std::string> read_file = {argv[2]};
    std::cerr<<" Starting ... \n";
    GFAReader pf(argv[1], 30);
    pf.parseFile();

    pf.mapContig2Pos();

    std::ofstream ct("contig_table.tsv");
    for (auto const & ent : pf.contig2pos) {
        ct << ent.first;
        for (auto const & ent2 : ent.second) {
            ct << '\t' << ent2.transcript_id << ',' << ent2.pos << ',' <<
ent2.orien;
        }
        ct << '\n';
    }
    ct.close();
}
*/

} // namespace pufferfish
