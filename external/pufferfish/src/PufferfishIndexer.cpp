#include "FastxParser.hpp"
#include <cmath>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>
#include <sstream>
#include <bitset>
#include <chrono> 
#include <cerrno>
#include <cstring>
#include <memory>
#include <thread>
#include <cereal/archives/binary.hpp>
//#include <unistd.h>
#include "ghc/filesystem.hpp"

#include "ProgOpts.hpp"
#include "CanonicalKmer.hpp"
#include "PufferfishBinaryGFAReader.cpp"
#include "PufferFS.hpp"
#include "PufferfishIndex.hpp"
#include "ScopedTimer.hpp"
#include "Util.hpp"
#include "PufferfishConfig.hpp"
#include "cereal/archives/json.hpp"
#include "jellyfish/mer_dna.hpp"
#include "rank9b.hpp"
#include "spdlog/spdlog.h"
#include "Kmer.hpp" // currently requires k <= 32
#include "compact_vector/compact_vector.hpp"

namespace kmers = combinelib::kmers;

uint64_t swap_uint64(uint64_t val) {
  val = ((val << 8) & 0xFF00FF00FF00FF00ULL) |
        ((val >> 8) & 0x00FF00FF00FF00FFULL);
  val = ((val << 16) & 0xFFFF0000FFFF0000ULL) |
        ((val >> 16) & 0x0000FFFF0000FFFFULL);
  return (val << 32) | (val >> 32);
}

// adapted from :
// http://stackoverflow.com/questions/34875315/implementation-my-own-list-and-iterator-stl-c
class ContigKmerIterator {
public:
  using self_type = ContigKmerIterator;
  using value_type = uint64_t;
  using reference = value_type&;
  using pointer = value_type*;
  using iterator_category = std::forward_iterator_tag;
  using difference_type = int64_t;

  ContigKmerIterator() = delete;

  ContigKmerIterator(compact::vector<uint64_t, 2>* storage, compact::vector<uint64_t, 1>* rank,
                     uint8_t k, uint64_t startAt)
      : storage_(storage), rank_(rank), k_(k), curr_(startAt) {
    if (curr_ + k_ <= rank_->size()) {
      //nextValidPosition_();
      mer_.fromNum(storage_->get_int(2*curr_, 2*k_));
      // mer_.word__(0) = storage_->get_int(2 * curr_, 2 * k_);
    }
    // rcMer_ = mer_.get_reverse_complement();
  }

  bool advanceToValid() {
    nextValidPosition_();
    return (curr_ + k_ <= rank_->size());
  }

  ContigKmerIterator(const ContigKmerIterator& other) { 
    storage_ = other.storage_;
    rank_ = other.rank_;
    k_ = other.k_;
    curr_ = other.curr_;
    mer_ = other.mer_;
    word_ = other.word_;
  }

  ContigKmerIterator&
  operator=(ContigKmerIterator& other) { //}= default;
                                         // storage, 
    // uint8_t k, uint64_t startAt) :
    storage_ = other.storage_;
    rank_ = other.rank_;
    k_ = other.k_;
    curr_ = other.curr_;
    mer_ = other.mer_;
    // rcMer_ = other.rcMer_;
    word_ = other.word_;
    return *this;
  }

  ContigKmerIterator operator++() {
    ContigKmerIterator i = *this;
    advance_();
    return i;
  }

  ContigKmerIterator operator++(int) {
    advance_();
    return *this;
  }

  reference operator*() {
    // word_ = (mer_.word(0) < rcMer_.word(0)) ? mer_.word(0) : rcMer_.word(0);
    word_ = mer_.getCanonicalWord();
    return word_;
  }

  difference_type pos() { return curr_; }

  bool isCanonical(){
	  return mer_.fwWord() == mer_.getCanonicalWord() ;
  }

  bool isEndKmer() {
	  size_t endPos = curr_ + k_ - 1;
	  return ((*rank_)[endPos] == 1);
  }

  pointer operator->() {
    word_ = mer_.getCanonicalWord(); //(mer_.word(0) < rcMer_.word(0)) ?
                                     // mer_.word(0) : rcMer_.word(0);
    return &word_;
  }
  bool operator==(const self_type& rhs) { return curr_ == rhs.curr_; }

  bool operator!=(const self_type& rhs) { return curr_ != rhs.curr_; }

  bool operator<(const self_type& rhs) { return curr_ < rhs.curr_; }

  bool operator<=(const self_type& rhs) { return curr_ <= rhs.curr_; }

  inline bool isValid() {
    size_t endPos = curr_ + k_ - 1;
    return !(endPos + 1 >= rank_->size() or (*rank_)[endPos] == 1);
  }

  inline bool isInvalid() {
    return !isValid();
  }

private:

  void nextValidPosition_() {
    size_t endPos = curr_ + k_ - 1;
    if (endPos < rank_->size()) {
      // See if we cross a rank boundary
      bool crossesBoundary{false};
      bool isNextValid{false};
      size_t boundaryIndex{0};
      for (size_t i = curr_; i < endPos; ++i) {
        if ((*rank_)[i] == 1) {
          crossesBoundary = true;
          boundaryIndex = i;
          isNextValid = (i + k_) < rank_->size();
          break;
        }
      }

      // If so, that's the start of the next valid k-mer
      // if that position is valid
      if (crossesBoundary) {
        if (isNextValid) {
          curr_ = boundaryIndex + 1;
        } else {
          // if that position is invalid, then go to the end.
          goto endPos;
        }
      }
      // At this point, either curr_ points to the next valid
      // start position, or we have skipped over to the endPos label.
      mer_.fromNum(storage_->get_int(2*curr_, 2*k_));
      return;
    }

  endPos:
    // Fallthrough if we couldn't find a valid position.
    mer_.fromNum(storage_->get_int(2*(rank_->size() - k_), 2*k_));
    curr_ = storage_->size() - k_ + 1;
  }

  void advance_() {
    size_t endPos = curr_ + k_ - 1;
    if (endPos + 1 < rank_->size() and (*rank_)[endPos] == 1) {
      curr_ += k_;
      mer_.fromNum(storage_->get_int(2*curr_, 2*k_));
    } else {
      if (curr_ + k_ < rank_->size()) {
        int c = (*storage_)[curr_ + k_];
        mer_.shiftFw(c);
        ++curr_;
      } else {
        mer_.fromNum(storage_->get_int(2*(rank_->size() - k_), 2*k_));
        curr_ = storage_->size() - k_ + 1;
        //curr_ = rank_->size();
      }
    }
  }
  compact::vector<uint64_t, 2>* storage_{nullptr};
  compact::vector<uint64_t, 1>* rank_{nullptr};
  uint8_t k_{0};
  uint64_t curr_{0};
  CanonicalKmer mer_;
  uint64_t word_{0};
};

int pufferfishTest(pufferfish::TestOptions& testOpts) {
  (void)testOpts;
  std::cerr << "this command is not yet implemented\n";
  return 1;
}

void computeSampledPositionsLossy(size_t tlen, uint32_t k, int32_t sampleSize, std::vector<size_t>& sampledInds){
  // Let's start out simple, we always keep the first & last k-mers in a unipath.
  // If the unipath is longer than sampleSize, we keep intermediate samples as well.
  sampledInds.clear();
  auto numOfKmers = tlen - k;
  size_t lastSampled = 0;
  while (lastSampled < numOfKmers) {
    sampledInds.push_back(lastSampled) ;
    auto next_samp = std::min(lastSampled + sampleSize, numOfKmers) ;
    lastSampled = next_samp;
  }
  if (lastSampled ==  numOfKmers) {
    sampledInds.push_back(numOfKmers);
  }
}

void computeSampledPositions(size_t tlen, uint32_t k, int sampleSize, std::vector<size_t>& sampledInds){
  sampledInds.clear() ;
  auto numOfKmers = tlen - k;
  size_t lastCovered = 0 ;
  for(size_t j = 0 ; j <= numOfKmers; j++){
    if(j > lastCovered){
      auto next_samp = std::min(j + sampleSize/2 - 1 ,numOfKmers) ;
      sampledInds.push_back(next_samp) ;
      lastCovered = next_samp + sampleSize/2 + 1;
    }
  }
  if(lastCovered == numOfKmers)
    sampledInds.push_back(lastCovered) ;
}

std::string packedToString(compact::vector<uint64_t, 2>& seqVec, uint64_t offset, uint32_t len) {
  std::stringstream s;
  for (size_t i = offset; i < offset + len; ++i) {
    auto c = seqVec[i];
    s << kmers::charForCode(c);
  }
  auto st = s.str();
  return st;
}

enum class NextSampleDirection : uint8_t { FORWARD = 0, REVERSE=1 };

uint32_t getEncodedExtension(compact::vector<uint64_t, 2>& seqVec, uint64_t firstSampPos, uint64_t distToSamplePos,
                             uint32_t maxExt, NextSampleDirection dir) {
  uint32_t encodedNucs{0};
  uint32_t bitsPerCode{2};
  std::vector<uint32_t> charsToExtend;
  size_t i = 0;
  for (; i < distToSamplePos; ++i) {
    if ( firstSampPos + i >= seqVec.size() ) {
      std::cerr << "seqVec.size() " << seqVec.size() << ", looking for index " << firstSampPos + i << "\n";
      std::cerr << "dist to sample is " << distToSamplePos << ", i = " << i << "\n";
    }
    auto c = seqVec[firstSampPos + i];
    charsToExtend.push_back(c);
  }
  if (dir == NextSampleDirection::REVERSE) {
    std::reverse(charsToExtend.begin(), charsToExtend.end());
  }

  for (size_t j = 0; j < charsToExtend.size(); ++j) {
    auto c = charsToExtend[j];
    encodedNucs |= (c << (bitsPerCode  * (maxExt - j - 1)));
  }
  return encodedNucs;
}

template <typename VecT>
void dumpCompactToFile(VecT& v, std::string fname) {
  std::ofstream bfile(fname, std::ios::binary);
  v.serialize(bfile);
  bfile.close();
}

int fixFastaMain(std::vector<std::string>& args,
                 std::vector<uint32_t>& refIdExtension,
                 std::vector<std::pair<std::string, uint16_t>>& shortRefsNameLen,
                 std::shared_ptr<spdlog::logger> logger, bool hasFeatures);
int buildGraphMain(std::vector<std::string>& args);
int dumpGraphMain(std::vector<std::string>& args);
uint64_t getNumDistinctKmers(unsigned kmlen, const std::string& ifile);

bool copySigArchive(cereal::JSONInputArchive& sigArch, cereal::JSONOutputArchive& indexDesc) {
  std::string seqHash256;
  std::string nameHash256;
  std::string seqHash512;
  std::string nameHash512;
  std::string decoySeqHash256;
  std::string decoyNameHash256;
  uint64_t numberOfDecoys{0};
  uint64_t firstDecoyIndex{std::numeric_limits<uint64_t>::max()};
  bool keep_duplicates{false};

  sigArch( cereal::make_nvp("SeqHash", seqHash256) );
  sigArch( cereal::make_nvp("NameHash", nameHash256) );
  sigArch( cereal::make_nvp("SeqHash512", seqHash512) );
  sigArch( cereal::make_nvp("NameHash512", nameHash512) );
  sigArch( cereal::make_nvp("DecoySeqHash", decoySeqHash256) );
  sigArch( cereal::make_nvp("DecoyNameHash", decoyNameHash256) );
  sigArch( cereal::make_nvp("num_decoys", numberOfDecoys));
  sigArch( cereal::make_nvp("first_decoy_index", firstDecoyIndex));
  sigArch( cereal::make_nvp("keep_duplicates", keep_duplicates));

  indexDesc( cereal::make_nvp("SeqHash", seqHash256) );
  indexDesc( cereal::make_nvp("NameHash", nameHash256) );
  indexDesc( cereal::make_nvp("SeqHash512", seqHash512) );
  indexDesc( cereal::make_nvp("NameHash512", nameHash512) );
  indexDesc( cereal::make_nvp("DecoySeqHash", decoySeqHash256) );
  indexDesc( cereal::make_nvp("DecoyNameHash", decoyNameHash256) );
  indexDesc( cereal::make_nvp("num_decoys", numberOfDecoys));
  indexDesc( cereal::make_nvp("first_decoy_index", firstDecoyIndex));
  indexDesc( cereal::make_nvp("keep_duplicates", keep_duplicates));
  return true;
}

/*
void process_mem_usage(double& vm_usage, double& resident_set)
{
    vm_usage     = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}
*/

int pufferfishIndex(pufferfish::IndexOptions& indexOpts) {
  uint32_t k = indexOpts.k;
  std::vector<std::string> rfiles = indexOpts.rfile;
  std::string rfile;
  std::string outdir = indexOpts.outdir;
  std::string gfa_file = indexOpts.outdir;
  bool buildEdgeVec = indexOpts.buildEdgeVec;
  bool buildEqCls = indexOpts.buildEqCls;
  bool keepFixedFasta = indexOpts.keep_fixed_fasta;
  std::vector<uint32_t> refIdExtensions;
  std::vector<std::pair<std::string, uint16_t>> shortRefsNameLen;

  // If the user included the '/' in the output directory path, remove
  // it here
  if (outdir.back() == '/') {
    outdir.pop_back();
  }

  size_t tlen{0};
  size_t numKmers{0};
  size_t nread{0};
  CanonicalKmer::k(k);

  if (ghc::filesystem::exists(outdir.c_str())) {
      if (!ghc::filesystem::is_directory(outdir.c_str())) {
          auto console = spdlog::stderr_color_mt("console");
          console->error("{} exists as a file. Cannot create a directory of the same name.", outdir.c_str());
          std::exit(1);
      }
  } else {
      ghc::filesystem::create_directories(outdir.c_str());
  }

  std::string logPath = outdir + "/ref_indexing.log";
  auto fileSink = std::make_shared<spdlog::sinks::simple_file_sink_st>(logPath);
  auto consoleSink = std::make_shared<spdlog::sinks::ansicolor_stderr_sink_mt>();
  consoleSink->set_color(spdlog::level::warn, consoleSink->magenta);
  auto consoleLog = spdlog::create("puff::index::stderrLog", {consoleSink});
  auto fileLog = spdlog::create("puff::index::fileLog", {fileSink});
  std::vector<spdlog::sink_ptr> sinks{consoleSink, fileSink};
  auto jointLog = spdlog::create("puff::index::jointLog", std::begin(sinks), std::end(sinks));

  /*if (puffer::fs::MakePath(outdir.c_str()) != 0) {
      std::cerr << "\nyup that's it\n";
    jointLog->error(std::strerror(errno));
    std::exit(1);
  }*/

  // running fixFasta
  {
    jointLog->info("Running fixFasta");
    std::vector<std::string> args;
//    args.push_back("fixFasta");
    if (!indexOpts.decoy_file.empty()) {
      args.push_back("--decoys");
      args.push_back(indexOpts.decoy_file);
    }
    if (!indexOpts.header_sep.empty()) {
      args.push_back("--headerSep");
      args.push_back(indexOpts.header_sep);
    }
    if (indexOpts.keep_duplicates) {
      args.push_back("--keepDuplicates");
    }
    if (indexOpts.expect_transcriptome) {
      args.push_back("--expectTranscriptome");
    }
    if (indexOpts.noclip_polya) {
      args.push_back("--noClip");
    }
    args.push_back("--klen");
    args.push_back(std::to_string(k));
    args.push_back("--input");
    args.insert(args.end(), rfiles.begin(), rfiles.end());
    args.push_back("--output");
    args.push_back(outdir+"/ref_k"+std::to_string(k)+"_fixed.fa");

    int ffres = fixFastaMain(args, refIdExtensions, shortRefsNameLen, jointLog, indexOpts.featuresRef);

    if (ffres != 0) {
        jointLog->error("The fixFasta phase failed with exit code {}", ffres);
        std::exit(ffres);
    }
    // replacing rfile with the new fixed fasta file
    rfile = outdir+"/ref_k"+std::to_string(k)+"_fixed.fa";
  }

  //std::this_thread::sleep_for (std::chrono::seconds(10));
  //double vm, rss;
  //process_mem_usage(vm, rss);
  //std::cerr << "\n\n after fix fasta \n";
  //std::cerr << "VM: " << vm << "; RSS: " << rss << "\n\n";

  // If the filter size isn't set by the user, estimate it with ntCard
  if (indexOpts.filt_size == -1){
    jointLog->info("Filter size not provided; estimating from number of distinct k-mers");
    auto nk = getNumDistinctKmers(k, rfile);
    double p = 0.001;
    double k = 5.0;
    double logp_k = std::log(p) / k;
    double r = (-k) / std::log(1.0 - std::exp(logp_k));
    indexOpts.filt_size = static_cast<int32_t>(std::ceil(std::log2(std::ceil(nk * r))));
    jointLog->info("ntHll estimated {} distinct k-mers, setting filter size to 2^{}", nk, indexOpts.filt_size);
    /*
    lgp_k = l($p) / $k
    r = (-$k) / l(1 - e(lgp_k))
    ceil(l(ceil($n * r))/l(2))
    */
  }

  {
    std::vector<std::string> args;
    args.push_back("twopaco");
    args.push_back("-k");
    args.push_back(std::to_string(k));
    args.push_back("-t");
    args.push_back(std::to_string(indexOpts.p));
    args.push_back("-f");
    args.push_back(std::to_string(indexOpts.filt_size));
    args.push_back("--outfile");
    args.push_back(outdir+"/tmp_dbg.bin");
    args.push_back("--tmpdir");

    std::string twopaco_tmp_path = indexOpts.twopaco_tmp_dir;
    // if the tmp path wasn't set, then use a subdirectory 
    // of the index directory (that we will later remove).
    if (twopaco_tmp_path.empty()) {
      twopaco_tmp_path = outdir + "/twopaco_tmp";
    } 

    // create the tmp directory if we need to (and can). Complain and exit
    // if the user passed an existing file as the target path. 
    if (ghc::filesystem::exists(twopaco_tmp_path.c_str())) {
        if (!ghc::filesystem::is_directory(twopaco_tmp_path.c_str())) {
            jointLog->error("{} exists as a file. Cannot create a directory of the same name.", twopaco_tmp_path.c_str());
            jointLog->flush();
            std::exit(1);
        }
    } else {
        ghc::filesystem::create_directories(twopaco_tmp_path.c_str());
    }
    args.push_back(twopaco_tmp_path);

    args.push_back(rfile);
    buildGraphMain(args);

    // cleanup tmp
    ghc::filesystem::remove_all(twopaco_tmp_path);
  }

  //std::this_thread::sleep_for (std::chrono::seconds(10));
  //double vm, rss;
  //process_mem_usage(vm, rss);
  //std::cerr << "\n\n after build graph \n";
  //std::cerr << "VM: " << vm << "; RSS: " << rss << "\n\n";

  {
    std::vector<std::string> args;
    args.push_back("graphdump");
    args.push_back("-k");
    args.push_back(std::to_string(k));
    args.push_back("-s");
    args.push_back(rfile);
    args.push_back("-f");
    args.push_back("binPufferized");
    args.push_back(outdir+"/tmp_dbg.bin");
    args.push_back("-p");
    args.push_back(outdir);
    dumpGraphMain(args);

    // cleanup what we no longer need
    ghc::filesystem::path outpath{outdir};
    ghc::filesystem::path tmpDBG = outdir / ghc::filesystem::path{"tmp_dbg.bin"};
    if (ghc::filesystem::exists(tmpDBG)) {
      ghc::filesystem::remove(tmpDBG);
    }
  }

  //std::this_thread::sleep_for (std::chrono::seconds(10));
  //process_mem_usage(vm, rss);
  //std::cerr << "\n\n after graph dump \n";
  //std::cerr << "VM: " << vm << "; RSS: " << rss << "\n\n";


  jointLog->info("Starting the Pufferfish indexing by reading the GFA binary file.");
  pufferfish::BinaryGFAReader pf(outdir.c_str(), k - 1, buildEqCls, buildEdgeVec, jointLog);
  pf.parseFile();
  pf.mapContig2Pos();
  pf.serializeContigTable(outdir, shortRefsNameLen, refIdExtensions);
  {
    auto& cnmap = pf.getContigNameMap();
    for (auto& kv : cnmap) {
      const auto& r1 = kv.second;
      tlen += r1.length;
      numKmers += r1.length - k + 1;
      ++nread;
    }
    jointLog->info("# segments = {:n}", nread);
    jointLog->info("total length = {:n}", tlen);
  }

  // parse the reference list and store the strings in a 2bit-encoded vector
  bool keepRef = (rfile.size() > 0);
  if (keepRef) {
    auto &refIds = pf.getRefIDs();
    auto &refLengths = pf.getRefLengths();
    std::vector<uint64_t> refAccumLengths(refLengths.size());
    std::unordered_map<std::string, uint64_t> refIdMap;
    uint64_t prev = 0;
    for (uint64_t i = 0; i < refIds.size(); i++) {
      if (refIdMap.find(refIds[i]) != refIdMap.end()) {
        jointLog->error("Two references with the same name but different sequences: {}. "
                       "We require that all input records have a unique name "
                       "up to the first whitespace character.", refIds[i]);
        std::exit(1);
      }
      refIdMap[refIds[i]] = i;
      refAccumLengths[i] = prev + refLengths[i];
      prev = refAccumLengths[i];
//      std::cerr << i << ":" << refLengths[i] << "\n";
    }
    //compact 2bit vector
//    std::cerr << "\nrefAccumLengths.size():" << refAccumLengths.size() << "\n";
//    std::cerr << refAccumLengths.back() << "\n";
    compact::vector<uint64_t, 2> refseq(refAccumLengths.back());
    refseq.clear_mem();

    // go over all the reference files
    jointLog->info("Reading the reference files ...");
    std::vector<std::string> ref_files = {rfile};
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(ref_files, 1, 1);
    parser.start();
    auto rg = parser.getReadGroup();
    // read the reference sequences and encode them into the refseq int_vector

    while (parser.refill(rg)) {
      for (auto &rp : rg) {
        stx::string_view seqv(rp.seq);
        auto &refIdx = refIdMap[rp.name];
        auto offset = refIdx == 0 ? 0 : refAccumLengths[refIdx - 1];
        //std::cerr << "ref from [" << offset << ", " << offset + seqv.length() << "]\n";
        pf.encodeSeq(refseq, offset, seqv);
      }
    }

    parser.stop();
    // store the 2bit-encoded references
    // store reference accumulative lengths
    std::string accumLengthsFilename = outdir + "/refAccumLengths.bin";
    std::ofstream ral(accumLengthsFilename);
    cereal::BinaryOutputArchive ralAr(ral);
    ralAr(refAccumLengths);

    // store reference sequences
    std::ofstream seqFile(outdir + "/refseq.bin", std::ios::binary);
    refseq.serialize(seqFile);
    seqFile.close();
  }
  pf.clearContigTable();
  // At this point we should definitely not need path.bin anymore, so get rid of it
  ghc::filesystem::path tmpPath = ghc::filesystem::path{outdir} / ghc::filesystem::path{"path.bin"};
  if (!ghc::filesystem::remove(tmpPath)) {
    jointLog->warn("Could not seem to remove temporary file {}.", tmpPath.string());
  }

  // now we know the size we need --- create our bitvectors and pack!
  size_t w = std::log2(tlen) + 1;
  jointLog->info("positional integer width = {:n}", w);

  auto& seqVec = pf.getContigSeqVec();
  auto& rankVec = pf.getRankVec();
  auto& edgeVec = pf.getEdgeVec() ;

 /* compact::vector<uint64_t, 1> rankVec(tlen);
  for(uint64_t i=0; i<tlen; i++) rankVec[i]=0;//if(rankVec[i]!=0) {std::cerr<<"Not zero\n"; break;}
  auto& cnmap = pf.getContigNameMap();
  for (auto& kv : cnmap) {
    rankVec[kv.second.offset + kv.second.length - 1] = 1;
  }*/
  size_t nkeys{numKmers};
  size_t numContigs{pf.getContigNameMap().size()};

  jointLog->info("seqSize = {:n}", seqVec.size());
  jointLog->info("rankSize = {:n}", rankVec.size());

  jointLog->info("edgeVecSize = {:n}", edgeVec.size());

  jointLog->info("num keys = {:n}", nkeys);
  ContigKmerIterator kb(&seqVec, &rankVec, k, 0);
  ContigKmerIterator ke(&seqVec, &rankVec, k, seqVec.size() - k + 1);

#ifdef PUFFER_DEBUG
  auto ks = kb;
  size_t nkeyIt{0};
  for (; ks < ke; ++ks) {
    nkeyIt++;
  }
  jointLog->info("num keys (iterator)= {:n}", nkeyIt);
#endif // PUFFER_DEBUG
 
  typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
  typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;

  auto keyIt = boomphf::range(kb, ke);
  std::unique_ptr<boophf_t> bphf = std::make_unique<boophf_t>(outdir, nkeys, keyIt, indexOpts.p, 3.5); // keys.size(), keys, 16);
  jointLog->info("mphf size = {} MB", (bphf->totalBitSize() / 8) / std::pow(2, 20));

/*  std::ofstream seqFile(outdir + "/seq.bin", std::ios::binary);
  seqVec.serialize(seqFile);
  seqFile.close();

  std::ofstream rankFile(outdir + "/rank.bin", std::ios::binary);
  rankVec.serialize(rankFile);
  rankFile.close();*/

  bool haveEdgeVec = edgeVec.size() > 0;
  if (haveEdgeVec) {
    std::ofstream edgeFile(outdir + "/edge.bin", std::ios::binary);
    edgeVec.serialize(edgeFile);
    edgeFile.close();
  }

  // if using quasi-dictionary idea (https://arxiv.org/pdf/1703.00667.pdf)
  //uint32_t hashBits = 4;
  if (!indexOpts.isSparse and !indexOpts.lossySampling) {  
    // if using quasi-dictionary idea (https://arxiv.org/pdf/1703.00667.pdf)
    compact::ts_vector<uint64_t> posVec(w, nkeys);
    posVec.clear_mem();
    {

      struct ContigVecChunk {
        uint64_t s;
        uint64_t e;
      };

      // Build position table in parallel. We have up to indexOpts.p threads
      // so divide the contig array into even chunks.
      auto nthread = indexOpts.p;
      double chunkSizeFrac = seqVec.size() / static_cast<double>(nthread);
      // reduce the number of threads until chunks are big enough;
      while (chunkSizeFrac < 8192 and nthread > 1) {
        nthread /= 2;
        chunkSizeFrac = seqVec.size() / static_cast<double>(nthread);
      }

      auto chunkSize = static_cast<uint64_t>(std::ceil(chunkSizeFrac));
      jointLog->info("chunk size = {:n}", chunkSize);

      std::vector<ContigVecChunk> chunks;
      chunks.reserve(nthread);
      size_t itOffset{0};
      for (size_t i = 0; i < nthread; ++i) {
        ContigKmerIterator startIt(&seqVec, &rankVec, k, itOffset);
        startIt.advanceToValid();
        uint64_t s = startIt.pos();
        if (i == nthread - 1) {
          itOffset = seqVec.size() - k + 1;
        } else {
          itOffset = s + chunkSize;
        }
        ContigKmerIterator endIt(&seqVec, &rankVec, k, itOffset);
        endIt.advanceToValid();
        uint64_t e = endIt.pos();
        chunks.push_back({s,e});
        jointLog->info("chunk {} = [{:n}, {:n})", i, s, e);
      }

      auto fillPos = [&seqVec, &rankVec, k, &bphf, &jointLog, &posVec](ContigVecChunk chunk) -> void {

        ContigKmerIterator kb1(&seqVec, &rankVec, k, chunk.s);
        ContigKmerIterator ke1(&seqVec, &rankVec, k, chunk.e);
        for (; kb1 < ke1; ++kb1) {
          auto idx = bphf->lookup(*kb1); // fkm.word(0));
          if (idx >= posVec.size()) {
            std::cerr<<*kb1<<"\n";
            jointLog->info("seq size = {:n}, idx = {:n}, pos size = {:n}",
                          seqVec.size(), idx, posVec.size());
            std::cerr<<*kb1<<"\n";
          }
          // ContigKmerIterator::value_type mer = *kb1;
          // if using quasi-dictionary idea (https://arxiv.org/pdf/1703.00667.pdf)
          //posVec[idx] = (kb1.pos() << hashBits) | (mer & 0xF);
          posVec[idx] = kb1.pos();
          // validate
#ifdef PUFFER_DEBUG
          uint64_t kn = seqVec.get_int(2*kb1.pos(), 2*k);
          CanonicalKmer sk;
          sk.fromNum(kn);
          if (sk.isEquivalent(*kb1) == KmerMatchType::NO_MATCH) {
            my_mer r;
            r.word__(0) = *kb1;
            jointLog->error("I thought I saw {}, but I saw {} --- pos {}", sk.to_str(), r.toStr(), kb1.pos());
          }
#endif
        }
      };

      std::vector<std::thread> workers;
      workers.reserve(chunks.size());
      for (auto chunk : chunks) {
        workers.push_back(std::thread(fillPos, chunk));
      }
      for (auto& w : workers) {
        w.join();
      }
      jointLog->info("finished populating pos vector");

    }

    jointLog->info("writing index components");
    /** Write the index **/



    std::ofstream descStream(outdir + "/info.json");
    {
      cereal::JSONOutputArchive indexDesc(descStream);
      std::string sampStr = "dense";
      std::vector<std::string> refGFA{outdir};
      indexDesc(cereal::make_nvp("index_version", pufferfish::indexVersion));
      indexDesc(cereal::make_nvp("reference_gfa", refGFA));
      indexDesc(cereal::make_nvp("sampling_type", sampStr));
      indexDesc(cereal::make_nvp("k", k));
      indexDesc(cereal::make_nvp("num_kmers", nkeys));
      indexDesc(cereal::make_nvp("num_contigs", numContigs));
      indexDesc(cereal::make_nvp("seq_length", tlen));
      indexDesc(cereal::make_nvp("have_ref_seq", keepRef));
      indexDesc(cereal::make_nvp("have_edge_vec", haveEdgeVec));

      std::ifstream sigStream(outdir + "/ref_sigs.json");
      cereal::JSONInputArchive sigArch(sigStream);
      copySigArchive(sigArch, indexDesc);
      sigStream.close();
    }
    descStream.close();

    std::ofstream posFile(outdir + "/pos.bin", std::ios::binary);
    posVec.serialize(posFile);
    posFile.close();
    std::ofstream hstream(outdir + "/mphf.bin");
    bphf->save(hstream);
    hstream.close();
    jointLog->info("finished writing dense pufferfish index");

  } else if (indexOpts.isSparse) { // sparse index; it's GO time!
    int extensionSize = indexOpts.extensionSize;
    int sampleSize = 2 * extensionSize + 1;


    // Note: the compact_vector constructor does not
    // init mem to 0, so we do that with the clear_mem() function.
    compact::vector<uint64_t, 1> presenceVec(nkeys);
    presenceVec.clear_mem();

    size_t sampledKmers{0};
    std::vector<size_t> sampledInds;
    std::vector<size_t> contigLengths;
    //fill up optimal positions
    {
      auto& cnmap = pf.getContigNameMap() ;
      size_t ncontig = cnmap.size();
      std::vector<size_t> sampledInds ;
      for(size_t i = 0; i < ncontig; ++i) {//}auto& kv : cnmap){
        const auto& r1 = cnmap[i];
        sampledInds.clear();
        computeSampledPositions(r1.length, k, sampleSize, sampledInds) ;
        sampledKmers += sampledInds.size() ;
        contigLengths.push_back(r1.length) ;
      }
      jointLog->info("# sampled kmers = {:n}", sampledKmers) ;
      jointLog->info("# skipped kmers = {:n}", numKmers - sampledKmers) ;
    }

    //fill up the vectors
    uint32_t extSymbolWidth = 2;
    uint32_t extWidth = std::log2(extensionSize);
    jointLog->info("extWidth = {}", extWidth);

    compact::vector<uint64_t> auxInfo(extSymbolWidth*extensionSize, (numKmers-sampledKmers));
    auxInfo.clear_mem();

    compact::vector<uint64_t> extSize(extWidth, (numKmers-sampledKmers));
    extSize.clear_mem();

    compact::vector<uint64_t, 1> direction(numKmers - sampledKmers) ;
    direction.clear_mem();

    compact::vector<uint64_t, 1> canonicalNess(numKmers - sampledKmers);
    canonicalNess.clear_mem();

    compact::vector<uint64_t> samplePosVec(w, sampledKmers);
    samplePosVec.clear_mem();

  // new presence Vec
    size_t i = 0 ;
    std::unordered_set<uint64_t> indices;
  {
    jointLog->info("\nFilling presence Vector");

    ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
    ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);
    size_t contigId{0};

    //debug flags
    int loopCounter = 0;

    // walk over the entire contig array:
    // compute the sampled positions for each contig
    // fill in the corresponding values in presenceVec
    while(kb1 != ke1){
        sampledInds.clear();
        auto clen = contigLengths[contigId];
        computeSampledPositions(clen, k, sampleSize, sampledInds) ;
        contigId++;
        loopCounter++ ;

        my_mer r;
        auto zeroPos = kb1.pos();
        auto skipLen = kb1.pos() - zeroPos;
        auto nextSampIter = sampledInds.begin();
        //auto prevSamp = *nextSampIter;
        //bool didSample = false;
        bool done = false;

        for (size_t j = 0; j < clen - k + 1; ++kb1, ++j) {
          skipLen = kb1.pos() - zeroPos;
          if (!done and skipLen == static_cast<decltype(skipLen)>(*nextSampIter)) {
            auto idx = bphf->lookup(*kb1);
            presenceVec[idx] = 1 ;
            indices.insert(idx);
            i++ ;
            //didSample = true;
            //prevSamp = *nextSampIter;
            ++nextSampIter;
            if (nextSampIter == sampledInds.end()) {
              done = true;
            }
          }
          //didSample = false;
        }
        if (nextSampIter != sampledInds.end()) {
          jointLog->info("I didn't sample {}, samples for contig {}", std::distance(nextSampIter, sampledInds.end()), contigId - 1);
          jointLog->info("last sample is {}" , sampledInds.back());
          jointLog->info("contig length is {}" , contigLengths[contigId-1]);
        }
    }

    jointLog->info("i = {:n}, sampled kmers = {:n}, loops = {:n}, contig array = {:n}",
                  i, sampledKmers, loopCounter, contigLengths.size());
  }

  rank9b realPresenceRank(presenceVec.get(), presenceVec.size());
  jointLog->info("num ones in presenceVec = {:n}, i = {:n}, indices.size() = {:n}", realPresenceRank.rank(presenceVec.size()-1), i, indices.size());

  //bidirectional sampling
  {

    ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
    ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);

    size_t contigId{0} ;
    //size_t coveredKeys{0} ;
    size_t totalKmersIshouldSee{0} ;

    // For every valid k-mer (i.e. every contig)
    while(kb1 != ke1){
      sampledInds.clear();
      auto clen = contigLengths[contigId];
      auto thisContigLength = clen;
      computeSampledPositions(clen, k, sampleSize, sampledInds) ;
      totalKmersIshouldSee += (thisContigLength - k + 1);

      contigId++ ;

      //size_t skip = 0 ;

      my_mer r;

      auto zeroPos = kb1.pos();
      auto nextSampIter = sampledInds.begin();
      auto prevSampIter = sampledInds.end();
      auto skipLen = kb1.pos() - zeroPos;
      NextSampleDirection sampDir = NextSampleDirection::FORWARD;
      bool done = false;
      for (size_t j = 0; j < clen - k + 1; ++kb1, ++j) {
          int64_t nextSampPos = (nextSampIter != sampledInds.end()) ? *nextSampIter : -1;
          int64_t prevSampPos = (prevSampIter != sampledInds.end()) ? *prevSampIter : -1;
          uint64_t distToNext = (nextSampPos >= 0) ? nextSampPos - j : std::numeric_limits<uint64_t>::max();
          uint64_t distToPrev = (prevSampPos >= 0) ? j - prevSampPos : std::numeric_limits<uint64_t>::max();

          if (distToNext == std::numeric_limits<uint64_t>::max() and
              distToPrev == std::numeric_limits<uint64_t>::max()) {
            jointLog->error("Could not find valid sample position, should not happen!");
            std::exit(1);
          }

          sampDir = (distToNext < distToPrev) ? NextSampleDirection::FORWARD : NextSampleDirection::REVERSE;
          skipLen = kb1.pos() - zeroPos;
          // If this is a sampled position
          if (!done and skipLen == static_cast<decltype(skipLen)>(*nextSampIter)) {
            prevSampIter = nextSampIter;
            ++nextSampIter;
            if (nextSampIter == sampledInds.end()) {
              done = true;
            }
            auto idx = bphf->lookup(*kb1);
            auto rank = (idx == 0) ? 0 : realPresenceRank.rank(idx);
            samplePosVec[rank] = kb1.pos();
          } else { // not a sampled position
            uint32_t ext = 0;
            size_t firstSampPos = 0;
            uint32_t extensionDist = 0;
            if (sampDir == NextSampleDirection::FORWARD) {
              firstSampPos = zeroPos + j + k;
              extensionDist = distToNext - 1;
              ext = getEncodedExtension(seqVec, firstSampPos, distToNext, extensionSize, sampDir);
            } else if (sampDir == NextSampleDirection::REVERSE) {
              firstSampPos = zeroPos + prevSampPos;
              extensionDist = distToPrev - 1;
              ext = getEncodedExtension(seqVec, firstSampPos, distToPrev, extensionSize, sampDir);
            } else {
              std::cerr << "Error during extension encoding, should not happen!\n";
              std::exit(1);
            }
            auto idx = bphf->lookup(*kb1);
            auto rank = (idx == 0) ? 0 : realPresenceRank.rank(idx);

            int64_t target_idx = (idx - rank);
            if ( target_idx > static_cast<int64_t>(canonicalNess.size())) { jointLog->warn("target_idx = {}, but canonicalNess.size = {}", target_idx, canonicalNess.size()); }
            canonicalNess[idx - rank] = kb1.isCanonical();

            if ( target_idx > static_cast<int64_t>(extSize.size())) { jointLog->warn("target_idx = {}, but extSize.size = {}", target_idx, extSize.size()); }
            extSize[idx - rank] = extensionDist;

            if ( target_idx > static_cast<int64_t>(auxInfo.size())) { jointLog->warn("target_idx = {}, but auxInfo.size = {}", target_idx, auxInfo.size()); }
            auxInfo[idx - rank] = ext;

            if ( target_idx > static_cast<int64_t>(direction.size())) { jointLog->warn("target_idx = {}, but direction.size = {}", target_idx, direction.size()); }
            direction[idx - rank] = (sampDir == NextSampleDirection::FORWARD) ? 1 : 0;
          }
        }
    }


  }


  /** Write the index **/
  std::ofstream descStream(outdir + "/info.json");
  {
    cereal::JSONOutputArchive indexDesc(descStream);
    std::string sampStr = "sparse";
    std::vector<std::string> refGFA{outdir};
    indexDesc(cereal::make_nvp("index_version", pufferfish::indexVersion));
    indexDesc(cereal::make_nvp("reference_gfa", refGFA));
    indexDesc(cereal::make_nvp("sampling_type", sampStr));
    indexDesc(cereal::make_nvp("sample_size", sampleSize));
    indexDesc(cereal::make_nvp("extension_size", extensionSize));
    indexDesc(cereal::make_nvp("k", k));
    indexDesc(cereal::make_nvp("num_kmers", nkeys));
    indexDesc(cereal::make_nvp("num_sampled_kmers",sampledKmers));
    indexDesc(cereal::make_nvp("num_contigs", numContigs));
    indexDesc(cereal::make_nvp("seq_length", tlen));
    indexDesc(cereal::make_nvp("have_ref_seq", keepRef));
    indexDesc(cereal::make_nvp("have_edge_vec", haveEdgeVec));

    std::ifstream sigStream(outdir + "/ref_sigs.json");
    cereal::JSONInputArchive sigArch(sigStream);
    copySigArchive(sigArch, indexDesc);
    sigStream.close();
  }
  descStream.close();

  std::ofstream hstream(outdir + "/mphf.bin");
  dumpCompactToFile(presenceVec, outdir+"/presence.bin");
  dumpCompactToFile(samplePosVec, outdir + "/sample_pos.bin");
  dumpCompactToFile(auxInfo, outdir + "/extension.bin");
  dumpCompactToFile(extSize, outdir + "/extensionSize.bin");
  dumpCompactToFile(canonicalNess, outdir + "/canonical.bin");
  dumpCompactToFile(direction, outdir + "/direction.bin");
  bphf->save(hstream);
  hstream.close();

  } else { // lossy sampling index
    int32_t sampleSize = static_cast<int32_t>(indexOpts.lossy_rate);
    compact::vector<uint64_t, 1> presenceVec(nkeys);
    presenceVec.clear_mem();

    size_t sampledKmers{0};
    std::vector<size_t> sampledInds;
    std::vector<size_t> contigLengths;
    //fill up optimal positions
    {
      auto& cnmap = pf.getContigNameMap() ;
      std::vector<size_t> sampledInds ;
      for(auto& kv : cnmap){
        auto& r1 = kv.second ;
        sampledInds.clear();
        computeSampledPositionsLossy(r1.length, k, sampleSize, sampledInds) ;
        sampledKmers += sampledInds.size() ;
        contigLengths.push_back(r1.length) ;
      }
      jointLog->info("# sampled kmers = {:n}", sampledKmers) ;
      jointLog->info("# skipped kmers = {:n}", numKmers - sampledKmers) ;
    }

    compact::vector<uint64_t> samplePosVec(w, sampledKmers);
    samplePosVec.clear_mem();



    // new presence Vec
    {
      {
      jointLog->info("\nFilling presence vector");
      size_t i = 0 ;
      ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
      ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);
      size_t contigId{0};
      int loopCounter = 0;
      while(kb1 != ke1){
        sampledInds.clear();
        auto clen = contigLengths[contigId];
        computeSampledPositionsLossy(clen, k, sampleSize, sampledInds) ;
        contigId++;
        loopCounter++ ;

        my_mer r;
        auto zeroPos = kb1.pos();
        auto skipLen = kb1.pos() - zeroPos;
        auto nextSampIter = sampledInds.begin();
        bool done = false;
        
        for (size_t j = 0; j < clen - k + 1; ++kb1, ++j) {
          skipLen = kb1.pos() - zeroPos;
          if (!done and skipLen == static_cast<decltype(skipLen)>(*nextSampIter)) {
            auto idx = bphf->lookup(*kb1);
            presenceVec[idx] = 1 ;
            samplePosVec[i] = kb1.pos();
            i++;
            ++nextSampIter;
            if (nextSampIter == sampledInds.end()) {
              done = true;
            }
          }
        }
        if (nextSampIter != sampledInds.end()) {
          jointLog->info("I didn't sample {:n}, samples for contig {:n}", std::distance(nextSampIter, sampledInds.end()), contigId - 1);
          jointLog->info("last sample is {:n}", sampledInds.back());
          jointLog->info("contig length is {:n}", contigLengths[contigId-1]);
        }
      }
      }

      {
        jointLog->info("\nFilling sampled position vector");
      size_t i = 0 ;
      ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
      ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);
      size_t contigId{0};
      int loopCounter = 0;
      rank9b realPresenceRank(presenceVec.get(), presenceVec.size());
      while(kb1 != ke1){
        sampledInds.clear();
        auto clen = contigLengths[contigId];
        computeSampledPositionsLossy(clen, k, sampleSize, sampledInds) ;
        contigId++;
        loopCounter++ ;

        my_mer r;
        auto zeroPos = kb1.pos();
        auto skipLen = kb1.pos() - zeroPos;
        auto nextSampIter = sampledInds.begin();
        bool done = false;
        
        for (size_t j = 0; j < clen - k + 1; ++kb1, ++j) {
          skipLen = kb1.pos() - zeroPos;
          if (!done and skipLen == static_cast<decltype(skipLen)>(*nextSampIter)) {
            auto idx = bphf->lookup(*kb1);
            auto rank = (idx == 0) ? 0 : realPresenceRank.rank(idx);
            samplePosVec[rank] = kb1.pos();
            ++i;
            ++nextSampIter;
            if (nextSampIter == sampledInds.end()) {
              done = true;
            }
          }
        }
        if (nextSampIter != sampledInds.end()) {
          jointLog->info("I didn't sample {:n}, samples for contig {}", std::distance(nextSampIter, sampledInds.end()), contigId - 1);
          jointLog->info("last sample is {:n}", sampledInds.back());
          jointLog->info("contig length is {:n}", contigLengths[contigId-1]);
        }
      }
      jointLog->info("i = {:n}, sampled kmers = {:n}, loops = {:n}, contig array = {:n}",
                    i, sampledKmers, loopCounter, contigLengths.size());

    }
    }

    /** Write the index **/
    std::ofstream descStream(outdir + "/info.json");
    {
      cereal::JSONOutputArchive indexDesc(descStream);
      std::string sampStr = "lossy";
      std::vector<std::string> refGFA{outdir};
      indexDesc(cereal::make_nvp("index_version", pufferfish::indexVersion));
      indexDesc(cereal::make_nvp("reference_gfa", refGFA));
      indexDesc(cereal::make_nvp("sampling_type", sampStr));
      indexDesc(cereal::make_nvp("sample_size", sampleSize));
      indexDesc(cereal::make_nvp("k", k));
      indexDesc(cereal::make_nvp("num_kmers", nkeys));
      indexDesc(cereal::make_nvp("num_sampled_kmers",sampledKmers));
      indexDesc(cereal::make_nvp("num_contigs", numContigs));
      indexDesc(cereal::make_nvp("seq_length", tlen));
      indexDesc(cereal::make_nvp("have_ref_seq", keepRef));
      indexDesc(cereal::make_nvp("have_edge_vec", haveEdgeVec));

      std::ifstream sigStream(outdir + "/ref_sigs.json");
      cereal::JSONInputArchive sigArch(sigStream);
      copySigArchive(sigArch, indexDesc);
      sigStream.close();
    }
    descStream.close();

    std::ofstream hstream(outdir + "/" + pufferfish::util::MPH);
    dumpCompactToFile(presenceVec, outdir + "/presence.bin");
    dumpCompactToFile(samplePosVec, outdir + "/sample_pos.bin");
    bphf->save(hstream);
    hstream.close();
  }

  // cleanup the fixed.fa file
  if (!keepFixedFasta) { ghc::filesystem::remove(rfile); }
  ghc::filesystem::remove(outdir + "/ref_sigs.json");
  return 0;
}
