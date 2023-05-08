#ifndef __PUFFERFISH_PROG_OPTS_HPP__
#define __PUFFERFISH_PROG_OPTS_HPP__
#include <cstdint>
#include <string>
#include <vector>

namespace pufferfish {
class IndexOptions {
public:
  uint32_t k{31};
  uint32_t p{16};
  std::string gfa_file;
  std::string cfile;
  std::vector<std::string> rfile;
  std::string outdir;
  std::string decoy_file{""};
  std::string header_sep{""};
  bool noclip_polya{false};
  bool isSparse{false};
  uint32_t extensionSize{4};
  uint32_t sampleSize{9};
  bool lossySampling{false};
  bool keep_fixed_fasta{false};
  bool keep_duplicates{false};
  uint32_t lossy_rate{5};
  int32_t filt_size{-1};
  bool buildEdgeVec{false};
  bool buildEqCls{false};
  bool featuresRef{false};
  bool expect_transcriptome{false};
  std::string twopaco_tmp_dir{""};
};

class ExamineOptions {
public:
  std::string index_dir{""};
  std::string fasta_out{""};
  std::string kmer_freq_out{""};
};

class TestOptions {
public:
};

enum StatType {
    ctab, motif
};
class StatsOptions {
public:
    StatType statType{ctab};
    std::string indexDir;
};

class KmerQueryOptions {
  public:
  std::string indexDir;
  std::vector<std::string> queryFiles;
  uint32_t num_threads{16};
};

class ValidateOptions {
public:
  std::string indexDir;
  std::string refFile;
  std::string gfaFileName ;
};

class AlignmentOpts{
public:

  std::string indexDir;
  std::string read1;
  std::string read2;
  std::string unmatedReads;
  bool singleEnd{false};
  uint32_t numThreads{1};
  uint32_t maxNumHits{200};
  uint32_t maxSpliceGap{100};
  uint32_t maxFragmentLength{1000};
  double scoreRatio{0.6};
  double consensusFraction{0.65};
  uint32_t altSkip{5};
  std::string outname;
  double quasiCov{0.0};
  bool pairedEnd{false};
  bool noOutput{false};
  bool sensitive{false};
  bool strictCheck{false};
  bool fuzzy{false};
  bool consistentHits{false};
  bool quiet{false};
  bool writeOrphans{false} ;
  bool justMap{false};
  bool krakOut{false};
  bool salmonOut{false};
  bool radOut{false};
  bool noDiscordant{false};
  bool noOrphan{false};
  bool noDovetail{false};
  bool compressedOutput{false};
  bool verbose{false};
  bool validateMappings{true};
  bool bestStrata{false};
  bool writeQualities{false};
  int32_t gapOpenPenalty{5};
  int32_t gapExtendPenalty{3};
  int32_t matchScore{2};
  int32_t mismatchScore{-4};
  uint32_t refExtendLength{20};
  double minScoreFraction{0.65};
  bool fullAlignment{false};
  bool heuristicChaining{true};
  bool genomicReads{false};
  std::string genesNamesFile{""};
  std::string rrnaFile{""};
  bool filterGenomics{false};
  bool filterMicrobiomBestScore{false};
  bool filterMicrobiom{false};
  bool filterRrna{false};
  bool primaryAlignment{false};
  bool listOfReads{false};
  uint32_t maxAllowedRefsPerHit{1000};
  bool recoverOrphans{false};
  bool mimicBt2Default{false};
  bool mimicBt2Strict{false};
  bool allowOverhangSoftclip{false};
  bool allowSoftclip{false};
  bool useAlignmentCache{true};
  uint32_t alignmentStreamLimit{10000};
  double preMergeChainSubThresh{0.9};
  double postMergeChainSubThresh{0.9};
  double orphanChainSubThresh{1};
};
}

#endif
