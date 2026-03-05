#include <cstdio>
#include <iostream>
#include <random>
#include <unistd.h>
#include <unordered_map>

#include "FastxParser.hpp"
#include "jellyfish/mer_dna.hpp"

#include "FASTAParser.hpp"
#include "SalmonOpts.hpp"
#include "SalmonStringUtils.hpp"
#include "Transcript.hpp"

FASTAParser::FASTAParser(const std::string& fname) : fname_(fname) {}

void FASTAParser::populateTargets(std::vector<Transcript>& refs,
                                  SalmonOpts& sopt) {
  using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

  using std::string;
  using std::unordered_map;

  unordered_map<string, size_t> nameToID;
  for (auto& ref : refs) {
    nameToID[ref.RefName] = ref.id;
  }

  // Separators for the header (default ' ' and '\t')
  // If we have the gencode flag, then add '|'.
  std::string sepStr = " \t";
  if (sopt.gencodeRef) {
    sepStr += '|';
  }

  std::vector<std::string> readFiles{fname_};
  size_t maxReadGroup{1000}; // Number of files to read simultaneously
  size_t concurrentFile{1};  // Number of reads in each "job"

  single_parser parser(readFiles, 1, 1, maxReadGroup);
  parser.start();

  constexpr char bases[] = {'A', 'C', 'G', 'T'};
  // Create a random uniform distribution
  constexpr const uint64_t randseed{271828}; 
  std::default_random_engine eng(randseed);
  std::uniform_int_distribution<> dis(0, 3);
  uint64_t numNucleotidesReplaced{0};

  // All header names we encounter in the fasta file
  std::unordered_set<std::string> fastaNames;

  auto rg = parser.getReadGroup();
  while (parser.refill(rg)) {
    for (auto& read : rg) {
      std::string& header = read.name;
      std::string name = header.substr(0, header.find_first_of(sepStr));

      if (fastaNames.find(name) != fastaNames.end()){
        sopt.jointLog->error("Transcript {} appears twice in the transcript FASTA file. "
                             "Duplicate transcripts or transcripts with the same name are not allowed.",
                             name);
        sopt.jointLog->flush();
        std::exit(1);
      }

      fastaNames.insert(name);

      auto it = nameToID.find(name);
      if (it == nameToID.end()) {
        sopt.jointLog->warn("Transcript {} appears in the reference but did "
                            "not appear in the BAM",
                            name);
      } else {

        std::string& seq = read.seq;
        size_t readLen = seq.length();

        refs[it->second].setSAMSequenceOwned(
            salmon::stringtools::encodeSequenceInSAM(seq.c_str(), readLen));

        // Replace non-ACGT bases
        for (size_t b = 0; b < readLen; ++b) {
          seq[b] = ::toupper(seq[b]);
          int c = jellyfish::mer_dna::code(seq[b]);
          // Replace non-ACGT bases with pseudo-random bases
          if (jellyfish::mer_dna::not_dna(c)) {
            char rbase = bases[dis(eng)];
            c = jellyfish::mer_dna::code(rbase);
            seq[b] = rbase;
            ++numNucleotidesReplaced;
          }
        }

        // allocate space for the new copy
        char* seqCopy = new char[seq.length() + 1];
        std::strcpy(seqCopy, seq.c_str());
        refs[it->second].setSequenceOwned(seqCopy, sopt.gcBiasCorrect,
                                          sopt.reduceGCMemory);
        // seqCopy will only be freed when the transcript is destructed!
      }
    }
  }

  parser.stop();

  // Check that every sequence present in the BAM header was also present in the
  // transcriptome fasta.
  bool missingTxpError{false};
  for (auto& kv : nameToID) {
    auto& name = kv.first;
    if (fastaNames.find(name) == fastaNames.end()) {
      sopt.jointLog->critical("Transcript {} appeared in the BAM header, but "
                              "was not in the provided FASTA file",
                              name);
      missingTxpError = true;
    }
  }

  if (missingTxpError) {
    sopt.jointLog->critical(
        "Please provide a reference FASTA file that includes all targets "
        "present in the BAM header\n"
        "If you have access to the genome FASTA and GTF used for alignment \n"
        "consider generating a transcriptome fasta using a command like: \n"
        "gffread -w output.fa -g genome.fa genome.gtf\n"
        "you can find the gffread utility at "
        "(http://ccb.jhu.edu/software/stringtie/gff.shtml)");
    sopt.jointLog->flush();
    std::exit(1);
  }

  sopt.jointLog->info(
      "replaced {:n} non-ACGT nucleotides with random nucleotides",
      numNucleotidesReplaced);
}
