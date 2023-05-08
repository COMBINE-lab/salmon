//
// Pufferfish - An efficient compacted dBG index
//
// Copyright (C) 2017 Rob Patro, Fatemeh Almodaresi, Hirak Sarkar
//
// This file is part of pufferfish.
//
// pufferfish is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pufferfish is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pufferfish.  If not, see <http://www.gnu.org/licenses/>.
//

#include "clipp.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <clocale>
#include <ghc/filesystem.hpp>
#include <cereal/archives/json.hpp>

#include "PufferfishConfig.hpp"
#include "ProgOpts.hpp"
#include "Util.hpp"
//#include "IndexHeader.hpp"

int pufferfishIndex(pufferfish::IndexOptions& indexOpts); // int argc, char* argv[]);
int pufferfishTest(pufferfish::TestOptions& testOpts);    // int argc, char* argv[]);
int pufferfishValidate(pufferfish::ValidateOptions& validateOpts); // int argc, char* argv[]);
int pufferfishTestLookup(pufferfish::ValidateOptions& lookupOpts); // int argc, char* argv[]);
int pufferfishAligner(pufferfish::AlignmentOpts& alignmentOpts) ;
int pufferfishExamine(pufferfish::ExamineOptions& examineOpts);
int pufferfishStats(pufferfish::StatsOptions& statsOpts);
int pufferfishKmerQuery(pufferfish::KmerQueryOptions& kqueryOpts);

int main(int argc, char* argv[]) {
  using namespace clipp;
  using std::cout;
  std::setlocale(LC_ALL, "en_US.UTF-8");


  enum class mode {help, index, validate, lookup, kquery, align, examine, stat};
  mode selected = mode::help;

  std::map<std::string, clipp::parameter> cmd_map = {
    {"align", command("align").set(selected, mode::align)},
    {"index", command("index").set(selected, mode::index)},
    {"validate", command("validate").set(selected, mode::validate)},
    {"lookup", command("lookup").set(selected, mode::lookup)},
    {"kquery", command("kquery").set(selected, mode::kquery)},
    {"examine", command("examine").set(selected, mode::examine)},
    {"stat", command("stat").set(selected, mode::stat)}
  };

  pufferfish::AlignmentOpts alignmentOpt ;
  pufferfish::IndexOptions indexOpt;
  //TestOptions testOpt;
  pufferfish::ValidateOptions validateOpt;
  pufferfish::ValidateOptions lookupOpt;
  pufferfish::KmerQueryOptions kmerQueryOpt;
  pufferfish::ExamineOptions examineOpt;
  pufferfish::StatsOptions statOpt;

  auto ensure_file_exists = [](const std::string& s) -> bool {
      bool exists = ghc::filesystem::exists(s);
      if (!exists) {
        std::string e = "The required input file " + s + " does not seem to exist.";
        throw std::runtime_error{e};
      }
      return true;
  };

  auto ensure_index_exists = [](const std::string& s) -> bool {
      bool exists = ghc::filesystem::exists(s);
      if (!exists) {
        std::string e = "The index directory " + s + " does not seem to exist.";
        throw std::runtime_error{e};
      }
      bool isDir = ghc::filesystem::is_directory(s);
      if (!isDir) {
          std::string e = s + " is not a directory containing index files.";
          throw std::runtime_error{e};
      }
  
	  for (auto & elem : {pufferfish::util::INFO,
                          pufferfish::util::MPH,
                          pufferfish::util::SEQ,
                          pufferfish::util::RANK,
                          pufferfish::util::CTABLE,
                          pufferfish::util::REFSEQ,
                          pufferfish::util::REFNAME,
                          pufferfish::util::REFLENGTH,
                          pufferfish::util::REFACCUMLENGTH}) {
          if (!ghc::filesystem::exists(s+"/"+elem)) {
              std::string e = "Index is incomplete. Missing file ";
              e+=elem;
              throw std::runtime_error{e};
          }
      }
      std::string indexType;
      {
          std::ifstream infoStream(s + "/" + pufferfish::util::INFO);
          cereal::JSONInputArchive infoArchive(infoStream);
          infoArchive(cereal::make_nvp("sampling_type", indexType));
          infoStream.close();
          if (indexType == "dense") {
              if (!ghc::filesystem::exists(s+"/"+pufferfish::util::POS)) {
                  std::string e = "Index is incomplete. Missing file " + std::string(pufferfish::util::POS);
                  throw std::runtime_error{e};
              }
          } else if (indexType == "sparse") {
              for (auto & elem : {pufferfish::util::EXTENSION,
                                  pufferfish::util::EXTENSIONSIZE,
                                  pufferfish::util::SAMPLEPOS}) {
                  if (!ghc::filesystem::exists(s+"/"+elem)) {
                      std::string e = "Index is incomplete. Missing file ";
                      e+=elem;
                      throw std::runtime_error{e};
                  }
              }
          }
      }
 
      return true;
  };



  auto indexMode = (
                    command("index").set(selected, mode::index),
                    (required("-r", "--ref") & values(ensure_file_exists, "ref_file", indexOpt.rfile)) % "path to the reference fasta file",
                    (required("-o", "--output") & value("output_dir", indexOpt.outdir)) % "directory where index is written",
                    //(required("-g", "--gfa") & value("gfa_file", indexOpt.gfa_file)) % "path to the GFA file",
                    (option("--expectTranscriptome").set(indexOpt.expect_transcriptome) % 
                    "expect (non-decoy) sequences to be transcripts rather than genomic contigs"),
                    (option("--headerSep") & value("sep_strs", indexOpt.header_sep)) %
                    "Instead of a space or tab, break the header at the first "
                    "occurrence of this string, and name the transcript as the token before "
                    "the first separator (default = space & tab)",
                    (option("--keepFixedFasta").set(indexOpt.keep_fixed_fasta) % "Retain the fixed fasta file (without short transcripts and duplicates, clipped, etc.) generated during indexing"),
                    (option("--keepDuplicates").set(indexOpt.keep_duplicates) % "Retain duplicate references in the input"),
                    (option("-d", "--decoys") & value("decoy_list", indexOpt.decoy_file)) %
                    "Treat these sequences as decoys that may be sequence-similar to some known indexed reference",
                    (option("-n", "--noClip").set(indexOpt.noclip_polya)) % "Don't clip poly-A tails from the ends of target sequences",
                    (option("-f", "--filt-size") & value("filt_size", indexOpt.filt_size)) % "filter size to pass to TwoPaCo when building the reference dBG",
                    (option("--tmpdir") & value("twopaco_tmp_dir", indexOpt.twopaco_tmp_dir)) % "temporary work directory to pass to TwoPaCo when building the reference dBG",
                    (option("-k", "--klen") & value("kmer_length", indexOpt.k))  % "length of the k-mer with which the dBG was built (default = 31)",
                    (option("-p", "--threads") & value("threads", indexOpt.p))  % "total number of threads to use for building MPHF (default = 16)",
                    (option("-l", "--build-edges").set(indexOpt.buildEdgeVec, true) % "build and record explicit edge table for the contaigs of the ccdBG (default = false)"),
                    (option("-q", "--build-eqclses").set(indexOpt.buildEqCls, true) % "build and record equivalence classes (default = false)"),
                    (((option("-s", "--sparse").set(indexOpt.isSparse, true)) % "use the sparse pufferfish index (less space, but slower lookup)",
                     ((option("-e", "--extension") & value("extension_size", indexOpt.extensionSize)) % "length of the extension to store in the sparse index (default = 4)")) |
                     ((option("-x", "--lossy-rate").set(indexOpt.lossySampling, true)) & value("lossy_rate", indexOpt.lossy_rate) % "use the lossy sampling index with a sampling rate of x (less space and fast, but lower sensitivity)"))
                    );

  // Examine properties of the index
  auto examineMode = (
                      command("examine").set(selected, mode::examine),
                      (required("-i", "--index") & value("index", examineOpt.index_dir)) % "pufferfish index directory",
                      (option("--dump-fasta") & value("fasta_out", examineOpt.fasta_out)) %
                      "dump the reference sequences in the index in the provided fasta file",
                      (option("--dump-kmer-freq") & value("kmer_freq_out", examineOpt.kmer_freq_out)) %
                      "dump the frequency histogram of k-mers"
                      );
  /*
  auto testMode = (
                   command("test").set(selected, mode::test)
                   );
  */
  auto validateMode = (
                       command("validate").set(selected, mode::validate),
                       (required("-i", "--index") & value("index", validateOpt.indexDir)) % "directory where the pufferfish index is stored");/*,
                       (required("-r", "--ref") & value("ref", validateOpt.refFile)) % "fasta file with reference sequences",
                       (required("-g", "--gfa") & value("gfa", validateOpt.gfaFileName)) % "GFA file name needed for edge table validation"
                       );
                                                                                                                                              */
  auto lookupMode = (
                     command("lookup").set(selected, mode::lookup),
                     (required("-i", "--index") & value("index", lookupOpt.indexDir)) % "directory where the pufferfish index is stored",
                     (required("-r", "--ref") & value("ref", lookupOpt.refFile)) % "fasta file with reference sequences"
                     );
  
  auto kmerQueryMode = (
                      command("kquery").set(selected, mode::kquery),
                     (required("-i", "--index") & value("index", kmerQueryOpt.indexDir)) % "directory where the pufferfish index is stored",
                     (required("-q", "--query") & values("ref", kmerQueryOpt.queryFiles)) % "fasta file with reference sequences",
                     (option("-p", "--threads") & value("threads", kmerQueryOpt.num_threads))  % "total number of threads to use for mapping k-mers"
                      );

  std::string statType = "ctab";
  auto statMode = (
                    command("stat").set(selected, mode::stat),
                    (option("-t", "--type") & value("statType", statType)) % "statType (options:ctab, motif)",
                    (required("-i", "--index") & value("index", statOpt.indexDir)) % "directory where the pufferfish index is stored");
  if (statType == "ctab") {
      //std::cerr << statType << "\n";
      statOpt.statType = pufferfish::StatType::ctab;
  } else if (statType == "motif") {
      //std::cerr << statType << "\n";
      statOpt.statType = pufferfish::StatType::motif;
  }

  std::string throwaway;
  auto isValidRatio = [](const char* s) -> void {
    float r{0.0};
    std::string sv(s);
    try {
      r = std::stof(sv);
    } catch (std::exception& e) {
      std::string m = "Could not convert " + sv + " to a valid ratio\n";
      throw std::domain_error(m);
    }
    if (r <= 0 or r > 1) {
      std::string m = "The --scoreRatio you provided was " + sv + ", it must be in (0,1]\n";
      throw std::domain_error(m);
    }
  };

  auto alignMode = (
                    command("align").set(selected, mode::align),
                    (required("-i", "--index") & value(ensure_index_exists, "index", alignmentOpt.indexDir)) % "Directory where the Pufferfish index is stored",
                    (
                      (
                        ((required("--mate1", "-1") & value("mate 1", alignmentOpt.read1)) % "Path to the left end of the read files"),
                        ((required("--mate2", "-2") & value("mate 2", alignmentOpt.read2)) % "Path to the right end of the read files")
                      ) 
                      |
                      ((required("--read").set(alignmentOpt.singleEnd, true) & value("reads", alignmentOpt.unmatedReads)) % "Path to single-end read files")
                    ),
                    (option("-b", "--batchOfReads").set(alignmentOpt.listOfReads, true)) % "Is each input a file containing the list of reads? (default=false)",


                    (option("--coverageScoreRatio") & value("score ratio", alignmentOpt.scoreRatio).call(isValidRatio)) % "Discard mappings with a coverage score < scoreRatio * OPT (default=0.6)",
                    (option("-t", "--threads") & value("num threads", alignmentOpt.numThreads)) % "Specify the number of threads (default=8)",
                    (option("-m", "--just-mapping").set(alignmentOpt.justMap, true)) % "don't attempt alignment validation; just do mapping",
                    (
                      (required("--noOutput").set(alignmentOpt.noOutput, true)) % "Run without writing SAM file"
                      |
                      (required("-o", "--outdir") & value("output file", alignmentOpt.outname)) % "Output file where the alignment results will be stored"
                    ),
                    (option("--allowSoftclip").set(alignmentOpt.allowSoftclip, true) % "Allow soft-clipping at start and end of alignments"),
                    (option("--allowOverhangSoftclip").set(alignmentOpt.allowOverhangSoftclip, true) % "Allow soft-clipping part of a read that overhangs the reference (the regular --allowSoftclip flag overrides this one)"),
                    (option("--maxSpliceGap") & value("max splice gap", alignmentOpt.maxSpliceGap)) % "Specify maximum splice gap that two uni-MEMs should have",
                    (option("--maxFragmentLength") & value("max frag length", alignmentOpt.maxFragmentLength)) % 
                            "Specify the maximum distance between the last uni-MEM of the left and first uni-MEM of the right end of the read pairs (default:1000)",
                    (option("--noOrphans").set(alignmentOpt.noOrphan, true)) % "Write Orphans flag",
                    (option("--writeQualities").set(alignmentOpt.writeQualities, true)) % "Write quality strings in output SAM",
                    (option("--orphanRecovery").set(alignmentOpt.recoverOrphans, true)) % "Recover mappings for the other end of orphans using alignment",
                    (option("--noDiscordant").set(alignmentOpt.noDiscordant, true)) % "Write Orphans flag",
                    (option("--noDovetail").set(alignmentOpt.noDovetail, true)) % "Disallow dovetail alignment for paired end reads",
		            (option("-z", "--compressedOutput").set(alignmentOpt.compressedOutput, true)) % "Compress (gzip) the output file",
                    (
                      (option("-r", "--radOut").set(alignmentOpt.radOut, true)) % "Write output in the format required for krakMap"
                      |
                      (option("-k", "--krakOut").set(alignmentOpt.krakOut, true)) % "Write output in the format required for krakMap"
                      |
                      (option("-p", "--pam").set(alignmentOpt.salmonOut, true)) % "Write output in the format required for salmon"
                    ),
					(option("--verbose").set(alignmentOpt.verbose, true)) % "Print out auxilary information to trace program's flow",
                    (option("--fullAlignment").set(alignmentOpt.fullAlignment, true)) % "Perform full alignment instead of gapped alignment",
                    (option("--heuristicChaining").set(alignmentOpt.heuristicChaining, true)) % "Whether or not perform only 2 rounds of chaining",
                    (option("--bestStrata").set(alignmentOpt.bestStrata, true)) % "Keep only the alignments with the best score for each read",
					(option("--genomicReads").set(alignmentOpt.genomicReads, true)) % "Align genomic dna-seq reads instead of RNA-seq reads",
					(option("--primaryAlignment").set(alignmentOpt.primaryAlignment, true).set(alignmentOpt.bestStrata, true)) % "Report at most one alignment per read",
                    (option("--filterGenomics").set(alignmentOpt.filterGenomics, true) & value("genes names file", alignmentOpt.genesNamesFile)) % 
                         "Path to the file containing gene IDs. Filters alignments to the IDs listed in the file. Used to filter genomic reads while aligning to both genome and transcriptome."
                         "A read will be reported with only the valid gene ID alignments and will be discarded if the best alignment is to an invalid ID"
                         "The IDs are the same as the IDs in the fasta file provided for the index construction phase",
                    (option("--filterBestScoreMicrobiome").set(alignmentOpt.filterMicrobiomBestScore, true) & value("genes ID file", alignmentOpt.genesNamesFile)) % "Path to the file containing gene IDs. Same as option \"filterGenomics\" except that a read will be discarded if aligned equally best to a valid and invalid gene ID.",
                    (option("--filterMicrobiome").set(alignmentOpt.filterMicrobiom, true) & value("genes ID file", alignmentOpt.rrnaFile)) % "Path to the file containing gene IDs. Same as option \"filterGenomics\" except that a read will be discarded if an invalid gene ID is in the list of alignments.",
                    (option("--bt2DefaultThreshold").set(alignmentOpt.mimicBt2Default, true)) % "mimic the default threshold function of Bowtie2 which is t = -0.6 -0.6 * read_len",
                    (option("--minScoreFraction") & value("minScoreFraction", alignmentOpt.minScoreFraction)) % "Discard alignments with alignment score < minScoreFraction * max_alignment_score for that read (default=0.65)",
                    (option("--consensusFraction") & value("consensus fraction", alignmentOpt.consensusFraction)) % "The fraction of mems, relative to the reference with "
                    "the maximum number of mems, that a reference must contain in order "
                    "to move forward with computing an optimal chain score (default=0.65)",
                    (option("--gapOpenPenalty") & value("gap open penalty", alignmentOpt.gapOpenPenalty)) % " gap open penalty (default 5)",
                    (option("--gapExtendPenalty") & value("gap extend penalty", alignmentOpt.gapExtendPenalty)) % "gap extend penalty (default 3)",
                    (option("--mismatchScore") & value("gap extend score", alignmentOpt.mismatchScore)) % " mismatch score, should be negative or at least smaller than match score to make sense (default -4)",
                    (option("--matchScore") & value("match score", alignmentOpt.matchScore)) % " match score (default 2)",
                    (option("--altSkip") & value("alternative k-mer skip", alignmentOpt.altSkip)) % "Set the value of k-mer skipping, skipping happens if a mis-match is encountered at the time of querying k-kmers (default 5)",
                    (option("--preMergeChainSubThresh") & value("pre merge chain threshold", alignmentOpt.preMergeChainSubThresh)) % "The threshold of sub-optimal chains, compared to the best chain on a giventarget, that will be retained and passed to the next phase of mapping.  Specifically, if the best chain for a read (or read-end in paired-end mode) to target t has score X_t, then all chains for this read wit score >= X_t * preMergeChainSubThresh will be retained and passed to subsequent mapping phases.  This value must be in the range [0, 1] (default 0.9)",
                    (option("--postMergeChainSubThresh") & value("post merge chain threshold", alignmentOpt.postMergeChainSubThresh)) % "The threshold of sub-optimal chain pairs, compared to the best chain pair on a given target, that will be retained and passed to the next phase of mapping.  This is different than preMergeChainSubThresh, because this is applied to pairs of chains (from the ends of paired-end reads) after merging (i.e. after checking concordancy constraints etc.).  Specifically, if the best chain pair to target t has score X_t, then all chain pairs for this read pair with score >= X_t * postMergeChainSubThresh will be retained and passed to subsequent mapping phases.  This value must be in the range [0, 1]. Note: This option is only meaningful for paired-end libraries, and is ignored for single-end libraries. (default 0.9)",
                    (option("--orphanChainSubThresh") & value("orphan chain threshold", alignmentOpt.orphanChainSubThresh)) % "This threshold sets a global sub-optimality threshold for chains corresponding to orphan mappings.  That is, if the merging procedure results in no concordant mappings then only orphan mappings with a chain score >= orphanChainSubThresh * bestChainScore will be retained and passed to subsequent mapping phases.  This value must be in the range [0, 1]. Note: This option is only meaningful for paired-end libraries, and is ignored for single-end libraries. (default 1)",
                    (option("--noAlignmentCache").set(alignmentOpt.useAlignmentCache, false)) % "Do not use the alignment cache during the alignment."
  );

  auto cli = (
              (indexMode | validateMode | lookupMode | kmerQueryMode | alignMode | examineMode | statMode | command("help").set(selected,mode::help) ),
              option("-v", "--version").call([]{std::cout << "version " << pufferfish::version << "\n"; std::exit(0);}).doc("show version"));

  decltype(parse(argc, argv, cli)) res;
  try {
    res = parse(argc, argv, cli);
  } catch (std::exception& e) {
    std::cout << "\n\nParsing command line failed with exception: " << e.what() << "\n";
    std::cout << "\n\n";
    std::cout << make_man_page(cli, pufferfish::progname);
    return 1;
  }


  if(res) {
    switch(selected) {
    case mode::index: pufferfishIndex(indexOpt);  break;
    case mode::validate: pufferfishValidate(validateOpt);  break;
    case mode::lookup: pufferfishTestLookup(lookupOpt); break;
    case mode::kquery: pufferfishKmerQuery(kmerQueryOpt); break;
    case mode::align: pufferfishAligner(alignmentOpt); break;
    case mode::examine: pufferfishExamine(examineOpt); break;
    case mode::stat: pufferfishStats(statOpt); break;
    case mode::help: std::cout << make_man_page(cli, pufferfish::progname); break;
    }
  } else {
    std::string nl = "The valid commands to pufferfish are : ";
    bool first = true;
    for (auto& kv : cmd_map) {
      nl += (first ? "{" : " ") + kv.first + ",";
      first = false;
    }
    nl.pop_back();
    nl += "}";

    auto b = res.begin();
    auto e = res.end();
    // if there was a command provided
    if (std::distance(b,e) > 0) {
      if (b->arg() == "index") {
        std::cout << make_man_page(indexMode, pufferfish::progname);
      } else if (b->arg() == "validate") {
        std::cout << make_man_page(validateMode, pufferfish::progname);
      } else if (b->arg() == "lookup") {
        std::cout << make_man_page(lookupMode, pufferfish::progname);
      } else if (b->arg() == "align") {
        std::cout << make_man_page(alignMode, pufferfish::progname);
      } else {
        std::cout << "There is no command \"" << b->arg() << "\"\n";
        std::cout << nl << '\n';
        return 1;
      }
    } else {
      std::cout << nl << '\n';
      return 1;
    }
  }

  return 0;
}
