#include <boost/thread/thread.hpp>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <thread>
#include <vector>

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/range/irange.hpp>

#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_sort.h"

#include "GenomicFeature.hpp"
#include "SalmonIndex.hpp"
#include "SalmonUtils.hpp"
#include "Transcript.hpp"
#include "pufferfish/ProgOpts.hpp"
#include "spdlog/fmt/fmt.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/spdlog.h"

// Cool way to do this from
// http://stackoverflow.com/questions/108318/whats-the-simplest-way-to-test-whether-a-number-is-a-power-of-2-in-c
bool isPowerOfTwo(uint32_t n) { return (n > 0 and (n & (n - 1)) == 0); }

int salmonIndex(int argc, const char* argv[]) {

  using std::string;
  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;

  bool useStreamingParser = true;

  string indexTypeStr = "quasi";
  string decoyList;
  uint32_t saSampInterval = 1;
  uint32_t auxKmerLen = 0;
  uint32_t numThreads{2};
  bool useQuasi{false};
  bool perfectHash{false};
  bool gencodeRef{false};
  bool featuresRef{false};
  bool keepDuplicates{false};

  pufferfish::IndexOptions idxOpt;

  po::options_description generic("Command Line Options");
  generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("transcripts,t", po::value<string>()->required(), "Transcript fasta file.")
    ("kmerLen,k",
      po::value<uint32_t>(&idxOpt.k)->default_value(31)->required(),
      "The size of k-mers that should be used for the quasi index.")
    ("index,i", po::value<string>()->required(), "salmon index.")
    ("gencode", po::bool_switch(&gencodeRef)->default_value(false),
      "This flag will expect the input transcript fasta to be in GENCODE "
      "format, and will split the transcript name at the first \'|\' character. "
      "These reduced names will be used in the output and when looking for these "
      "transcripts in a gene to transcript "
      "GTF.")
    ("features", po::bool_switch(&featuresRef)->default_value(false),
     "This flag will expect the input reference to be in the tsv file "
     "format, and will split the feature name at the first \'tab\' character. "
     "These reduced names will be used in the output and when looking for the "
     "sequence of the features."
     "GTF.")
    ("keepDuplicates",po::bool_switch(&idxOpt.keep_duplicates)->default_value(false),
      "This flag will disable the default indexing behavior of "
       "discarding sequence-identical duplicate "
       "transcripts.  If this flag is passed, then duplicate "
       "transcripts that appear in the input will be "
       "retained and quantified separately.")
    ("threads,p", po::value<uint32_t>(&idxOpt.p)->default_value(2)->required(),
      "Number of threads to use during indexing.")
    ("keepFixedFasta",po::bool_switch(&idxOpt.keep_fixed_fasta)->default_value(false),
      "Retain the fixed fasta file (without short transcripts and duplicates, "
      "clipped, etc.) generated during indexing")
    ("filterSize,f", po::value<int32_t>(&idxOpt.filt_size)->default_value(-1),
      "The size of the Bloom filter that will be used by TwoPaCo during indexing. "
      "The filter will be of size 2^{filterSize}. The default value of -1 "
      "means that the filter size will be automatically set based on the number of "
      "distinct k-mers in the input, as estimated by nthll.")
    ("tmpdir",po::value<std::string>(&idxOpt.twopaco_tmp_dir)->default_value(""),
      "The directory location that will be used for TwoPaCo temporary files; "
      "it will be created if need be and be removed prior to indexing completion. "
      "The default value will cause a (temporary) subdirectory of the salmon index "
      "directory to be used for this purpose.")
    ("sparse", po::bool_switch(&idxOpt.isSparse)->default_value(false),
      "Build the index using a sparse sampling of k-mer positions "
      "This will require less memory (especially during quantification), but "
      "will take longer to construct and can slow down mapping / alignment")
    ("decoys,d", po::value<string>(&idxOpt.decoy_file),
      "Treat these sequences ids from the reference as the decoys that may "
      "have sequence homologous to some known transcript. for example in "
      "case of the genome, provide a list of chromosome name --- one per line")
    ("type",
      po::value<string>(&indexTypeStr)->default_value("puff")->required(),
      "The type of index to build; the only option is \"puff\" in this version "
      "of salmon.");

  po::variables_map vm;
  int ret = 0;
  try {

    po::store(po::command_line_parser(argc, argv).options(generic).run(), vm);

    if (vm.count("help")) {
      auto hstring = R"(
Index
==========
Creates a salmon index.
)";
      std::cerr << hstring << std::endl;
      std::cerr << generic << std::endl;
      std::exit(0);
    }
    po::notify(vm);

    if (indexTypeStr == "fmd") {
      fmt::MemoryWriter errWriter;
      errWriter
          << "Error: FMD indexing is not supported in this version of salmon.";
      throw(std::logic_error(errWriter.str()));
    }
    if (indexTypeStr == "quasi") {
      fmt::MemoryWriter errWriter;
      errWriter << "Error: RapMap-based indexing is not supported in this "
                   "version of salmon.";
      throw(std::logic_error(errWriter.str()));
    }
    if (indexTypeStr != "puff") {
      fmt::MemoryWriter errWriter;
      errWriter
          << "Error: If explicitly provided, the index type must be \"puff\"."
          << "You passed [" << indexTypeStr << "], but this is not supported.";
      throw(std::logic_error(errWriter.str()));
    }

    bool usePuff = (indexTypeStr == "puff");
    string transcriptFile = vm["transcripts"].as<string>();
    bfs::path indexDirectory(vm["index"].as<string>());

    // Check that the transcriptome file exists
    if (!bfs::exists(transcriptFile)) {
      std::cerr << "The file [" << transcriptFile
                << "] provided for the transcriptome "
                << "does not appear to exist.";
      std::exit(1);
    }
    // and is not a directory
    if (bfs::is_directory(transcriptFile)) {
      std::cerr << "The provided transcriptome argument [" << transcriptFile
                << "] appears to be a directory; "
                << "please provide a file.";
      std::exit(1);
    }

    if (!bfs::exists(indexDirectory)) {
      std::cerr << "index [" << indexDirectory << "] did not previously exist "
                << " . . . creating it\n";
      bfs::create_directories(indexDirectory);
    }

    bfs::path logPath = indexDirectory / "pre_indexing.log";
    size_t max_q_size = 2097152;
    // spdlog::set_async_mode(max_q_size);

    auto fileSink = std::make_shared<spdlog::sinks::simple_file_sink_mt>(
        logPath.string(), true);
    auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
    auto consoleLog = spdlog::create("consoleLog", {consoleSink});
    auto fileLog = spdlog::create("fLog", {fileSink});
    auto jointLog = spdlog::create("jLog", {fileSink, consoleSink});

    idxOpt.rfile = {transcriptFile};
    fmt::MemoryWriter infostr;

    std::unique_ptr<SalmonIndex> sidx = nullptr;
    // Build a quasi-mapping index
    if (usePuff) {
      idxOpt.outdir = indexDirectory.string();
      if (idxOpt.k == 0) {
        jointLog->info(
            "You cannot have a k-mer length of 0 with the pufferfish index.");
        jointLog->info("Setting to the default value of 31.");
        idxOpt.k = 31;
      }

      // give the user a warning if they are not using any decoy file
      if (idxOpt.decoy_file.empty()) {
        jointLog->warn("The salmon index is being built without any decoy sequences.  It is recommended that "
                      "decoy sequence (either computed auxiliary decoy sequence or the genome of the organism) "
                      "be provided during indexing. Further details can be found at "
                      "https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode.");
      }

      if (gencodeRef) {
        idxOpt.header_sep = "|";
      }

      if (featuresRef) { idxOpt.featuresRef = true; }
      
      // by default we expect to be indexing transcriptome references 
      // so set the expect_transcriptome flag 
      idxOpt.expect_transcriptome = true;

      sidx.reset(new SalmonIndex(jointLog, SalmonIndexType::PUFF));
    } else {
      jointLog->error("This version of salmon does not support FMD or "
                      "RapMap-based indexing.");
      return 1;
    }

    jointLog->info("building index");
    sidx->build(indexDirectory, idxOpt);
    jointLog->info("done building index");
    // If we want to build the auxiliary k-mer index, do it here.
    /*
    uint32_t k = 15;
    buildAuxKmerIndex(outputPrefix, k, jointLog);
    */

  } catch (po::error& e) {
    std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (const spdlog::spdlog_ex& ex) {
    std::cerr << "logger failed with : [" << ex.what() << "]. Exiting.\n";
    ret = 1;
  } catch (std::exception& e) {
    std::cerr << "Exception : [" << e.what() << "]\n";
    std::cerr << argv[0] << " index was invoked improperly.\n";
    std::cerr << "For usage information, try " << argv[0]
              << " index --help\nExiting.\n";
    ret = 1;
  }
  return ret;
}
