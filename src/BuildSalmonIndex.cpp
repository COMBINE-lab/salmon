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
#include "tbb/task_scheduler_init.h"

#include "GenomicFeature.hpp"
#include "SalmonIndex.hpp"
#include "SalmonUtils.hpp"
#include "Transcript.hpp"
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
  bool keepDuplicates{false};

  po::options_description generic("Command Line Options");
  generic.add_options()("version,v", "print version string")(
      "help,h", "produce help message")("transcripts,t",
                                        po::value<string>()->required(),
                                        "Transcript fasta file.")(
      "kmerLen,k",
      po::value<uint32_t>(&auxKmerLen)->default_value(31)->required(),
      "The size of k-mers that should be used for the quasi index.")(
      "index,i", po::value<string>()->required(), "salmon index.")(
      "gencode", po::bool_switch(&gencodeRef)->default_value(false),
      "This flag will expect the input transcript fasta to be in GENCODE "
      "format, and will split "
      "the transcript name at the first \'|\' character.  These reduced names "
      "will be used in the "
      "output and when looking for these transcripts in a gene to transcript "
      "GTF.")("keepDuplicates",
              po::bool_switch(&keepDuplicates)->default_value(false),
              "This flag will disable the default indexing behavior of "
              "discarding sequence-identical duplicate "
              "transcripts.  If this flag is passed, then duplicate "
              "transcripts that appear in the input will be "
              "retained and quantified separately.")(
      "threads,p",
      po::value<uint32_t>(&numThreads)->default_value(2)->required(),
      "Number of threads to use (only used for computing bias features)")(
      "perfectHash", po::bool_switch(&perfectHash)->default_value(false),
      "[quasi index only] Build the index using a perfect hash rather than a "
      "dense hash.  This "
      "will require less memory (especially during quantification), but will "
      "take longer to construct")(
      "decoys,d",
      po::value<string>(&decoyList),
      "Treat these sequences as decoys that may have sequence homologous to some "
      "known transcript"
      )(
      "type",
      po::value<string>(&indexTypeStr)->default_value("quasi")->required(),
      "The type of index to build; the only option is \"quasi\" in this version of salmon.");

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
      errWriter << "Error: FMD indexing is not supported in this version of salmon.";
      throw(std::logic_error(errWriter.str()));
    }
    if (indexTypeStr != "quasi") {
      fmt::MemoryWriter errWriter;
      errWriter << "Error: If explicitly provided, the index type must be \"quasi\"."
                << "You passed [" << indexTypeStr << "], but this is not supported.";
      throw(std::logic_error(errWriter.str()));
    }

    bool useQuasi = (indexTypeStr == "quasi");
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

    bool haveDecoys{false};
    if (vm.count("decoys")) {
      haveDecoys = true;
    }

    bfs::path logPath = indexDirectory / "indexing.log";
    size_t max_q_size = 2097152;
    // spdlog::set_async_mode(max_q_size);

    auto fileSink = std::make_shared<spdlog::sinks::simple_file_sink_mt>(
        logPath.string(), true);
    auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
    auto consoleLog = spdlog::create("consoleLog", {consoleSink});
    auto fileLog = spdlog::create("fLog", {fileSink});
    auto jointLog = spdlog::create("jLog", {fileSink, consoleSink});

    std::vector<std::string> transcriptFiles = {transcriptFile};
    fmt::MemoryWriter infostr;

    bfs::path outputPrefix;
    std::unique_ptr<std::vector<std::string>> argVec(new std::vector<std::string>);
    fmt::MemoryWriter optWriter;

    std::unique_ptr<SalmonIndex> sidx = nullptr;
    // Build a quasi-mapping index
    if (useQuasi) {
      outputPrefix = indexDirectory;
      argVec->push_back("dummy");
      argVec->push_back("-k");

      if (auxKmerLen == 0) {
        jointLog->info(
            "You cannot have a k-mer length of 0 with the quasi-index.");
        jointLog->info("Setting to the default value of 31.");
        auxKmerLen = 31;
      }

      optWriter << auxKmerLen;
      argVec->push_back(optWriter.str());
      optWriter.clear();

      argVec->push_back("-t");
      argVec->push_back(transcriptFile);
      argVec->push_back("-i");
      argVec->push_back(outputPrefix.string());

      argVec->push_back("-x");
      optWriter << numThreads;
      argVec->push_back(optWriter.str());
      optWriter.clear();

      if (perfectHash) {
        argVec->push_back("--perfectHash");
      }
      if (gencodeRef) {
        argVec->push_back("-s");
        argVec->push_back("\"|\"");
      }
      if (keepDuplicates) {
        argVec->push_back("--keepDuplicates");
      }
      if (haveDecoys) {
        argVec->push_back("--decoys");
        argVec->push_back(decoyList);
      }

      sidx.reset(new SalmonIndex(jointLog, SalmonIndexType::QUASI));
    } else {
      jointLog->error("This version of salmon does not support FMD indexing.");
      return 1;
    }

    jointLog->info("building index");
    sidx->build(indexDirectory, *(argVec.get()), auxKmerLen);
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
