#include <boost/thread/thread.hpp>

#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <atomic>
#include <chrono>
#include <thread>
#include <functional>
#include <memory>
#include <cassert>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>


#include "cereal/types/vector.hpp"
#include "cereal/archives/binary.hpp"

#include "jellyfish/config.h"
#include "jellyfish/err.hpp"
#include "jellyfish/misc.hpp"
#include "jellyfish/jellyfish.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"


#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>
#include <boost/filesystem.hpp>

#include "tbb/parallel_for_each.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_sort.h"

#include "Transcript.hpp"
#include "SalmonUtils.hpp"
#include "SalmonIndex.hpp"
#include "GenomicFeature.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/details/format.h"

using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

typedef jellyfish::stream_manager<char**>                stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> sequence_parser;

extern "C" {
int bwa_index(int argc, char* argv[]);
}

int computeBiasFeatures(
    std::vector<std::string>& transcriptFiles,
    boost::filesystem::path outFilePath,
    bool useStreamingParser,
    size_t numThreads);

// Cool way to do this from
// http://stackoverflow.com/questions/108318/whats-the-simplest-way-to-test-whether-a-number-is-a-power-of-2-in-c
bool isPowerOfTwo(uint32_t n) {
  return (n > 0 and (n & (n-1)) == 0);
}

int salmonIndex(int argc, char* argv[]) {

    using std::string;
    namespace bfs = boost::filesystem;
    namespace po = boost::program_options;

    bool useStreamingParser = true;

    uint32_t saSampInterval = 1;
    uint32_t auxKmerLen = 0;
    uint32_t maxThreads = std::thread::hardware_concurrency();
    uint32_t numThreads;
    bool useQuasi{false};

    po::options_description generic("Command Line Options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("transcripts,t", po::value<string>()->required(), "Transcript fasta file.")
    ("kmerLen,k", po::value<uint32_t>(&auxKmerLen)->default_value(31)->required(),
                    "The size of k-mers that should be used for the quasi index.")
    ("index,i", po::value<string>()->required(), "Salmon index.")
    ("threads,p", po::value<uint32_t>(&numThreads)->default_value(maxThreads)->required(),
                            "Number of threads to use (only used for computing bias features)")
    ("quasi,q", po::bool_switch(&useQuasi)->default_value(false), "Build quasi-mapping index instead of FMD-based index")
    ("sasamp,s", po::value<uint32_t>(&saSampInterval)->default_value(1)->required(),
                            "The interval at which the suffix array should be sampled. "
                            "Smaller values are faster, but produce a larger index. "
                            "The default should be OK, unless your transcriptome is huge. "
			    "This value should be a power of 2.")
    ;

    po::variables_map vm;
    int ret = 1;
    try {

        po::store(po::command_line_parser(argc, argv).options(generic).run(), vm);

        if ( vm.count("help") ) {
            auto hstring = R"(
Index
==========
Creates a salmon index.
)";
            std::cout << hstring << std::endl;
            std::cout << generic << std::endl;
            std::exit(1);
        }
        po::notify(vm);

        uint32_t sasamp = vm["sasamp"].as<uint32_t>();
        if (!isPowerOfTwo(sasamp) and !useQuasi) {
          fmt::MemoryWriter errWriter;
          errWriter << "Error: The suffix array sampling interval must be "
                   "a power of 2. The value provided, " << sasamp << ", is not.";
          throw(std::logic_error(errWriter.str()));
        }

        string transcriptFile = vm["transcripts"].as<string>();
        bfs::path indexDirectory(vm["index"].as<string>());


        if (!bfs::exists(indexDirectory)) {
            std::cerr << "index [" << indexDirectory << "] did not previously exist "
                      << " . . . creating it\n";
            bfs::create_directory(indexDirectory);
        }

        bfs::path logPath = indexDirectory / "indexing.log";
        size_t max_q_size = 2097152;
        //spdlog::set_async_mode(max_q_size);

        auto fileSink = std::make_shared<spdlog::sinks::simple_file_sink_mt>(logPath.string(), true);
        auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
        auto consoleLog = spdlog::create("consoleLog", {consoleSink});
        auto fileLog = spdlog::create("fileLog", {fileSink});
        auto jointLog = spdlog::create("jointLog", {fileSink, consoleSink});

        // First, compute the transcript features in case the user
        // ever wants to bias-correct his / her results
        // NOTE: Currently, we're using the same bias correction technique here that
        // we use in Sailfish. In the future, test more "traditional" bias correction
        // techniques to see if we should adopt them instead
        bfs::path transcriptBiasFile(indexDirectory); transcriptBiasFile /= "bias_feats.txt";

        std::vector<std::string> transcriptFiles = {transcriptFile};
        fmt::MemoryWriter infostr;
        infostr << "computeBiasFeatures( {";
        for (auto& tf : transcriptFiles) {
            infostr << "[" << tf << "] ";
        }
        infostr << ", " << transcriptBiasFile.c_str() << ", " << useStreamingParser << ", " << numThreads << ")\n";
        jointLog->info() << infostr.str();
        computeBiasFeatures(transcriptFiles, transcriptBiasFile, useStreamingParser, numThreads);
        // ==== finished computing bias fetures

        bfs::path outputPrefix;
        std::unique_ptr<std::vector<char const*>> argVec(new std::vector<char const*>);
	    fmt::MemoryWriter optWriter;

        std::unique_ptr<SalmonIndex> sidx = nullptr;
        // Build a quasi-mapping index
        if (useQuasi) {
            outputPrefix = indexDirectory;
            argVec->push_back("dummy");
            argVec->push_back("-k");

            if (auxKmerLen == 0) {
                jointLog->info("You cannot have a k-mer length of 0 with the quasi-index.");
                jointLog->info("Setting to the default value of 31.");
                auxKmerLen = 31;
            }

            optWriter << auxKmerLen;
            argVec->push_back(optWriter.str().c_str());
            argVec->push_back("-t");
            argVec->push_back(transcriptFile.c_str());
            argVec->push_back("-i");
            argVec->push_back(outputPrefix.string().c_str());
            sidx.reset(new SalmonIndex(jointLog, IndexType::QUASI));
        } else {
            // Build the FMD-based index
            bfs::path outputPrefix = indexDirectory / "bwaidx";
            argVec->push_back("index");
            argVec->push_back("-s");
	        optWriter << vm["sasamp"].as<uint32_t>();
            argVec->push_back(optWriter.str().c_str());
            argVec->push_back("-p");
            argVec->push_back(outputPrefix.string().c_str());
            argVec->push_back(transcriptFile.c_str());
            sidx.reset(new SalmonIndex(jointLog, IndexType::FMD));
	    // Disable the auxiliary k-mer index for now
	    auxKmerLen = 0;
        }

        jointLog->info("building index");
	    sidx->build(indexDirectory, *(argVec.get()), auxKmerLen);
        jointLog->info("done building index");
        // If we want to build the auxiliary k-mer index, do it here.
        /*
        uint32_t k = 15;
        buildAuxKmerIndex(outputPrefix, k, jointLog);
        */

    } catch (po::error &e) {
        std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (const spdlog::spdlog_ex& ex) {
        std::cout << "logger failed with : [" << ex.what() << "]. Exiting.\n";
    } catch (std::exception& e) {
        std::cerr << "Exception : [" << e.what() << "]\n";
        std::cerr << argv[0] << " index was invoked improperly.\n";
        std::cerr << "For usage information, try " << argv[0] << " index --help\nExiting.\n";
    }
    return ret;
}

