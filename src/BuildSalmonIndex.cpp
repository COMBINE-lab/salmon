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
#include "GenomicFeature.hpp"
#include "format.h"
#include "spdlog/spdlog.h"

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

/*
bool buildAuxKmerIndex(boost::filesystem::path& outputPrefix, uin32_t k,
                       std::shared_ptr<spdlog::logger>& logger){
    namespace bfs = boost::filesystem;

    bfs::path indexPath = outputPrefix / "bwaindex";
    // Load the bwa index
    bwaidx_t *idx_{nullptr};
    {
        logger->info("Reading BWT index from file");
        if ((idx_ = bwa_idx_load(indexPath.string().c_str(), BWA_IDX_BWT|BWA_IDX_BNS|BWA_IDX_PAC)) == 0) {
            logger->error("Couldn't open index [{}] --- ", indexPath);
            logger->error("Please make sure that 'salmon index' has been run successfully");
            std::exit(1);
        }
    }

    size_t numRecords = idx_->bns->n_seqs;
    std::vector<Transcript> transcripts_tmp;
    { // Load transcripts from file

            logger->info("Index contained {} targets; loading them", numRecords);
            //transcripts_.resize(numRecords);
            for (auto i : boost::irange(size_t(0), numRecords)) {
                uint32_t id = i;
                char* name = idx_->bns->anns[i].name;
                uint32_t len = idx_->bns->anns[i].len;
                // copy over the length, then we're done.
                transcripts_tmp.emplace_back(id, name, len);
            }

            std::sort(transcripts_tmp.begin(), transcripts_tmp.end(),
                    [](const Transcript& t1, const Transcript& t2) -> bool {
                    return t1.id < t2.id;
                    });
            double alpha = 0.005;
            char nucTab[256];
            nucTab[0] = 'A'; nucTab[1] = 'C'; nucTab[2] = 'G'; nucTab[3] = 'T';
            for (size_t i = 4; i < 256; ++i) { nucTab[i] = 'N'; }

            // Load the transcript sequence from file
            for (auto& t : transcripts_tmp) {
                transcripts_.emplace_back(t.id, t.RefName.c_str(), t.RefLength, alpha);
                // from BWA
                uint8_t* rseq = nullptr;
                int64_t tstart, tend, compLen, l_pac = idx_->bns->l_pac;
                tstart  = idx_->bns->anns[t.id].offset;
                tend = tstart + t.RefLength;
                rseq = bns_get_seq(l_pac, idx_->pac, tstart, tend, &compLen);
                if (compLen != t.RefLength) {
                    fmt::print(stderr,
                               "For transcript {}, stored length ({}) != computed length ({}) --- index may be corrupt. exiting\n",
                               t.RefName, compLen, t.RefLength);
                    std::exit(1);
                }

                std::string seq(t.RefLength, ' ');
                if (rseq != 0) {
                    for (size_t i = 0; i < compLen; ++i) {
                        seq[i] = rseq[i];
                                 //nst_nt4_table[static_cast<int>(nucTab[rseq[i]])];
                    }
                }
                transcripts_.back().Sequence = salmon::stringtools::encodeSequenceInSAM(seq.c_str(), t.RefLength);
                free(rseq);
                // end BWA code
            }
            // Since we have the de-coded reference sequences, we no longer need
            // the encoded sequences, so free them.
            free(idx_->pac); idx_->pac = nullptr;
            transcripts_tmp.clear();
            // ====== Done loading the transcripts from file
}
*/

int salmonIndex(int argc, char* argv[]) {

    using std::string;
    namespace bfs = boost::filesystem;
    namespace po = boost::program_options;

    bool useStreamingParser = true;

    uint32_t saSampInterval = 1;
    uint32_t maxThreads = std::thread::hardware_concurrency();
    uint32_t numThreads;

    po::options_description generic("Command Line Options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("transcripts,t", po::value<string>()->required(), "Transcript fasta file.")
    ("index,i", po::value<string>()->required(), "Salmon index.")
    ("threads,p", po::value<uint32_t>(&numThreads)->default_value(maxThreads)->required(),
                            "Number of threads to use (only used for computing bias features)")
    ("sasamp,s", po::value<uint32_t>(&saSampInterval)->default_value(1)->required(),
                            "The interval at which the suffix array should be sampled. "
                            "Smaller values are faster, but produce a larger index. "
                            "The default should be OK, unless your transcriptome is huge.")
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

        fmt::MemoryWriter optWriter;
        optWriter << vm["sasamp"].as<uint32_t>();

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

        bfs::path outputPrefix = indexDirectory / "bwaidx";

        std::vector<char const*> bwaArgVec{ "index",
                                    "-s",
                                    optWriter.str().c_str(),
                                    "-p",
                                    outputPrefix.string().c_str(),
                                    transcriptFile.c_str() };

        char* bwaArgv[] = { const_cast<char*>(bwaArgVec[0]),
                            const_cast<char*>(bwaArgVec[1]),
                            const_cast<char*>(bwaArgVec[2]),
                            const_cast<char*>(bwaArgVec[3]),
                            const_cast<char*>(bwaArgVec[4]),
                            const_cast<char*>(bwaArgVec[5]) };
        int bwaArgc = 6;

        ret = bwa_index(bwaArgc, bwaArgv);

        jointLog->info("done building BWT Index");

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

