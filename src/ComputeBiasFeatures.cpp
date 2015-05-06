
/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Salmon.

    Salmon is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Salmon is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Salmon.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/


#include <boost/thread/thread.hpp>

#include <iostream>
#include <vector>
#include <array>
#include <atomic>
#include <thread>

#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"
#include "jellyfish/mer_dna.hpp"

#include "tbb/concurrent_queue.h"

#include <boost/range/irange.hpp>
#include <boost/filesystem.hpp>

#include "CommonTypes.hpp"

// holding 2-mers as a uint64_t is a waste of space,
// but using Jellyfish makes life so much easier, so
// we'll live with it for now.
using Kmer = uint64_t;
using Sailfish::TranscriptFeatures;
namespace bfs = boost::filesystem;

template <typename ParserT>
bool computeBiasFeaturesHelper(ParserT& parser,
                               tbb::concurrent_bounded_queue<TranscriptFeatures>& featQueue,
                               size_t& numComplete, size_t numThreads) {

    using stream_manager = jellyfish::stream_manager<char**>;
    using sequence_parser = jellyfish::whole_sequence_parser<stream_manager>;

    size_t merLen = 2;
    Kmer lshift(2 * (merLen - 1));
    Kmer masq((1UL << (2 * merLen)) - 1);
    std::atomic<size_t> readNum{0};

    size_t numActors = numThreads;
    std::vector<std::thread> threads;
    auto tstart = std::chrono::steady_clock::now();

    for (auto i : boost::irange(size_t{0}, numActors)) {
        threads.push_back(std::thread(
	        [&featQueue, &numComplete, &parser, &readNum, &tstart, lshift, masq, merLen, numActors]() -> void {

                size_t cmlen, numKmers;
                jellyfish::mer_dna_ns::mer_base_dynamic<uint64_t> kmer(merLen);

                // while there are transcripts left to process
                while (true) { //producer.nextRead(s)) {
                    sequence_parser::job j(parser);
                    // If this job is empty, then we're done
                    if (j.is_empty()) { return; }

                    for (size_t i=0; i < j->nb_filled; ++i) {
                        ++readNum;
                        if (readNum % 100 == 0) {
                            auto tend = std::chrono::steady_clock::now();
                            auto sec = std::chrono::duration_cast<std::chrono::seconds>(tend-tstart);
                            auto nsec = sec.count();
                            auto rate = (nsec > 0) ? readNum / sec.count() : 0;
                            std::cerr << "processed " << readNum << " transcripts (" << rate << ") transcripts/s\r\r";
                        }

                        // we iterate over the entire read
                        const char* start     = j->data[i].seq.c_str();
                        uint32_t readLen      = j->data[i].seq.size();
                        const char* const end = start + readLen;

                        TranscriptFeatures tfeat{};

                        // reset all of the counts
                        numKmers = 0;
                        cmlen = 0;
                        kmer.polyA();

                        // the maximum number of kmers we'd have to store
                        uint32_t maxNumKmers = (readLen >= merLen) ? readLen - merLen + 1 : 0;
                        if (maxNumKmers == 0) { featQueue.push(tfeat); continue; }

                        // The transcript name
                        std::string fullHeader(j->data[i].header);
                        tfeat.name = fullHeader.substr(0, fullHeader.find(' '));
                        tfeat.length = readLen;
                        auto nfact = 1.0 / readLen;

                        // iterate over the read base-by-base
                        size_t offset{0};
                        size_t numChars{j->data[i].seq.size()};
                        while (offset < numChars) {
                            auto c = jellyfish::mer_dna::code(j->data[i].seq[offset]);
                            kmer.shift_left(c);
                            if (jellyfish::mer_dna::not_dna(c)) {
                                cmlen = 0;
                                ++offset;
                                continue;
                            }
                            if (++cmlen >= merLen) {
                                size_t twomer = kmer.get_bits(0, 2*merLen);
                                tfeat.diNucleotides[twomer]++;
                                switch(c) {
                                    case jellyfish::mer_dna::CODE_G:
                                    case jellyfish::mer_dna::CODE_C:
                                        tfeat.gcContent += nfact;
                                        break;
                                }
                            }
                            ++offset;
                        } // end while

                        char lastBase = j->data[i].seq.back();
                        auto c = jellyfish::mer_dna::code(lastBase);
                        switch(c) {
                            case jellyfish::mer_dna::CODE_G:
                            case jellyfish::mer_dna::CODE_C:
                                tfeat.gcContent += nfact;
                                break;
                        }

                        featQueue.push(tfeat);
                    } // end job
                } // end while(true)
            } // end lambda
            ));

        } // actor loop

        for (auto& t : threads) { t.join(); ++numComplete; }
        return true;
}

int computeBiasFeatures(
    std::vector<std::string>& transcriptFiles,
    bfs::path outFilePath,
    bool useStreamingParser,
    size_t numThreads) {

    using std::string;
    using std::vector;
    using std::cerr;

    size_t numActors = numThreads;
    size_t numComplete = 0;
    tbb::concurrent_bounded_queue<TranscriptFeatures> featQueue;

    std::ofstream ofile(outFilePath.string());

    auto outputThread = std::thread(
         [&ofile, &numComplete, &featQueue, numActors]() -> void {
             TranscriptFeatures tf{};
             while( numComplete < numActors or !featQueue.empty() ) {
                 while(featQueue.try_pop(tf)) {
                     ofile << tf.name << '\t';
                     ofile << tf.length << '\t';
                     ofile << tf.gcContent << '\t';
                     for (auto i : boost::irange(size_t{0}, tf.diNucleotides.size())) {
                         ofile << tf.diNucleotides[i];
                         char end = (i == tf.diNucleotides.size() - 1) ? '\n' : '\t';
                         ofile << end;
                     }
                 }
                 boost::this_thread::sleep_for(boost::chrono::milliseconds(100));

             }
             ofile.close();
         });

    for( auto rf : transcriptFiles) {
        std::cerr << "readFile: " << rf << ", ";
    }
    std::cerr << "\n";

    for (auto& readFile : transcriptFiles) {
        std::cerr << "file " << readFile << ": \n";

        //namespace bfs = boost::filesystem;
        //bfs::path filePath(readFile);

        char* pc = new char[readFile.size() + 1];
        std::strcpy(pc, readFile.c_str());
        char* fnames[] = {pc};

        // Create a jellyfish parser
        const int concurrentFile{1};

        using stream_manager = jellyfish::stream_manager<char**>;
        using sequence_parser = jellyfish::whole_sequence_parser<stream_manager>;
        stream_manager streams(fnames, fnames + 1, concurrentFile);

        size_t maxReadGroupSize{100};
        sequence_parser parser(1*numActors, maxReadGroupSize, concurrentFile, streams);
        computeBiasFeaturesHelper<sequence_parser>(
                parser, featQueue, numComplete, numActors);
        delete pc;
    }

    std::cerr << "\n";
    outputThread.join();
    return 0;
}
