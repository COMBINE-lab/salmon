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


#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <sstream>
#include <memory>
#include <functional>
#include <unordered_map>
#include <mutex>
#include <thread>
#include <chrono>
#include <iomanip>

#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"

#include <boost/range/irange.hpp>
#include "ezETAProgressBar.hpp"
#include "LookUpTableUtils.hpp"


namespace LUTTools {

/**
 *  \brief Dump the k-mer memberships vector to the file fname
 **/
void dumpKmerEquivClasses(
                          const std::vector<KmerID>& memberships,
                          const std::string& fname) {

  std::ofstream ofile(fname, std::ios::binary);

  // Write the length of the vector to the file
  size_t vecLen{memberships.size()};
  ofile.write(reinterpret_cast<const char*>(&vecLen), sizeof(vecLen));

  // Write the actual vector to the file
  ofile.write(reinterpret_cast<const char*>(&memberships.front()), sizeof(memberships.front()) * vecLen);

  // Close the file
  ofile.close();
}


std::vector<KmerID> readKmerEquivClasses(const std::string& fname) {
  std::ifstream ifile(fname, std::ios::binary);
  size_t vecLen{0};
  ifile.read(reinterpret_cast<char*>(&vecLen), sizeof(vecLen));

  std::vector<KmerID> memberships(vecLen, std::numeric_limits<KmerID>::max());
  ifile.read(reinterpret_cast<char*>(&memberships.front()), sizeof(memberships.front()) * vecLen);

  ifile.close();
  return memberships;
}

void dumpKmerLUT(
    std::vector<TranscriptList> &transcriptsForKmerClass,
    const std::string &fname) {

    tbb::parallel_for_each( transcriptsForKmerClass.begin(), transcriptsForKmerClass.end(),
    [&]( TranscriptList & t ) {
        std::sort(t.begin(), t.end());
    });

    std::ofstream ofile(fname, std::ios::binary);
    // number of kmers
    auto numk = transcriptsForKmerClass.size();
    ofile.write(reinterpret_cast<const char *>(&numk), sizeof(numk));
    for (auto i : boost::irange(size_t(0), numk)) {
        // write the vector's size
        auto s = transcriptsForKmerClass[i].size();
        ofile.write(reinterpret_cast<const char *>(&s), sizeof(s));
        // write the vector's contents
        if ( s > 0 ) {
            ofile.write(reinterpret_cast<const char *>(&transcriptsForKmerClass[i][0]), s * sizeof(transcriptsForKmerClass[i][0]));
        }
    }
    // close the output file
    ofile.close();
}

void readKmerLUT(
    const std::string &fname,
    std::vector<TranscriptList> &transcriptsForKmer) {

    std::ifstream ifile(fname, std::ios::binary);
    // get the size of the vector from file
    size_t numk = 0;
    ifile.read(reinterpret_cast<char *>(&numk), sizeof(numk));
    // resize the vector now
    transcriptsForKmer.resize(numk);

    for (auto i : boost::irange(size_t(0), numk)) {
        // read the vector's size
        size_t numTran = 0;
        ifile.read(reinterpret_cast<char *>(&numTran), sizeof(numTran));
        // read the vector's contents
        if ( numTran > 0 ) {
            transcriptsForKmer[i].resize(numTran);
            ifile.read(reinterpret_cast<char *>(&transcriptsForKmer[i][0]), numTran * sizeof(TranscriptID));
        }
    }

    ifile.close();
}


void writeTranscriptInfo (TranscriptInfo *ti, std::ofstream &ostream) {
    size_t numKmers = ti->kmers.size();
    size_t recordSize = sizeof(ti->transcriptID) +
                        sizeof(ti->geneID) +
                        sizeof(ti->name.length()) +
                        ti->name.length() +
                        sizeof(ti->length) +
                        sizeof(numKmers) +
                        sizeof(KmerID) * numKmers;

    ostream.write(reinterpret_cast<const char *>(&recordSize), sizeof(recordSize));
    ostream.write(reinterpret_cast<const char *>(&ti->transcriptID), sizeof(ti->transcriptID));
    ostream.write(reinterpret_cast<const char *>(&ti->geneID), sizeof(ti->geneID));
    auto l = ti->name.length();
    ostream.write(reinterpret_cast<const char *>(&l), sizeof(l));
    ostream.write(reinterpret_cast<const char *>(ti->name.c_str()), l);
    ostream.write(reinterpret_cast<const char *>(&ti->length), sizeof(ti->length));
    ostream.write(reinterpret_cast<const char *>(&numKmers), sizeof(numKmers));
    ostream.write(reinterpret_cast<const char *>(&ti->kmers[0]), numKmers * sizeof(KmerID));
}

std::unique_ptr<TranscriptInfo> readTranscriptInfo(std::ifstream &istream) {
    std::unique_ptr<TranscriptInfo> ti(new TranscriptInfo);
    size_t recordSize = 0;
    istream.read(reinterpret_cast<char *>(&recordSize), sizeof(recordSize));
    istream.read(reinterpret_cast<char *>(&ti->transcriptID), sizeof(ti->transcriptID));
    istream.read(reinterpret_cast<char *>(&ti->geneID), sizeof(ti->geneID));
    size_t slen = 0;
    istream.read(reinterpret_cast<char *>(&slen), sizeof(slen));
    char *name = new char[slen+1];
    name[slen] = '\0';
    istream.read(name, slen);
    ti->name = std::string(name);
    // read the transcript's length
    istream.read(reinterpret_cast<char *>(&ti->length), sizeof(ti->length));
    size_t numKmers = 0;
    istream.read(reinterpret_cast<char *>(&numKmers), sizeof(numKmers));
    ti->kmers = std::vector<KmerID>(numKmers, 0);
    istream.read(reinterpret_cast<char *>(&ti->kmers[0]), sizeof(KmerID)*numKmers);
    return ti;
}

/*
std::vector<std::unique_ptr<TranscriptInfo>> getTranscriptsFromFile(const std::string &tlutfname,
        const std::vector<Offset> &offsets,
std::pair<CTVecIt, CTVecIt> be) {

    Offset INVALID = std::numeric_limits<Offset>::max();
    std::cerr << "before creating transcript vector" << std::endl;
    auto batchSize = std::distance(be.first, be.second);
    std::cerr << "batch size = " << batchSize << "\n";
    std::vector<std::unique_ptr<TranscriptInfo>> transcripts;//(batchSize);
    //transcripts.reserve();
    std::cerr << "after creating transcript vector " << std::endl;

    std::cerr << "opening input file" << std::endl;
    std::ifstream ifile(tlutfname, std::ios::binary);
    std::cerr << "opened" << std::endl;
    size_t idx = 0;
    std::cerr << "iterating through batch" << std::endl;
    for (auto it = be.first; it != be.second; ++it) {
        std::cerr << "getting offset for transcript " << it->tid << std::endl;
        auto offset = offsets[it->tid];
        if ( offset == INVALID ) {
            std::cerr << "encountered invalid offset for transcript " << it->tid << std::endl;
            std::exit(1);
        }
        std::cerr << "seeking to " << offset << std::endl;
        ifile.seekg(offset);
        std::cerr << "reading in transcript from file" << std::endl;
        transcripts.emplace_back(readTranscriptInfo(ifile));
        ++idx;
        // Check that the transcript we read is the one we expected
        if (it->tid != transcripts.back()->transcriptID) {
            std::cerr << "read the wrong transcript!\n";
            std::exit(1);
        }
    }

    ifile.close();
    std::cerr << "done reading all transcripts from this batch" << std::endl;

    return std::move(transcripts);
}
*/

std::vector<Offset> buildTLUTIndex(const std::string &tlutfname, size_t numTranscripts) {

    std::cerr << "creating vector of " << numTranscripts << " offsets \n";
    Offset INVALID = std::numeric_limits<Offset>::max();
    std::vector<Offset> offsets(numTranscripts, INVALID);
    std::cerr << "done\n";

    std::cerr << "opening file\n";
    std::ifstream ifile(tlutfname, std::ios::binary);
    std::cerr << "done\n";

    size_t numRecords {0};
    std::cerr << "reading numRecords\n";
    ifile.read(reinterpret_cast<char *>(&numRecords), sizeof(numRecords));
    std::cerr << "numRecords = " << numRecords << std::endl;
    std::cerr << "numTranscripts = " << numTranscripts << std::endl;

    std::cerr << "building transcript lookup table index" << std::endl;
    ez::ezETAProgressBar pb(numRecords);
    pb.start();
    Offset offset = sizeof(numRecords);
    for (auto recNum : boost::irange(size_t(0), numRecords)) {
        size_t recordSize = 0;
        TranscriptID tid = 0;
        ifile.read(reinterpret_cast<char *>(&recordSize), sizeof(recordSize));
        ifile.read(reinterpret_cast<char *>(&tid), sizeof(tid));
        offsets[tid] = offset;
        offset += recordSize + sizeof(recordSize);
        ifile.seekg(offset);
        ++pb;
    }

    ifile.close();
    return offsets;
}
}
