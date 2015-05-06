/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Sailfish.

    Sailfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sailfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sailfish.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/


#ifndef __MULTITHREADED_BAM_PARSER__
#define __MULTITHREADED_BAM_PARSER__

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <thread>
#include <atomic>
#include <iostream>

#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>
#include "tbb/concurrent_queue.h"

namespace bfs = boost::filesystem;

struct ReadSeq {
    char* seq = nullptr;
    size_t len = 0;
    char* name = nullptr;
    size_t nlen = 0;
};

class MultithreadedBAMParser {
public:
    MultithreadedBAMParser( std::vector<bfs::path>& files );
    ~MultithreadedBAMParser();
    bool start();
    bool nextAlignment(BamTools::BamAlignment*& seq);
    void finishedWithAlignment(BamTools::BamAlignment*& s);

private:
    std::vector<bfs::path>& inputStreams_;
    bool parsing_;
    std::thread* parsingThread_;
    tbb::concurrent_bounded_queue<BamTools::BamAlignment*> alnQueue_, seqContainerQueue_;
    BamTools::BamAlignment* alnStructs_;
    const size_t queueCapacity_ = 5000000;
};

//#include "Parser.cpp"

#endif // __MULTITHREADED_BAM_PARSER__
