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


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

int performBiasCorrection(
        bfs::path featureFile,
        bfs::path expressionFile,
        double estimatedReadLength,
        double kmersPerRead,
        uint64_t mappedKmers,
        uint32_t merLen,
        bfs::path outputFile,
        size_t numThreads);

int main(int argc, char* argv[]) {
        bfs::path featureFile(argv[1]);
        bfs::path expressionFile(argv[2]);
        double estimatedReadLength = atod(argv[3]);
        double kmersPerRead = atod(argv[4]);
        uint64_t mappedKmers = atol(argv[5]);
        uint32_t mappedKmers = atoi(argv[6]);
        bfs::path outputFile(argv[7]);
        size_t numThreads = atoi(argv[8]);

        performBiasCorrection(featureFile, expressionFile, estimatedReadLength, kmersPerRead,
                              mappedKmers, merLen, outputFile, numThreads);
}
