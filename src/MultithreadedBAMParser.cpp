/**
>HEADER
    Copyright (c) 2014 Rob Patro robp@cs.cmu.edu

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


#include "MultithreadedBAMParser.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <thread>
#include <atomic>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>
#include "tbb/concurrent_queue.h"


namespace bfs = boost::filesystem;
using BamTools::BamReader;
using BamTools::BamAlignment;

MultithreadedBAMParser::MultithreadedBAMParser(std::vector<bfs::path>& files) : inputStreams_(files),
        parsing_(false), parsingThread_(nullptr)
    {
        alnStructs_ = new BamTools::BamAlignment[queueCapacity_];
        alnQueue_.set_capacity(queueCapacity_);
        seqContainerQueue_.set_capacity(queueCapacity_);
        for (size_t i = 0; i < queueCapacity_; ++i) {
            seqContainerQueue_.push(&alnStructs_[i]);
        }
    }

MultithreadedBAMParser::~MultithreadedBAMParser() {
        parsingThread_->join();
        delete [] alnStructs_;
        delete parsingThread_;
    }

bool MultithreadedBAMParser::start() {
        if (!parsing_) {
            parsing_ = true;
            parsingThread_ = new std::thread([this](){

                BamTools::BamReader reader;

                for (auto file : this->inputStreams_) {
                    reader.Open(file.string());
                    std::cerr << "reading from " << file.native() << "\n";

                    BamAlignment* aln;
                    this->seqContainerQueue_.pop(aln);
                    while (reader.GetNextAlignment(*aln)) {
                        this->alnQueue_.push(aln);
                        this->seqContainerQueue_.pop(aln);
                    }
                    // close the file
                    reader.Close();
                }

                this->parsing_ = false;
            });
            return true;
        } else {
            return false;
        }

    }

bool MultithreadedBAMParser::nextAlignment(BamAlignment*& seq) {
        while(parsing_ or !alnQueue_.empty()) {
            if (alnQueue_.try_pop(seq)) { return true; }
        }
        return false;
}

void MultithreadedBAMParser::finishedWithAlignment(BamAlignment*& s) { seqContainerQueue_.push(s); }
