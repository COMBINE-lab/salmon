#include "FastxParser.hpp"
#include <cmath>
#include <chrono>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>

#include "CLI/Timer.hpp"
#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "PufferFS.hpp"
#include "Util.hpp"
#include "spdlog/spdlog.h"

#include "ProgOpts.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "Util.hpp"
#include "SAMWriter.hpp"

namespace kmers = combinelib::kmers;

template <typename IndexT>
int doPufferfishKmerQuery(IndexT& pi, 
                          fastx_parser::FastxParser<fastx_parser::ReadSeq>& parser,
                          std::mutex& iomut) { 
  int k = pi.k();
  (void)k;

  uint64_t processed = 0;
  uint64_t buff_size = 10000;
  std::string workstr;

  std::ostringstream osstream;
    

  size_t found = 0;
  size_t notFound = 0;
  size_t totalHits = 0;
  {
    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    size_t rn{0};
    pufferfish::util::QueryCache qc;

    CanonicalKmer km;
    auto rg = parser.getReadGroup();
    while (parser.refill(rg)) {
      // Here, rg will contain a chunk of read pairs
      // we can process.
      for (auto& rp : rg) {
        // kmer_pos = 0;
        /*
        if (rn % 500000 == 0) {
          std::cerr << "kn : " << rn << "\n";
          std::cerr << "found = " << found << ", notFound = " << notFound
                    << ", total hits = " << totalHits << "\n";
        }
        */
        ++rn;
        auto& r1 = rp.seq;

        km.fromStr(r1.data());
        auto phits = pi.getRefPos(km, qc);
        if (phits.empty()) {
          ++notFound;
          osstream << rp.name << "\t0\n";
        } else {
          ++found;
          size_t nh = phits.refRange.size();
          totalHits += nh;
          osstream << rp.name << '\t' << nh << '\n';
          for (auto& rh : phits.refRange) {
            auto ref_id = pi.getRefId(rh.transcript_id());
            auto h = phits.decodeHit(rh);
            osstream << ref_id << '\t' << h.pos << '\t' << (h.isFW ? 'f' : 'r') << "\n";
          }
        }

        if (processed >= buff_size) {
            std::string o = osstream.str();
            iomut.lock();
            std::cout << o;
            iomut.unlock();
            osstream.clear();
            osstream.str("");
            processed = 0;
        }
        processed += 1;
      }

    }
  }
  
  // dump any remaining output
  std::string o = osstream.str();
  iomut.lock();
  std::cout << o;
  iomut.unlock();
  osstream.clear(); 
  osstream.str("");

  iomut.lock();
  std::cerr << "found = " << found << ", not found = " << notFound << "\n";
  std::cerr << "total hits = " << totalHits << "\n";
  iomut.unlock();
  return 0;
}


int pufferfishKmerQuery(pufferfish::KmerQueryOptions& kqueryOpts) {
  std::ios_base::sync_with_stdio(false);
  auto indexDir = kqueryOpts.indexDir;
  std::string indexType;
  {
    std::ifstream infoStream(indexDir + "/info.json");
    cereal::JSONInputArchive infoArchive(infoStream);
    infoArchive(cereal::make_nvp("sampling_type", indexType));
    std::cerr << "Index type = " << indexType << '\n';
    infoStream.close();
  }

  std::mutex iomut;
  size_t nthread = kqueryOpts.num_threads;
  std::vector<std::string> &read_file = kqueryOpts.queryFiles;
  uint32_t np = 1;
  if ((read_file.size() > 1) and (nthread >= 6)) { np += 1; nthread -= 1;}
  fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(read_file, nthread, np);

  parser.start();

  std::vector<std::thread> workers;

  if (indexType == "sparse") { 
    PufferfishSparseIndex pi(kqueryOpts.indexDir);
    CanonicalKmer::k(pi.k());
    writeSAMHeader(pi, std::cout);

    for (size_t i = 0; i < nthread; ++i) {
        workers.push_back(std::thread([&pi, &parser, &iomut]() {
            doPufferfishKmerQuery(pi, parser, iomut);
        }));
    }

    for (auto& w : workers) {
      w.join();
    } 
  } else if (indexType == "dense") {
    PufferfishIndex pi(kqueryOpts.indexDir);
    CanonicalKmer::k(pi.k());
    writeSAMHeader(pi, std::cout);

    for (size_t i = 0; i < nthread; ++i) {
        workers.push_back(std::thread([&pi, &parser, &iomut]() {
            doPufferfishKmerQuery(pi, parser, iomut);
        }));
    }
    
    for (auto& w : workers) {
      w.join();
    } 
  }

  parser.stop();
  return 0;
}

