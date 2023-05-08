
#include <cstdint>
#include <vector>
#include <iostream>
#include <fstream>

#include "sparsepp/spp.h"
#include "FastxParser.hpp"
#include "Kmer.hpp"
#include "clipp.h"
#include "string_view.hpp"

struct pconfig {
    uint32_t k=31;
    std::string ref;
    std::string gfa;
};

using KmerT = combinelib::kmers::Kmer<32,1>;
enum class OccT : uint8_t { START, END, BOTH };

struct ContigInfo {
	OccT type;
	uint32_t id;
	uint32_t length;
};

struct KmerContigMap {
	spp::sparse_hash_map<uint64_t, OccT> kmers;
    //spp::sparse_hash_set<uint64_t> skmer;
    //spp::sparse_hash_set<uint64_t> ekmer;
};

struct UnitigMap {
    spp::sparse_hash_map<uint64_t, ContigInfo> kmers;
    //spp::sparse_hash_map<uint64_t, std::pair<uint32_t, uint32_t>> skmer;
    //spp::sparse_hash_map<uint64_t, std::pair<uint32_t, uint32_t>> ekmer;
};


KmerContigMap
getTerminalKmers(fastx_parser::FastxParser<fastx_parser::ReadSeq>& parser, uint32_t k) {
    KmerContigMap km;
    KmerT::k(k);
    KmerT start;
    KmerT end;
  // Get the read group by which this thread will
  // communicate with the parser (*once per-thread*)
  auto rg = parser.getReadGroup();

  while (parser.refill(rg)) {
    // Here, rg will contain a chunk of read pairs
    // we can process.
    for (auto& rp : rg) {
      auto& r1 = rp.seq;
      start.fromChars(r1.begin());
      end.fromChars(r1.end() - k);
      if (start == end) {
	      auto it = km.kmers.find(start.word(0));
	      if (it == km.kmers.end()) {
		      km.kmers[start.word(0)] = OccT::BOTH;
	      }
      } else {
      	auto sit = km.kmers.find(start.word(0));
	if (sit != km.kmers.end()) {
		sit->second = (sit->second == OccT::END or sit->second == OccT::BOTH) ? OccT::BOTH : OccT::START;
	} else {
		km.kmers[start.word(0)] = OccT::START;
	}
	auto eit = km.kmers.find(end.word(0));
	if (eit != km.kmers.end()) {
		eit->second = (eit->second == OccT::START or eit->second == OccT::BOTH) ? OccT::BOTH : OccT::END;
	} else {
		km.kmers[end.word(0)] = OccT::END;
	}
      }
      //km.skmer.insert(start.word(0));
      //km.ekmer.insert(end.word(0));
    }
  }
  return km;
}

uint32_t saveUnitig(stx::string_view seq, uint32_t id, std::ofstream& ofile) {
    ofile << "S\t" << id << "\t" << seq << '\n'; 
    return id + 1;
}

uint32_t createUnitig(UnitigMap& umap, uint32_t k, const std::string& header, stx::string_view seqin, 
                      uint32_t id, std::ofstream& ofile) {
  (void) header;
    stx::string_view seq;
    std::string str;
    if (seq.length() == k) {
        KmerT ut(seq.begin());
        auto cut = ut.getCanonical();
        if (ut == cut) {
            seq = seqin;
        } else {
          str = ut.getCanonical().toStr();
          seq = str;
        }
    } else {
        seq = seqin;
    }
    auto nid = saveUnitig(seq, id, ofile);

    KmerT canonFirst(seq.begin());
    canonFirst.canonicalize();
    /*
    if (umap.skmer.contains(canonFirst.word(0)) or umap.ekmer.contains(canonFirst.word(0))) {
        std::cerr << "Error: Initial kmer " << canonFirst.toStr() << " is repeated.\n";
        std::exit(1);
    }
    */  
    KmerT canonLast(seq.end() - k);
    canonLast.canonicalize();
    /*
    if (umap.skmer.contains(canonLast.word(0)) or umap.ekmer.contains(canonLast.word(0))) {
        std::cerr << "Error: Last kmer is repeated.\n";
        std::exit(1);
    }
    */
    /*auto lastIt = umap.kmers.find(canonLast.word(0));
    if ( lastIt != umap.kmers.end() ) {
	    lastIt->kind = true;
	    lastIt->endID = id;
	    lastIt->endLen = seq.length();
    } else {
	    umap.kmers[
    }
    */

    if (canonLast == canonFirst) {
	    ContigInfo info {OccT::BOTH, id, static_cast<uint32_t>(seq.length())};
    	umap.kmers.emplace(canonLast.word(0), info);
    } else {
	    ContigInfo infoe {OccT::END, id, static_cast<uint32_t>(seq.length())};
    	umap.kmers.emplace(canonLast.word(0), infoe); 
	    ContigInfo infos {OccT::START, id, static_cast<uint32_t>(seq.length())};
      umap.kmers.emplace(canonFirst.word(0), infos);
    }

    /*
    umap.ekmer.emplace(canonLast.word(0), std::make_pair(id, seq.length()));
    umap.skmer.emplace(canonFirst.word(0), std::make_pair(id, seq.length()));
    */
    return nid;
}

UnitigMap splitUnitigs(fastx_parser::FastxParser<fastx_parser::ReadSeq>& parser, KmerContigMap& km, uint32_t k, std::ofstream& ufile) {
  // Get the read group by which this thread will
  // communicate with the parser (*once per-thread*)
  UnitigMap umap;
  auto rg = parser.getReadGroup();
  uint32_t unitigId{0};
  KmerT mer;
  uint32_t numProcessed{0};
  while (parser.refill(rg)) {
    for (auto& rp : rg) {
        if (numProcessed % 100000 == 0 and numProcessed > 0) {
            std::cout << "created " << numProcessed << " unitigs\n";
        }
      auto& h = rp.name;
      auto& seq = rp.seq;
      stx::string_view seqview = seq;
      uint32_t prev{0};
      bool first{true};
      mer.fromChars(seq.begin());
      for (uint32_t i = 0; i < seq.length() - k + 1; ++i) {
          if (!first) {
              mer.append(seq[i+k-1]);
          }
          first = false;
          auto rcmer = mer.getRC();
	  auto merIt = km.kmers.find(mer.word(0));
	  auto rcmerIt = (mer == rcmer) ? merIt : km.kmers.find(rcmer.word(0));
	  bool startFW = false;
	  bool endFW = false;
	  bool startRC = false;
	  bool endRC = false;
	  if (merIt != km.kmers.end()) {
		  startFW = (merIt->second == OccT::START or merIt->second == OccT::BOTH) ? true : false;
		  endFW = (merIt->second == OccT::END or merIt->second == OccT::BOTH) ? true : false;
	  }
	  if (rcmerIt != km.kmers.end()) {
		  startRC = (rcmerIt->second == OccT::START or rcmerIt->second == OccT::BOTH) ? true : false;
		  endRC = (rcmerIt->second == OccT::END or rcmerIt->second == OccT::BOTH) ? true : false;
	  }

          if (startFW or endRC) {//km.skmer.contains(mer.word(0)) or km.ekmer.contains(rcmer.word(0))) {
              if (i + k - 1 - prev >= k) {
              uint64_t len = i + k - 1  - prev;
              unitigId = createUnitig(umap, k, h, seqview.substr(prev, len), unitigId, ufile);
              prev = i;
              }
          }
          if (endFW or startRC) {//km.ekmer.contains(mer.word(0)) or km.skmer.contains(rcmer.word(0))) {
              uint64_t len = i + k - prev;
              unitigId = createUnitig(umap, k, h, seqview.substr(prev, len), unitigId, ufile);
              prev = i + 1;
          }
      } 
      if (seq.length() - prev >= k) {
          unitigId = createUnitig(umap, k, h, seqview.substr(prev), unitigId, ufile);
      }
      ++numProcessed;
    }
  }
  return umap;
}

bool buildPaths(fastx_parser::FastxParser<fastx_parser::ReadSeq>& parser, UnitigMap& umap, uint32_t k, std::ofstream& ufile) {
    uint32_t pctr{0};
    KmerT mer;
    char ori = '*';
    auto rg = parser.getReadGroup();
    //const auto itEnd = umap.kmers.end();
    //const auto skend = umap.skmer.end();
    //const auto ekend = umap.ekmer.end();
    while (parser.refill(rg)) {
        for (auto& rp : rg) { 
            auto& ref = rp.seq;
            uint32_t i{0};
            ufile << "P\t" << rp.name << '\t';

            //bool first{true};
            while (i < ref.length() - k + 1) {
                mer.fromChars(ref.begin() + i); 
                auto nkmer = mer.getCanonical();
		auto kIt = umap.kmers.find(nkmer.word(0));
                //auto skIt = umap.skmer.find(nkmer.word(0));
                //auto ekIt = umap.ekmer.find(nkmer.word(0));
                std::decay<decltype(umap.kmers[nkmer.word(0)])>::type unitigInfo;
                if (kIt->second.type == OccT::BOTH) {//umap.skmer.contains(nkmer.word(0)) and umap.ekmer.contains(nkmer.word(0))) {
                    unitigInfo = kIt->second;//umap.skmer[nkmer.word(0)];
                    ori = (mer == nkmer) ? '+' : '-';
                } else if (kIt->second.type == OccT::START) {///umap.skmer.contains(nkmer.word(0))) {
                    unitigInfo = kIt->second;//umap.skmer[nkmer.word(0)];
                    ori = '+';
                } else if (kIt->second.type == OccT::END) {//umap.ekmer.contains(nkmer.word(0))) {
                    unitigInfo = kIt->second;//umap.ekmer[nkmer.word(0)];
                    ori = '-';
                } else {
                    std::cerr << pctr << " paths reconstructed.\n";
                    std::cerr << "ERROR: kmer is not found in the start or end of a unitig\n";
                    std::cerr << mer.toStr() << "   ,   " << nkmer.toStr() << "\n";
                    std::exit(1);
                }
                i += unitigInfo.length - k + 1;
                char comma = (i < ref.length() - k + 1) ? ',' : '\n';
                ufile << unitigInfo.id << ori << comma;
            }
            ++pctr;
	    //std::cerr << "wrote " << pctr << " paths \n";
        }
    }
    ufile << '\n';
    return true;
}

int main(int argc, char* argv[]) {
    using namespace clipp;
    pconfig config;

    auto cli = (
        required("-r", "--ref") & value("reference file", config.ref) % "reference fasta file",
        required("-u", "--uinitigs") & value("BCALM unitig file", config.gfa) % "unitig file",
        option("-k", "--klen") & value("length of k", config.k) % "length of kmer (default : 31)"
    );

    if(!parse(argc, argv, cli)) std::cout << make_man_page(cli, argv[0]);
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser({config.ref}, 1, 1);
    parser.start();
    auto kmap = getTerminalKmers(parser, config.k);
    parser.stop();
    std::cout << "have " << kmap.kmers.size() << " contigs\n";
    
    std::cout << "splitting unitigs\n";

    std::string ufname = config.gfa + ".pufferized.gfa";
    std::ofstream ufile(ufname);
    fastx_parser::FastxParser<fastx_parser::ReadSeq> uparser({config.gfa}, 1, 1);
    uparser.start();
    auto umap = splitUnitigs(uparser, kmap, config.k, ufile);
    uparser.stop();
    std::cerr << "done\n";
    std::cerr << "umap size = " << umap.kmers.size() << "\n";
    
    std::cerr << "start reconstructing the paths\n";
    fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser({config.ref}, 1, 1);
    rparser.start();
    buildPaths(rparser, umap, config.k, ufile);
    rparser.stop();
    std::cerr << "done!\n";
    return 0;
}
