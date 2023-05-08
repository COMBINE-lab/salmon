
#include "FastxParser.hpp"
#include <cmath>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>
#include <bitset>

#include "ProgOpts.hpp"
#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "PufferFS.hpp"
#include "ScopedTimer.hpp"
#include "Util.hpp"
#include "jellyfish/mer_dna.hpp"
#include "spdlog/spdlog.h"

#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "Util.hpp"
#include "PufferfishBinaryGFAReader.hpp"


uint8_t reverseBits(uint8_t b) {
  b = (b & 0xF0) >> 4 | (b & 0x0F) << 4;
  b = (b & 0xCC) >> 2 | (b & 0x33) << 2;
  b = (b & 0xAA) >> 1 | (b & 0x55) << 1;
  return b;
}

/*
enum class Direction : bool { FORWARD = 0, BACKWORD = 1 };

struct extension{
  char c;
  Direction dir ;
};

std::vector<extension> getEdges(uint8_t edgeVec){
  std::vector<extension> ext ;
  uint8_t mask = 1 ;
  std::vector<char> nuclmap = {'C','G','T','A','C','G','T','A'} ;
  for(uint8_t i=0;i<8;i++){
    if(edgeVec & (mask << i)){
      if(i<4)
        ext.push_back({nuclmap[i], Direction::FORWARD}) ;
      else
        ext.push_back({nuclmap[i], Direction::BACKWORD}) ;
    }
  }
  return ext ;

}
*/

template <typename IndexT>
bool checkKmer(IndexT& pi, CanonicalKmer& kb, pufferfish::util::QueryCache& qc, uint64_t kbi, uint32_t rn, uint32_t posWithinRef) {
  auto chits = pi.getRefPos(kb, qc);
  if (chits.empty()) {
    std::cerr << "k-mer " << kb.to_str() << " completely absent from the index\n\n";
    return false;
  }

  if (chits.refRange.size() < 0) {

    std::cerr << "There was a row in the contig table with an invalid range (size = " << chits.refRange.size() << "). "
                 "this means there is something wrong with the stored index.  Please report this issue on GitHub.\n";
    std::exit(1);
  }

  bool verbose = (kb.to_str() == "TGATATGATTTTTCTTCCTTTTTAAAATCTA");

  bool foundKmer{false};
  for (auto& rpos : chits.refRange) {
    auto refInfo = chits.decodeHit(rpos);
    if (verbose) {
    std::cout << "kmer : " << kb.to_str() << "\n";
    std::cout << "hit at ref : " << rpos.transcript_id() << " (name : " 
      << pi.getFullRefNames()[rpos.transcript_id()]
      << "), pos : " << refInfo.pos << "\n";
    }
    if (rpos.transcript_id() == rn) {
      if (refInfo.pos == posWithinRef) {
        if (refInfo.isFW) {
          foundKmer = (kb.fwWord() == kbi);
        } else {
          foundKmer = (kb.rcWord() == kbi);
        }
        break;
      }
    }
  }

  if (!foundKmer) {
    std::cerr << "couldn't find " << kb.to_str() << " where it actually occurs (" << posWithinRef << "), but found it at :\n";
    for (auto& rpos : chits.refRange) {
      auto refInfo = chits.decodeHit(rpos);
      std::cerr << "\ttr : " << rpos.transcript_id() << ", pos : " << refInfo.pos << "\n";
    }
    std::cerr << "\n";
  }
  return foundKmer;
}

/**
 * Check for internal consistency of the index
 **/

template <typename IndexT>
int doPufferfishInternalValidate(IndexT& pi, pufferfish::ValidateOptions& validateOpts) {
  (void)(validateOpts);
  int32_t k = pi.k();
  CanonicalKmer::k(k);

  auto console = spdlog::stderr_color_mt("console");

  auto& refSeq = pi.refseq_;

  // iterate over all reference sequences
//  const auto& refLengths = pi.getRefLengths();
  uint32_t rn{0};
  uint64_t gpos{0};
  uint64_t totalKmersSearched{0};
  uint64_t validCnt = pi.getIndexedRefCount();
  for (uint64_t i = 0; i < validCnt; i++) {
  
    auto refLen = pi.refLength(i);
    uint32_t posWithinRef{0};

    if (static_cast<int32_t>(refLen) <= k) {
      console->warn("reference sequence of length l = {} (< k = {}); not validating.", refLen, k);
      gpos += refLen;
      continue;
    }


    CanonicalKmer kb;
    uint64_t kbi = refSeq.get_int(2*gpos, 2*k);
    kb.fromNum(kbi);
    pufferfish::util::QueryCache qc;


    auto foundKmer = checkKmer(pi, kb, qc, kbi, rn, posWithinRef);
    if (!foundKmer) {
      console->error("Could not find k-mer ({}) occurring at position {} of reference {}, which is global position {}.", kb.to_str(), posWithinRef, rn, gpos);
      //return -1;
    }
    ++totalKmersSearched;

    gpos += k;
    for (uint32_t p = k; p < refLen; ++p) {
      ++posWithinRef;
      kb.shiftFw(static_cast<int>(refSeq[gpos]));
      kbi = kb.fwWord();

      auto foundKmer = checkKmer(pi, kb, qc, kbi, rn, posWithinRef);
      if (!foundKmer) {

        console->error("Could not find k-mer ({}) occurring at position {} of reference {}, which is global position {}.", kb.to_str(), posWithinRef, rn, gpos);
        //return -1;
      }

      ++totalKmersSearched;
      ++gpos;
      if (gpos % 10000000 == 0) {
        console->info("processed {} of {} positions.", gpos, refSeq.size() - (k*validCnt) + (validCnt));
      }
    }

    ++rn;
  }

  console->info("successfully looked up {} k-mer occurrences across {} reference sequences.", totalKmersSearched, rn);
  return 0;
}

template <typename IndexT>
int doPufferfishValidate(IndexT& pi, pufferfish::ValidateOptions& validateOpts) {
  //size_t k = pi.k();
  CanonicalKmer::k(pi.k());
  size_t found = 0;
  size_t notFound = 0;
  size_t correctPosCntr = 0;
  size_t incorrectPosCntr = 0;
  size_t numTrueTxp = 0;
  size_t numLostTxp = 0;

  if(validateOpts.gfaFileName.length() != 0){
    int k = pi.k() ;
    ScopedTimer st ;
    auto console = spdlog::stderr_color_mt("console");
    // NOTE: Should the false argument below be a command line option?
    // that is, should we consider building the edge vector here?
    pufferfish::BinaryGFAReader pf(validateOpts.gfaFileName.c_str(), k-1, false, false, console);
    pf.parseFile() ;

    auto& seq = pi.getSeq() ;
    //auto& edge = pi.getEdge() ;

    auto& contigid2seq = pf.getContigNameMap() ;

    //auto& paths = pf.getPaths() ;
    for(auto& ctg : contigid2seq){
      auto& ctgInfo = ctg.second ;

      uint64_t kbi, kei ;
      kbi = seq.get_int(2*ctgInfo.offset, 2*k) ;
      kei = seq.get_int(2*(ctgInfo.offset+ctgInfo.length-k), 2*k) ;
      CanonicalKmer kb,ke ;
      kb.fromNum(kbi) ;
      //check if rank is same as file order
      auto chits = pi.getRefPos(kb) ;
      if(chits.contigIdx_ != ctgInfo.fileOrder){
        std::cerr << "should not happen : " << "rank: "<< chits.contigIdx_ << " fileOrder: "<< ctgInfo.fileOrder << " offset " << ctgInfo.offset  << "\n" ;
        //std::exit(1) ;
      }

      ke.fromNum(kei) ;

      //std::cerr << "Start kmer " << kb.to_str() << "\n" ;
      //std::cerr << "End kmer " << ke.to_str() << "\n" ;

      //uint8_t edgeVec = edge[ctgInfo.fileOrder] ;
      //uint8_t loop = 0 ;
/*
      std::bitset<8> b = edgeVec ;
      std::cerr << "Parsing vector " << b << "\n" ;

      std::vector<pufferfish::util::extension> ext = pufferfish::util::getExts(edgeVec) ;
      for(auto& e : ext){
       
        auto kbtmp = kb ;
        auto ketmp = ke ;
        char c = e.c ;

        if(e.dir == pufferfish::util::Direction::FORWARD){
          ke.shiftFw(c) ;
          CanonicalKmer kt;
          kt.fromNum(ke.getCanonicalWord()) ;
          auto phits = pi.getRefPos(kt);
          if(phits.empty()){
            std::cerr<<" Should not happen " << ke.to_str() << " not found  " << (int)loop << " bit\n" ;
            std::exit(1) ;
          }
        }else{
          kb.shiftBw(c) ;
          CanonicalKmer kt;
          kt.fromNum(ke.getCanonicalWord()) ;
          auto phits = pi.getRefPos(kt);
          if(phits.empty()){
            std::cerr<<" Should not happen " << ke.to_str() << " not found  " << (int)loop << " bit\n" ;
            std::exit(1) ;
          }
        }

        kb = kbtmp ;
        ke = ketmp ;
      }
      */
      
      /*
      while(edgeVec){
        auto kbtmp = kb ;
        auto ketmp = ke ;
        uint8_t e = edgeVec & 0x01 ;
        edgeVec = edgeVec >> 1 ;
        switch(loop){
        case 0: if(e) {
            ke.shiftFw('C');
            auto phits = pi.getRefPos(ke);
            if(phits.empty()){
              std::cerr<<" Should not happen " << ke.to_str() << " not found " << (int)loop << " bit\n" ;
              std::exit(1) ;
            }
          }
          break ;
        case 1: if(e) {
            ke.shiftFw('G');
            auto phits = pi.getRefPos(ke);
            if(phits.empty()){
              std::cerr<<" Should not happen " << ke.to_str() << " not found " << (int)loop << " bit\n" ;
              std::exit(1) ;
            }
          }
          break ;
        case 2: if(e) {
            ke.shiftFw('T');
            CanonicalKmer kt;
            kt.fromNum(ke.getCanonicalWord()) ;
            auto phits = pi.getRefPos(kt);
            if(phits.empty()){
              std::cerr<<" Should not happen " << ke.to_str() << " not found  " << (int)loop << " bit\n" ;
              std::exit(1) ;
            }
          }
          break ;
        case 3: if(e) {
            auto kbackup = ke ;
            ke.shiftFw('A');
            CanonicalKmer kt;
            kt.fromNum(ke.getCanonicalWord()) ;
            auto phits = pi.getRefPos(kt);
            if(phits.empty()){
              std::cerr<<" Should not happen " << ke.to_str() << " not found  " << (int)loop << " bit\n" ;
              std::exit(1) ;
            }
          }
          break ;
        case 4: if(e) {
            kb.shiftBw('C');
            auto phits = pi.getRefPos(kb);
            if(phits.empty()){
              std::cerr<<" Should not happen " << kb.to_str() << "not found " << (int)loop << " bit\n" ;
              std::exit(1) ;
            }
          }
          break ;
        case 5: if(e) {
            kb.shiftBw('G');
            auto phits = pi.getRefPos(kb);
            if(phits.empty()){
              std::cerr<<" Should not happen " << kb.to_str() << "not found  " << (int)loop << " bit\n" ;
              std::exit(1) ;
            }
          }
          break ;
        case 6: if(e) {
            kb.shiftBw('T');
            auto phits = pi.getRefPos(kb);
            if(phits.empty()){
              std::cerr<<" Should not happen " << kb.to_str() << "not found " << (int)loop << " bit\n" ;
              std::exit(1) ;
            }
          }
          break ;
        case 7: if(e) {
            kb.shiftBw('A');
            auto phits = pi.getRefPos(kb);
            if(phits.empty()){
              std::cerr<<" Should not happen " << kb.to_str() << "not found " << (int)loop << " bit\n" ;
              std::exit(1) ;
            }
          }
          break ;

        }
        kb = kbtmp ;
        ke = ketmp ;
        loop++ ;
        }*/


      }
    
    std::cerr << "Found all edges \n" ;


  }

  
  {
    ScopedTimer st;
    std::vector<std::string> read_file = {validateOpts.refFile};
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(read_file, 1, 1);
    parser.start();
    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    size_t rn{0};
    // size_t kmer_pos{0};
    pufferfish::CanonicalKmerIterator kit_end;
    auto rg = parser.getReadGroup();
    while (parser.refill(rg)) {
      // Here, rg will contain a chunk of read pairs
      // we can process.
      for (auto& rp : rg) {
        // kmer_pos = 0;
        if (rn % 500000 == 0) {
          std::cerr << "rn : " << rn << "\n";
          std::cerr << "found = " << found << ", notFound = " << notFound
                    << "\n";
        }
        ++rn;
        auto& r1 = rp.seq;
        pufferfish::CanonicalKmerIterator kit1(r1);
        pufferfish::util::QueryCache qc;
        for (; kit1 != kit_end; ++kit1) {
          auto phits = pi.getRefPos(kit1->first, qc);
          if (phits.empty()) {
            ++notFound;
          } else {
            ++found;
            bool correctPos = false;
            bool foundTxp = false;
            bool cor = false;
            uint32_t clen = 0;
            std::vector<uint32_t> wrongPos;
            for (auto& rpos : phits.refRange) {
              if (pi.refName(rpos.transcript_id()) == rp.name) {
                foundTxp = true;
                auto refInfo = phits.decodeHit(rpos);
                if (static_cast<int>(refInfo.pos) == kit1->second) {
                  correctPos = true;
                } else {
                  cor = phits.contigOrientation_;
                  clen = phits.contigLen_;
                  wrongPos.push_back(refInfo.pos);
                }
              } else {
              }
            }
            if (foundTxp) {
              ++numTrueTxp;
              if (correctPos) {
                ++correctPosCntr;
              } else {
                std::cerr
                  << "txp = [" << rp.name << "], "
                  << "kmer = [" << kit1->first.to_str() << "], "
                    << "correct pos = " << kit1->second << ", found "
                    << pufferfish::util::str(wrongPos) << ", contig orientation = " << cor
                    << ", contig len = " << clen << "\n";
                   ///<< ", contig rank = " << pi.contigID(kit1->first) << '\n';
                  //<< ", contig rank = " << pi.contigID(kit1->first) << '\n';
                ++incorrectPosCntr;
              }
            } else {
              ++numLostTxp;
            }
          }
        }
      }
    }
  }
  std::cerr << "found = " << found << ", not found = " << notFound << "\n";
  std::cerr << "correctPos = " << correctPosCntr
            << ", incorrectPos = " << incorrectPosCntr << "\n";
  std::cerr << "corrextTxp = " << numTrueTxp << ", lostTxp = " << numLostTxp
            << "\n";

  return 0;
}


int pufferfishValidate(pufferfish::ValidateOptions& validateOpts) {

  auto indexDir = validateOpts.indexDir;
  std::string indexType;
    {
      std::ifstream infoStream(indexDir + "/info.json");
      cereal::JSONInputArchive infoArchive(infoStream);
      infoArchive(cereal::make_nvp("sampling_type", indexType));
      std::cerr << "Index type = " << indexType << '\n';
      infoStream.close();
    }

    /*compact::vector<uint64_t, 2> seq;
    seq.deserialize(validateOpts.indexDir+"seq.bin", false);

    uint64_t bits=100;
    for (uint64_t i = 10; i < bits; i++) {
        std::cerr << seq[i];
    }
    std::cerr << "\n";
    uint64_t i = 10;
    while (i < bits) {
      uint64_t end = bits-i>32?32:bits-i;
        auto w = seq.get_int(i*2, end*2);
        for (uint64_t j = 0; j < end*2; j+=2) {
            std::cerr << ((w >> j) & 0x3);
        }
        i+=32;
    }
    std::cerr << "\n";
    std::exit(3);
*/
    if (indexType == "sparse") { 
      PufferfishSparseIndex pi(validateOpts.indexDir);
      return doPufferfishInternalValidate(pi, validateOpts);
      //return doPufferfishValidate(pi, validateOpts);
    } else if (indexType == "dense") {
      PufferfishIndex pi(validateOpts.indexDir);
      return doPufferfishInternalValidate(pi, validateOpts);
      //return doPufferfishValidate(pi, validateOpts);
    }
    return 0;
}


