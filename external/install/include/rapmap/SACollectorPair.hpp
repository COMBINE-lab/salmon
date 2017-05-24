//
// RapMap - Rapid and accurate mapping of short reads to transcriptomes using
// quasi-mapping.
// Copyright (C) 2015, 2016 Rob Patro, Avi Srivastava, Hirak Sarkar
//
// This file is part of RapMap.
//
// RapMap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RapMap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RapMap.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef SA_COLLECTOR_PAIR_HPP
#define SA_COLLECTOR_PAIR_HPP

#include "RapMapSAIndex.hpp"
#include "RapMapUtils.hpp"
#include "SASearcher.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>

template <typename RapMapIndexT, typename MerT> class SACollectorPair {
public:
  using OffsetT = typename RapMapIndexT::IndexType;

  /** Disable NIP skipping **/
  void disableNIP() { disableNIP_ = true; }

  /** Enable NIP skipping --- default state **/
  void enableNIP() { disableNIP_ = false; }

  /** Require a coverage fraction of at least req for all reported mappings **/
  void setCoverageRequirement(double req) { covReq_ = req; }

  /** Get the current coverage requirement for mappings (0 means no requirement)
   * **/
  double getCoverageRequirement() const { return covReq_; }

  /** If any hit has a suffix array interval of this length or larger, just skip
   * it **/
  void setMaxInterval(OffsetT maxInterval) { maxInterval_ = maxInterval; }

  void setMMPThreshold(uint32_t mmpThreshold) { mmpThreshold_ = mmpThreshold; }

  /** Get the maximum allowable suffix array interval **/
  OffsetT getMaxInterval(OffsetT maxInterval) const { return maxInterval_; }

  /** Get/Set usage of strict-checking **/
  bool getStrictCheck() const { return strictCheck_; };
  void setStrictCheck(bool sc) { strictCheck_ = sc; }

  /** Construct an SACollector given an index **/
  SACollectorPair(RapMapIndexT* rmi)
      : rmi_(rmi), hashEnd_(rmi->khash.end()), disableNIP_(false),
        covReq_(0.0), maxInterval_(1000),
        strictCheck_(false) {}

  enum HitStatus { ABSENT = -1, UNTESTED = 0, PRESENT = 1 };
  // Record if k-mers are hits in the
  // fwd direction, rc direction or both
  //template <typename MerT>
  struct KmerDirScore {
    KmerDirScore(MerT kmerIn, int32_t kposIn,
                 HitStatus fwdScoreIn, HitStatus rcScoreIn)
        : kmer(kmerIn), kpos(kposIn), fwdScore(fwdScoreIn), rcScore(rcScoreIn) {
    }
    KmerDirScore() : kpos(0), fwdScore(UNTESTED), rcScore(UNTESTED) {}
    bool operator==(const KmerDirScore& other) const {
      return kpos == other.kpos;
    }
    bool operator<(const KmerDirScore& other) const {
      return kpos < other.kpos;
    }
    void print() {
      std::cerr << "{ " << kmer.to_str() << ", " << kpos << ", "
                << ((fwdScore) ? "PRESENT" : "ABSENT") << ", "
                << ((rcScore) ? "PRESENT" : "ABSENT") << "}\t";
    }
    //rapmap::utils::my_mer kmer;
    MerT kmer ;
    int32_t kpos;
    HitStatus fwdScore;
    HitStatus rcScore;
  };

  //template <typename MerT>
  bool operator()(std::string& read,
                  std::vector<rapmap::utils::SAIntervalHit<OffsetT>>& fwdSAInts,
                  std::vector<rapmap::utils::SAIntervalHit<OffsetT>>& rcSAInts,
                  SASearcher<RapMapIndexT>& saSearcher,
                  rapmap::utils::MateStatus mateStatus,
                  std::vector<uint32_t> goldenTids,
                  bool remap = false,
                  uint32_t mmpThreshold=15,
                  bool consistentHits = false) {

    using QuasiAlignment = rapmap::utils::QuasiAlignment;
    using MateStatus = rapmap::utils::MateStatus;
    using SAIntervalHit = rapmap::utils::SAIntervalHit<OffsetT>;

    auto& txpNames = rmi_->txpNames ;
    auto& rankDict = rmi_->rankDict;
    auto& txpStarts = rmi_->txpOffsets;
    auto& SA = rmi_->SA;
    //auto& khash = (remap)?rmi_->khash9:rmi_->khash;
    auto& khash = rmi_->khash;
    auto& text = rmi_->seq;
    auto salen = SA.size();
    if(remap) hashEnd_ = khash.end();
    auto readLen = read.length();
    auto maxDist = 1.5 * readLen;

    auto k = rapmap::utils::my_mer::k();
    if(remap){
         k = rapmap::utils::my_mer9::k() ;

    }

    //set mac interval
    auto newMaxInterval = SA.size()/std::pow(4,k);
    if(newMaxInterval > maxInterval_ && remap){
        setMaxInterval(SA.size()/std::pow(4,k));
    }else if(remap){
        setMaxInterval(newMaxInterval*2);
    }

    if(remap) setMMPThreshold(mmpThreshold);
    //std::cout << k << " " << int(remap) << "\n";
    //auto k = typename MerT::k() ;
    auto readStartIt = read.begin();
    auto readEndIt = read.end();

    auto rb = read.begin();
    auto re = rb + k;

    uint32_t fwdHit{0};
    uint32_t rcHit{0};

    size_t fwdCov{0};
    size_t rcCov{0};

    bool foundHit = false;
    bool isRev = false;

    //rapmap::utils::my_mer mer;
    //rapmap::utils::my_mer rcMer;

    MerT mer;
    MerT rcMer;

    bool useCoverageCheck{disableNIP_ and strictCheck_};

    // This allows implementing our heurisic for comparing
    // forward and reverse-complement strand matches
    std::vector<KmerDirScore> kmerScores;

    // Where we store the SA intervals for forward and rc hits
    //std::vector<SAIntervalHit> fwdSAInts;
    //std::vector<SAIntervalHit> rcSAInts;

    // Number of nucleotides to skip when encountering a homopolymer k-mer.
    OffsetT homoPolymerSkip = 1; // k / 2;

    // Iterator for k-mer and rc k-mer lookups
    auto merIt = hashEnd_;
    auto rcMerIt = hashEnd_;

    // The position of the k-mer in the read
    size_t pos{0};
    // The position of the next 'N' in the read
    size_t invalidPos{0};

    // Find a hit within the read
    // While we haven't fallen off the end
    while (re <= readEndIt) {

      // Get the k-mer at the current start position.
      // And make sure that it's valid (contains no Ns).
      pos = std::distance(readStartIt, rb);

      // See if this k-mer would contain an N
      // only check if we don't yet know that there are no remaining
      // Ns
      if (invalidPos != std::string::npos) {
        invalidPos = read.find_first_of("nN", pos);
        if (invalidPos <= pos + k) {
          rb = read.begin() + invalidPos + 1;
          re = rb + k;
          continue;
        }
      }

      // If the next k-bases are valid, get the k-mer and
      // reverse complement k-mer
      //mer = rapmap::utils::my_mer(read.c_str() + pos);
      mer = MerT(read.c_str() + pos);
      if (mer.is_homopolymer()) {
        rb += homoPolymerSkip;
        re += homoPolymerSkip;
        /* Walk base-by-base rather than skipping
        // If the first N is within k bases, then this k-mer is invalid
        if (invalidPos < pos + k) {
            // Skip to the k-mer starting at the next position
            // (i.e. right past the N)
            rb = read.begin() + invalidPos + 1;
            re = rb + k;
            // Go to the next iteration of the while loop
            continue;
        }
        */
        continue;
      }
      rcMer = mer.get_reverse_complement();

      // See if we can find this k-mer in the hash
      //if (mer)
      //  std::cout<<"mer null\n";
      //if(mer.Word(0))
        //  std::cout<<"mer.Word null\n";
      //std::cout<<(int)mer.Word(0) <<" mer.Word(0)\n";
      merIt = khash.find(mer.word(0));//get_bits(0, 2 * k));
      rcMerIt = khash.find(rcMer.word(0));//rcMer.get_bits(0, 2 * k));

      // If we can find the k-mer in the hash
      if (merIt != hashEnd_) {
        if (strictCheck_) {
          ++fwdHit;
          // If we also match this k-mer in the rc direction
          if (rcMerIt != hashEnd_) {
            ++rcHit;
            kmerScores.emplace_back(mer, pos, PRESENT, PRESENT);
          } else { // Otherwise it doesn't match in the rc direction
            kmerScores.emplace_back(mer, pos, PRESENT, ABSENT);
          }
        } else { // no strict check
          ++fwdHit;
          if (rcMerIt != hashEnd_) {
            ++rcHit;
          }
        }
      }

      // See if the reverse complement k-mer is in the hash
      if (rcMerIt != hashEnd_) {
        // The original k-mer didn't match in the foward direction
        if (!fwdHit) {
          ++rcHit;
          if (strictCheck_) {
            kmerScores.emplace_back(mer, pos, ABSENT, PRESENT);
          }
        }
      }

      // If we had a hit with either k-mer then we can
      // break out of this loop to look for the next informative position
      if (fwdHit + rcHit > 0) {
        foundHit = true;
        break;
      }
      ++rb;
      ++re;
    }

    // If we went the entire length of the read without finding a hit
    // then we can bail.
    if (!foundHit) {
      return false;
    }

    bool didCheckFwd{false};
    bool sigHit{false};
    // If we had a hit on the forward strand
    if (fwdHit) {
      didCheckFwd = true;
      getSAHits_(saSearcher,
                 read,             // the read
                 rb,               // where to start the search
                 &(merIt->second.interval), // pointer to the search interval
				 merIt->second.lcpLength, // starting lcpLength
                 merIt->second.safeLength, // safe length
                 fwdCov, fwdHit, rcHit, fwdSAInts, kmerScores, false, sigHit, remap);
    }


   //bool checkRC = useCoverageCheck ? (rcHit > 0) : (rcHit >= fwdHit);
    //TODO RC checked enabled
    bool checkRC = useCoverageCheck ? (rcHit > 0) : (rcHit >= fwdHit);
    // If we had a hit on the reverse complement strand
    // if(read=="CAGGCTGGAGTGCAGTGGCACGATCTTGGCTCACTGCAAGCTCCGCCTCCCAGGTTCACGTCATTCCCCTGCCAG" or read=="CTGGCAGGGGAATGACGTGAACCTGGGAGGCGGAGCTTGCAGTGAGCCAAGATCGTGCCACTGCACTCCAGCCTG")

    if (checkRC) {
      rapmap::utils::reverseRead(read, rcBuffer_);
      getSAHits_(saSearcher,
                 rcBuffer_,         // the read
                 rcBuffer_.begin(), // where to start the search
                 nullptr,           // pointer to the search interval
				 0, //starting lcpLength
				 0, //safe length
                 rcCov, rcHit, fwdHit, rcSAInts, kmerScores, true, sigHit, remap);
    }

    // Now, if we *didn't* check the forward strand at first, but we encountered
    // fwd hits
    // while looking at the RC strand, then check the fwd strand now
    //bool checkFwd = useCoverageCheck ? (fwdHit > 0) : (fwdHit >= rcHit);
    //TODO make check forward false
    bool checkFwd = useCoverageCheck ? (fwdHit > 0) : (fwdHit >= rcHit);
    //bool checkFwd = true ;
    if (!didCheckFwd and checkFwd) {
      didCheckFwd = true;
      getSAHits_(saSearcher,
                 read,         // the read
                 read.begin(), // where to start the search
                 nullptr,      // pointer to the search interval
				 0, //starting lcpLength
				 0, //safe length
                 fwdCov, fwdHit, rcHit, fwdSAInts, kmerScores, false, sigHit, remap);
    }


    if(remap){
        if(!sigHit){
            return false;
        }
    }

    /*

    if (strictCheck_) {
      // If we're computing coverage, then we can make use of that info here
      //useCoverageCheck = false;
      if (useCoverageCheck) {
        if (fwdCov > rcCov) {
          rcSAInts.clear();
        } else if (rcCov > fwdCov) {
          fwdSAInts.clear();
        }
      } else { // use the k-mer "spot check"
        // The first two conditions shouldn't happen
        // but I'm just being paranoid here
        if (fwdHit > 0 and rcHit == 0) {
          rcSAInts.clear();
        } else if (rcHit > 0 and fwdHit == 0) {
          fwdSAInts.clear();
        } else {
          std::sort(kmerScores.begin(), kmerScores.end());
          auto e = std::unique(kmerScores.begin(), kmerScores.end());
          // Compute the score for the k-mers we need to
          // test in both the forward and rc directions.
          int32_t fwdScore{0};
          int32_t rcScore{0};
          // For every kmer score structure
          // std::cerr << "[\n";
          for (auto kmsIt = kmerScores.begin(); kmsIt != e;
               ++kmsIt) { //: kmerScores) {
            auto& kms = *kmsIt;
            // If the forward k-mer is untested, then test it
            if (kms.fwdScore == UNTESTED) {
              auto merIt = khash.find(kms.kmer.word(0));//get_bits(0, 2 * k));
              kms.fwdScore = (merIt != hashEnd_) ? PRESENT : ABSENT;
            }
            // accumulate the score
            fwdScore += kms.fwdScore;

            // If the rc k-mer is untested, then test it
            if (kms.rcScore == UNTESTED) {
              rcMer = kms.kmer.get_reverse_complement();
              auto rcMerIt = khash.find(rcMer.word(0));//get_bits(0, 2 * k));
              kms.rcScore = (rcMerIt != hashEnd_) ? PRESENT : ABSENT;
            }
            // accumulate the score
            rcScore += kms.rcScore;
            // kms.print();
            // std::cerr << "\n";
          }
          // std::cerr << "]\n";
          // If the forward score is strictly greater
          // then get rid of the rc hits.
          if (fwdScore > rcScore) {
            rcSAInts.clear();
          } else if (rcScore > fwdScore) {
            // If the rc score is strictly greater
            // get rid of the forward hits
            fwdSAInts.clear();
          }
        }
      }
    }

    // Coverage requirements only make sense if
    // we have disabled NIP skipping.
    if (covReq_ > 0.0 and disableNIP_) {
      double fwdFrac{0.0};
      double rcFrac{0.0};
      if (fwdSAInts.size() > 0) {
        fwdFrac = fwdCov / static_cast<double>(readLen);
        if (fwdFrac < covReq_) {
          fwdSAInts.clear();
        }
      }
      if (rcSAInts.size() > 0) {
        rcFrac = rcCov / static_cast<double>(readLen);
        if (rcFrac < covReq_) {
          rcSAInts.clear();
        }
      }
    }

    auto fwdHitsStart = hits.size();
    // If we had > 1 forward hit
    if (fwdSAInts.size() > 1) {

      //if(remap) std::cout << "\nfwd Hit size " << fwdSAInts.size() << "\n";
        if(!remap){


        if(read == "GTGCAGTGGTGCAATCTCAGCACACTGCAACCTCTGCCGCCCGGGTTCAAGCAATTCTCCTGCCTCAGCCTCCCA" or read == "TGGGAGGCTGAGGCAGGAGAATTGCTTGAACCCGGGCGGCAGAGGTTGCAGTGTGCTGAGATTGCACCACTGCAC"){
          if(merIt != khash.end()){
              std::cout << "\n We are looking at the read \n" ;
              auto tSAInt = merIt->second ;
              for(auto i = tSAInt.interval.begin(); i != tSAInt.interval.end(); ++i){
                  auto iGlobalPos = rmi_->SA[i] ;
                  auto txpID = rmi_->transcriptAtPosition(iGlobalPos);
                  std::cout << rmi_->txpNames[txpID] << "\n";
              }
          }
          std::cout << "\n";
      }

            auto processedHits = rapmap::hit_manager::unionSAHits(
           fwdSAInts, *rmi_, readLen, consistentHits);
                 rapmap::hit_manager::collectHitsSimpleSAGreedy(processedHits, readLen, maxDist,
                                                                hits, false, mateStatus);


          //auto processedHits = rapmap::hit_manager::intersectSAHits(
              //fwdSAInts, *rmi_, readLen, consistentHits);
          //rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist, hits, false, mateStatus);
        }else{
            auto processedHits = rapmap::hit_manager::intersectReSAHits(fwdSAInts, goldenTids,*rmi_, readLen, consistentHits);
            rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist, hits, true, mateStatus);
        }

      //if(remap) std::cout << "\n" << hits.size() << "\n";

    } else if (fwdSAInts.size() == 1) { // only 1 hit!
      auto& saIntervalHit = fwdSAInts.front();
      auto initialSize = hits.size();
      for (OffsetT i = saIntervalHit.begin; i != saIntervalHit.end; ++i) {
        auto globalPos = SA[i];
        auto txpID = rmi_->transcriptAtPosition(globalPos);
        // the offset into this transcript
        auto pos = globalPos - txpStarts[txpID];
        int32_t hitPos = pos - saIntervalHit.queryPos;
        if(!remap){
            hits.emplace_back(txpID, hitPos, true, readLen,saIntervalHit.lcpLength);
            hits.back().mateStatus = mateStatus;
        }else{
            if(std::find(std::begin(goldenTids), std::end(goldenTids), txpID) != std::end(goldenTids)){
                hits.emplace_back(txpID, hitPos, true, readLen,saIntervalHit.lcpLength);
                hits.back().mateStatus = mateStatus;
            }
        }
      }
      // Now sort by transcript ID (then position) and eliminate
      // duplicates
      auto sortStartIt = hits.begin() + initialSize;
      auto sortEndIt = hits.end();
      std::sort(sortStartIt, sortEndIt,
                [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                  if (a.tid == b.tid) {
                    return a.pos < b.pos;
                  } else {
                    return a.tid < b.tid;
                  }
                });
      auto newEnd = std::unique(
          hits.begin() + initialSize, hits.end(),
          [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
            return a.tid == b.tid;
          });
      hits.resize(std::distance(hits.begin(), newEnd));
    }
    auto fwdHitsEnd = hits.size();

    auto rcHitsStart = fwdHitsEnd;
    // If we had > 1 rc hit
    if (rcSAInts.size() > 1) {
        if(!remap){
            auto processedHits = rapmap::hit_manager::unionSAHits(
                      rcSAInts, *rmi_, readLen, consistentHits);
            rapmap::hit_manager::collectHitsSimpleSAGreedy(processedHits, readLen, maxDist,
                                                                           hits, false, mateStatus);
          //auto processedHits = rapmap::hit_manager::intersectSAHits(
              //rcSAInts, *rmi_, readLen, consistentHits);
          //rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist,
                                                   //hits, false,  mateStatus);
        }else{
            auto processedHits = rapmap::hit_manager::intersectReSAHits(rcSAInts, goldenTids,*rmi_, readLen, consistentHits);
            rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist, hits, true, mateStatus);

        }
    } else if (rcSAInts.size() == 1) { // only 1 hit!
      auto& saIntervalHit = rcSAInts.front();
      auto initialSize = hits.size();
      for (OffsetT i = saIntervalHit.begin; i != saIntervalHit.end; ++i) {
        auto globalPos = SA[i];
        auto txpID = rmi_->transcriptAtPosition(globalPos);
        // the offset into this transcript
        auto pos = globalPos - txpStarts[txpID];
        int32_t hitPos = pos - saIntervalHit.queryPos;
        //hits.emplace_back(txpID, hitPos, false, readLen, saIntervalHit.lcpLength);
        //hits.back().mateStatus = mateStatus;
        if(!remap){
            hits.emplace_back(txpID, hitPos, false, readLen,saIntervalHit.lcpLength);
            hits.back().mateStatus = mateStatus;
        }else{
            if(std::find(std::begin(goldenTids), std::end(goldenTids), txpID) != std::end(goldenTids)){
                hits.emplace_back(txpID, hitPos, false, readLen,saIntervalHit.lcpLength);
                hits.back().mateStatus = mateStatus;
            }
        }
      }
      // Now sort by transcript ID (then position) and eliminate
      // duplicates
      auto sortStartIt = hits.begin() + rcHitsStart;
      auto sortEndIt = hits.end();
      std::sort(sortStartIt, sortEndIt,
                [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                  if (a.tid == b.tid) {
                    return a.pos < b.pos;
                  } else {
                    return a.tid < b.tid;
                  }
                });
      auto newEnd = std::unique(
          sortStartIt, sortEndIt,
          [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
            return a.tid == b.tid;
          });
      hits.resize(std::distance(hits.begin(), newEnd));
    }
    auto rcHitsEnd = hits.size();

    // If we had both forward and RC hits, then merge them
    if ((fwdHitsEnd > fwdHitsStart) and (rcHitsEnd > rcHitsStart)) {
      // Merge the forward and reverse hits
      std::inplace_merge(
          hits.begin() + fwdHitsStart, hits.begin() + fwdHitsEnd,
          hits.begin() + rcHitsEnd,
          [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
            return a.tid < b.tid;
          });
       std::sort(hits.begin(),hits.end(),
                 [](const QuasiAlignment& a, const QuasiAlignment& b)->bool {
                    return ((a.tid<b.tid) || (a.tid==b.tid && a.coverage>b.coverage));
            });
      // And get rid of duplicate transcript IDs
      auto newEnd = std::unique(
          hits.begin() + fwdHitsStart, hits.begin() + rcHitsEnd,
          [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
            return a.tid == b.tid;
          });
      hits.resize(std::distance(hits.begin(), newEnd));
    }


    //if(remap) std::cout << "\n" << hits.size() << "\n";
    // Return true if we had any valid hits and false otherwise.

    */
    return foundHit;
  }

private:
  // spot-check k-mers to see if there are forward or rc hits
  template <typename IteratorT>
  inline void
  spotCheck_(//rapmap::utils::my_mer mer,
             MerT mer,
             size_t pos, // the position of the k-mer on the read
             size_t readLen,
             IteratorT* merItPtr,           // nullptr if we haven't checked yet
             IteratorT* complementMerItPtr, // nullptr if we haven't checked yet
             bool isRC, // is this being called from the RC of the read
             uint32_t& strandHits, uint32_t& otherStrandHits,
             std::vector<KmerDirScore>& kmerScores,
             bool remap=false
             ) {
    IteratorT merIt = hashEnd_;
    IteratorT complementMerIt = hashEnd_;
    //auto& khash =(remap)?rmi_->khash9:rmi_->khash;
    auto& khash = rmi_->khash;
    //auto hashEnd_ = khash.end();
    auto k = rapmap::utils::my_mer::k();
    if(remap){
        k = rapmap::utils::my_mer9::k();
    }

    //if(remap) std::cout << "k"

    auto complementMer = mer.get_reverse_complement();

    if (merItPtr == nullptr) {
      // We haven't tested this, so do that here
      merIt = khash.find(mer.word(0));//get_bits(0, 2 * k));
    } else {
      // we already have this
      merIt = *merItPtr;
    }

    if (complementMerItPtr == nullptr) {
      // We haven't tested this, so do that here
      complementMerIt = khash.find(complementMer.word(0));//get_bits(0, 2 * k));
    } else {
      // we already have this
      complementMerIt = *complementMerItPtr;
    }

    HitStatus status{UNTESTED};
    HitStatus complementStatus{UNTESTED};

    if (merIt != hashEnd_) {
      ++strandHits;
      status = PRESENT;
    } else {
      status = ABSENT;
    }
    if (complementMerIt != hashEnd_) {
      ++otherStrandHits;
      complementStatus = PRESENT;
    } else {
      complementStatus = ABSENT;
    }

    HitStatus fwdStatus = isRC ? complementStatus : status;
    HitStatus rcStatus = isRC ? status : complementStatus;

    if (strictCheck_) {
      // If we're on the reverse complement strand, then
      // we have to adjust kmerPos to be with respect to the
      // forward strand.
      if (isRC) {
        auto kp = pos;
        pos = readLen - kp - k;
        mer = complementMer;
      }
      kmerScores.emplace_back(mer, pos, fwdStatus, rcStatus);
    }
  }
  /*
  // Attempts to find the next valid k-mer (a k-mer that doesn't contain an 'N' and is
  // not a homopolymer).  If no such k-mer exists within the read, then it returns false.
  inline bool getNextValidKmer_(std, size_t& pos, rapmap::utils::my_mer& mer) {
      bool validMer = mer.from_chars(read + pos);
      // if this kmer contains an 'N' then validMer is false, else true
  }
  */

  //template <typename MerT>
  inline void getSAHits_(
      SASearcher<RapMapIndexT>& saSearcher, std::string& read,
      std::string::iterator startIt,
      rapmap::utils::SAInterval<OffsetT>* startInterval, uint32_t lcpLengthIn, uint8_t safeLenIn, size_t& cov,
      uint32_t& strandHits, uint32_t& otherStrandHits,
      std::vector<rapmap::utils::SAIntervalHit<OffsetT>>& saInts,
      std::vector<KmerDirScore>& kmerScores,
      bool isRC, // true if read is the reverse complement, false otherwise
      bool& sigHit,
      bool remap=false
      ) {
    using SAIntervalHit = rapmap::utils::SAIntervalHit<OffsetT>;
    //auto& khash = (remap)?rmi_->khash9:rmi_->khash;
    auto& khash = rmi_->khash;

    //auto hashEnd_ = khash.end();
    //
    decltype(hashEnd_)* nullItPtr = nullptr;

    auto readLen = read.length();
    auto readStartIt = read.begin();
    auto readEndIt = read.end();
    OffsetT matchedLen{0};

    auto k = rapmap::utils::my_mer::k();
    if(remap){
        k = rapmap::utils::my_mer9::k() ;
    }

    auto skipOverlapMMP = k - 1;
    auto skipOverlapNIP = k - 1;
    OffsetT homoPolymerSkip = 1;//k / 2;

    auto rb = readStartIt;
    auto re = rb + k;
    OffsetT lb, ub;
    uint32_t lcpLength = lcpLengthIn;
    uint8_t safeLength = safeLenIn;
    size_t invalidPos{0};

    MerT mer, complementMer;
    auto merIt = hashEnd_;
    auto complementMerIt = hashEnd_;
    size_t pos{0};
    size_t sampFactor{1};
    bool lastSearch{false};
    size_t prevMMPEnd{0};
    bool validMer{true};
    
    OffsetT hybridSkip{std::max(size_t(1), readLen / 10)};

    // If we have some place to start that we have already computed
    // then use it.
    bool canSkipSetup{startInterval != nullptr};

    if (canSkipSetup) {
      rb = startIt;
      re = rb + k;
      pos = std::distance(readStartIt, rb);
      invalidPos = pos;
      lb = startInterval->begin();
      ub = startInterval->end();
      
      goto skipSetup;
    }

    while (re <= readEndIt) {
      // The distance from the beginning of the read to the
      // start of the k-mer
      pos = std::distance(readStartIt, rb);
      validMer = mer.from_chars(read.c_str() + pos);
      // Get the next valid k-mer at some position >= pos
      //validMer = getNextValidKmer_(read, pos, mer);
      //if (!validMer) { return; }

      // If this k-mer contains an 'N', then find the position
      // of this character and skip one past it.
      if (!validMer) {
        invalidPos = read.find_first_of("nN", pos);
        // If the first N is within k bases, then this k-mer is invalid
        if (invalidPos < pos + k) {
          // Skip to the k-mer starting at the next position
          // (i.e. right past the N)
          rb = read.begin() + invalidPos + 1;
          re = rb + k;
          // Go to the next iteration of the while loop
          continue;
        }
      }
      // If we got here, we have a k-mer without an 'N'

      // If this is a homopolymer, then skip it
      if (mer.is_homopolymer()) {
        rb += homoPolymerSkip;
        re += homoPolymerSkip;
        continue;
      }

      // If it's not a homopolymer, then get the complement
      // k-mer and query both in the hash.
      complementMer = mer.get_reverse_complement();
      merIt = khash.find(mer.word(0));//get_bits(0, 2 * k));

      if (merIt != hashEnd_) {
        spotCheck_(mer, pos, readLen, &merIt, nullItPtr, isRC, strandHits,
                   otherStrandHits, kmerScores, remap);

        lb = merIt->second.interval.begin();
        ub = merIt->second.interval.end();
        lcpLength = merIt->second.lcpLength ;
        safeLength = merIt->second.safeLength ;
      skipSetup:
        auto oldlb = lb;
        auto oldub = ub;


        if(readStartIt == startIt){

            lb = std::max(static_cast<OffsetT>(0), lb - 1);
            std::tie(lb, ub, matchedLen) =
                saSearcher.extendSearchNaive(lb, ub, k, rb, readEndIt);

            // **Note** : If we don't match the whole read, why go back and
            // re-do the extendSafe?  We know how long the matchedLen is, so
            // we should be able to just reduce it to min(matchedLen, safeLength)
            // and proceed without re-checking, right?

            if(matchedLen == readLen){
                matchedLen = readLen;
            } else {
                if(safeLength==k) { safeLength=k+1; }
                /** Assuming there isn't an OBO here, I think we can
                    replace the call to extendSafe in this block by
                    the code immediately below **/
                
		lb = oldlb;
		ub = oldub;
		auto newExtend =  std::min((uint32_t)matchedLen,(uint32_t) safeLength); 
                if (newExtend > k+1) {
                  matchedLen = newExtend;
                } else {
                  matchedLen = k+hybridSkip;
                }
                
                /**auto newExtend =  saSearcher.extendSafe(lb, ub, k, rb, readEndIt, safeLength);
                if(newExtend > k){
                    matchedLen = newExtend ;
                }else{
                    matchedLen = k+hybridSkip;
                }**/
            }
        } else {
            if(safeLength==k)
                safeLength=k+1;
            auto newExtend =  saSearcher.extendSafe(lb, ub, k, rb, readEndIt,safeLength );
            if(newExtend > k){
                matchedLen = newExtend;
            }else{
                matchedLen = k+hybridSkip;
            }
        }


        OffsetT diff = ub - lb;
        if (ub > lb and diff < maxInterval_) {
          uint32_t queryStart =
              static_cast<uint32_t>(std::distance(readStartIt, rb));

          auto bigk = rapmap::utils::my_mer::k() ;
          if(remap){
              if(matchedLen > mmpThreshold_){
                //std::cout << "k: " << k << " MMP Length: "<<matchedLen << "\n";
                saInts.emplace_back(lb, ub, matchedLen, queryStart, lcpLength,isRC);
              }
          }else{
                saInts.emplace_back(lb, ub, matchedLen, queryStart, lcpLength,isRC);

          }


        if(saInts.size() > 0 && remap && !sigHit){
            //std::cout << "\n remapped \n";
            sigHit=true;
        }

          size_t matchOffset = std::distance(readStartIt, rb);
          size_t correction = 0;

          // NOTE: prevMMPEnd points 1 position past the last *match* of the
          // previous MMP (i.e. it points to the *first mismatch*).  This is
          // why we ignore the case where prevMMPEnd == matchOffset, and why
          // we don't have to add 1 to correction.
          if (prevMMPEnd > matchOffset) {
            correction = prevMMPEnd - matchOffset;
          }
          // Update the coverage and position of the last MMP match
          cov += (matchedLen - correction);
          prevMMPEnd = matchOffset + matchedLen;

          // If we didn't end the match b/c we exhausted the query
          // test the mismatching k-mer in the other strand
          if (rb + matchedLen < readEndIt) {
            uint32_t kmerPos = static_cast<uint32_t>(
                std::distance(readStartIt, rb + matchedLen - skipOverlapMMP));
            bool validNucs = mer.from_chars(read.c_str() + kmerPos);
            if (validNucs) {
              /*
              // since the MMP *ended* before the end of the read, we assume
              // that the k-mer one past the MMP is a mismatch (i.e is ABSENT)
              // we avoid looking it up in spotCheck_ by simply passing a pointer
              // to the end of the k-mer hash, which will treat this mer as ABSENT.
              auto endItPtr = &hashEnd_;
              */
              // Even though the MMP *ended* before the end of the read, we're still
              // going to check the mismatching k-mer in both directions to ensure that
              // it doesn't appear somewhere else in the forward direction
              spotCheck_(mer, kmerPos, readLen, nullItPtr, nullItPtr, isRC,
                         strandHits, otherStrandHits, kmerScores, remap);
            }
          } // we didn't end the search by falling off the end
        // This hit was worth recording --- occurred fewer then maxInterval_
            // times
        }else{
            //if(diff > maxInterval_ && remap) std::cout<<"\nI got chucked out for diff: " <<diff<<" maxInterval_: "<<maxInterval_<<"\n";
        }

        // If we've previously declared that the search that just occurred was
        // our last, then we're done!
        if (lastSearch) {
          return;
        }

        // Otherwise, figure out how we should continue the search.
        auto mismatchIt = rb + matchedLen;
        // If we reached the end of the read, then we're done.
        if (mismatchIt >= readEndIt) {
          return;
        }

        auto remainingDistance = std::distance(mismatchIt, readEndIt);
        auto lce = disableNIP_ ? matchedLen
                               : saSearcher.lce(lb, ub - 1, matchedLen,
                                                remainingDistance);

        // Where we would jump if we just used the MMP
        auto skipMatch = mismatchIt - skipOverlapMMP;
        // if (skipMatch + k )
        // Where we would jump if we used the LCE
        auto skipLCE = rb + lce - skipOverlapNIP;
        // Pick the maximum of the two
        auto maxSkip = std::max(skipMatch, skipLCE);
        rb = rb+matchedLen-k+1;

        // If NIP skipping is *enabled*, and we got to the current position
        // by doing an LCE query, then we allow ourselves to *double check*
        // by querying the last k-mer in the read.
        // Otherwise, we just take the skip we're given.
        if (!disableNIP_ and (lce > matchedLen)) {
          if (readLen > k) {
            rb = std::min(readEndIt - k, rb);
          }
        }

        re = rb + k;

        // If the search ends at the end of the read, then
        // set the flag that we need not try after this.
        if (re == readEndIt) {
          lastSearch = true;
        }

      } else { // If we couldn't match this k-mer, move on to the next.

        // &merIt should point to the end of the k-mer hash,
        // complementMerItPtr is null because we want to spot-check the complement k-mer.
        spotCheck_(mer, pos, readLen, &merIt, nullItPtr, isRC, strandHits,
                   otherStrandHits, kmerScores, remap);
        rb += sampFactor;
        re = rb + k;
      }
    }
  }

  RapMapIndexT* rmi_;
  decltype(rmi_->khash.end()) hashEnd_;
  bool disableNIP_;
  double covReq_;
  OffsetT maxInterval_;
  uint32_t mmpThreshold_ ;
  bool strictCheck_;
  std::string rcBuffer_;
};

#endif // SA_COLLECTOR_PAIR_HPP
