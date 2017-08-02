//jello

#ifndef __SELECTIVE_ALIGNMENT_HPP__
#define __SELECTIVE_ALIGNMENT_HPP__

#include "RapMapSAIndex.hpp"
#include "RapMapUtils.hpp"
#include "SASearcher.hpp"
#include "EditDistance.hpp"
#include "edlib.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <set>

//#include "dtl/dtl.hpp"

static constexpr int8_t rc_table[128] = {
      78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 15
      78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 31
      78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 787
      78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 63
      78, 84, 78, 71, 78,  78,  78, 67, 78, 78, 78, 78,  78, 78, 78, 78, // 79
      78, 78,  78, 78,  65, 65, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 95
      78, 84, 78, 71, 78,  78,  78, 67, 78, 78, 78, 78,  78, 78, 78, 78, // 101
      78, 78,  78, 78,  65, 65, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78  // 127
  };

/*
template <typename OffsetT>
int32_t SECollector::hammingDist(QuasiAlignment& qa, std::string& read, std::string& txpSeq,  OffsetT maxOffset, int maxDist){
       int32_t hamming = 0;
       auto minOffset = (qa.pos < 0) ? -(qa.pos) : 0;
       hamming += minOffset;
       //int maxOffset = std::min(static_cast<int>(read.size()), static_cast<int>(transcript.RefLength - qa.pos));
       hamming += (read.size() - maxOffset);
       //const char* txpSeq = transcript.Sequence();
       //const char* readSeq{nullptr};
       auto readLen = read.length();
       if (qa.fwd) {
       for (int j=minOffset; j<maxOffset; ++j) {
            hamming += (read[j]!=txpSeq[qa.pos+j]);
                if (hamming > maxDist) { return hamming;  }
            }
        } else {
            for (int j=minOffset; j<maxOffset; ++j) {
             hamming += ((char)rc_table[(uint8_t)read[(readLen-j)-1]] != txpSeq[qa.pos+j]);
                if (hamming > maxDist) { return hamming;  }
            }
        }
        return hamming;

}
*/



template <typename RapMapIndexT> class SECollector{
public:
  using OffsetT = typename RapMapIndexT::IndexType;

  /** Construct an SECollector given an index **/
    SECollector(RapMapIndexT* rmi)
        : rmi_(rmi) {}

    //int32_t hammingDist(QuasiAlignment& q, std::string& read, std::string& seq,  Offset trancriptLen, int maxDist);

   template <typename ReadStructT>
       bool compute(ReadStructT& readT,
                  std::vector<rapmap::utils::QuasiAlignment>& hits,
                  int32_t editThreshold
                  ) {
           using QuasiAlignment = rapmap::utils::QuasiAlignment;
           using MateStatus = rapmap::utils::MateStatus;
           using kmerVal = rapmap::utils::kmerVal<OffsetT>;
	  //go through the lcp length of each SA Interval for each hit
	  //figure out if it needs to be changed


           if(hits.size() == 0){
               return false ;
           }

          //Mapping information
          //Loaed
          //std::cout << "\n HERE \n" ;
          std::set<int32_t> tidset ;

	  auto& txpStarts = rmi_->txpOffsets ;
	  auto& SA = rmi_->SA;
	  auto& concatText = rmi_->seq;
	  auto& txpLens = rmi_->txpLens ;
          auto& khash = rmi_->khash ;

          auto k = rapmap::utils::my_mer::k();


          auto& read = readT.seq;

          rapmap::utils::my_mer mer;

          //std::vector<SAIntervalHit> fwdSAInts ;
          //std::vector<SAIntervalHit> rcSAInts ;
          kmerVal tSAInt ;


          //start hit information
	  auto& startHit = hits.front();
	  auto lcpLength = startHit.lcpLength ;
	  auto readLen = read.length();

          //uint8_t editThreshold = readLen/2 ;
	  int startEditDistance = 0;

          bool skipLCPOpt{false};
          bool currentValid{false};
          std::string currentKmer;

        //TODO check if lcp is really returning same sequences
          std::string firsttidString ;

    /*
            auto& readName = readT.name;
            // If the read name contains multiple space-separated parts,
            // print only the first
            size_t splitPos = readName.find(' ');
            if (splitPos < readName.length()) {
                readName[splitPos] = '\0';
            } else {
                splitPos = readName.length();
            }

            // trim /1 from the pe read
            if (splitPos > 2 and readName[splitPos - 2] == '/') {
                readName[splitPos - 2] = '\0';
            }
    */

	  // There are two possible ways to go
	  // from here
	  // We have a large lcp so we will calculate alignment only once
	  // and gonna copy for rest of the things
	  // Nevertheless we have to align the first read to the transcript

	  // Align the startHit any way
	  // get transcript id
	  // sequence and read sequence
	  // and position

            if(!startHit.toAlign){
                  uint32_t txpID = startHit.tid ;
                  int32_t pos = startHit.pos;
                  int32_t startOffset, endOffset ;
                  OffsetT globalPos = txpStarts[txpID];
                  OffsetT thisTxpLen = txpLens[txpID];

                  auto overHangLeft = (pos < 0)?-(pos):0;
                  auto overHangRight = (pos+readLen > thisTxpLen)?(pos+readLen-thisTxpLen):0;

                  globalPos = (overHangLeft == 0)?(pos+globalPos):globalPos;
                  auto extend = (overHangLeft > 0)?(readLen - overHangLeft):readLen ;
                  extend = (overHangRight > 0)?(extend-overHangRight):extend;

                  /* ROB: get rid of the copy */
                  //auto thisTxpSeq = concatText.substr(globalPos, extend);
                  const char* thisTxpSeq = concatText.data() + globalPos;
                  int thisTargetLen = extend;


                  /* ROB : slight interface change */
                  //EdlibAlignResult startEdlibResult;
                  if(startHit.fwd){
                    ae_(read.c_str(), read.length(), thisTxpSeq, thisTargetLen, edlibNewAlignConfig((editThreshold+1)*3, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));
                  }else{
                    auto revRead = rapmap::utils::reverseComplement(read);
                    ae_(revRead.c_str(), read.length(), thisTxpSeq, thisTargetLen, edlibNewAlignConfig((editThreshold+1)*3, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));
                  }
                  auto& startEdlibResult = ae_.result();

                  startHit.editD = startEdlibResult.editDistance;

                  startEditDistance = startEdlibResult.editDistance;

                  //@debug purpose
                  auto readName = rapmap::utils::getReadName(readT) ;

                  /*
                  if(readName == "15006_1_6022_314_181"){
                      std::cout << "\nIn selective Alignment start  edit distance: "<< startEditDistance << " length: "<< thisTxpLen << " pos " << pos << " fwd "<<  startHit.fwd << "\n";
                  }*/

                  if(startEditDistance <= editThreshold ){
                    startHit.toAlign = true;
                    /* ROB: No CIGAR right now */
                    /*
                    char* cigar_ = edlibAlignmentToCigar(startEdlibResult.alignment, startEdlibResult.alignmentLength, EDLIB_CIGAR_STANDARD);
                    std::string cigar(cigar_);
                    startHit.cigar = cigar ;
                    */
                 }else{
                    startHit.toAlign = false ;
                 }

                /* Take the kmer from the transcript */
                currentKmer = concatText.substr(globalPos, k);
                firsttidString = concatText.substr(globalPos, 75);


                if(currentKmer.length() == k and currentKmer.find_first_of('$') == std::string::npos){
                    mer = currentKmer;
                    currentValid= true ;
                }else{
                    currentValid = false ;
                }
          }

          std::map<int,OffsetT> tidPos ;
          if(currentValid){
              auto bits = mer.word(0);
              auto hashIt = khash.find(bits);
              if(hashIt != khash.end()){
                  //std::cout << "\n Found \n";
                  tSAInt = hashIt->second;
                  lcpLength = hashIt->second.lcpLength ;
                  for(OffsetT i = tSAInt.interval.begin() ; i != tSAInt.interval.end(); ++i ){
                      auto iGlobalPos = SA[i];
                      auto txpID = rmi_->transcriptAtPosition(iGlobalPos);
                      tidset.insert(txpID);
                      tidPos[txpID] = iGlobalPos;
                }
              }else{
                  skipLCPOpt = true ;
              }

          }else{
              skipLCPOpt = true ;
          }




          if(hits.size() > 1){
            for(auto hitsIt= hits.begin()+1 ; hitsIt != hits.end() ; ++hitsIt){
                uint32_t txpID = hitsIt->tid ;
                auto search = tidset.find(txpID);
                int32_t pos = hitsIt->pos;
                int32_t startOffset, endOffset ;
                OffsetT globalPos = txpStarts[txpID];
                OffsetT thisTxpLen = txpLens[txpID];
                auto overHangLeft = (pos < 0)?-(pos):0;
                auto overHangRight = (pos+readLen > thisTxpLen)?(pos+readLen-thisTxpLen):0;

                globalPos = (overHangLeft == 0)?(pos+globalPos):globalPos;

                if(hitsIt->toAlign){
                    continue;
                }
                if ( (!skipLCPOpt) && (search != tidset.end()) && (lcpLength >= readLen) && (tidPos[txpID] == globalPos)){
                    if(startEditDistance <= editThreshold){
                                hitsIt->editD = startEditDistance;
                                hitsIt->toAlign = true;
                                //std::cout << " LCP Length " << lcpLength << "\n";

                                /* ROB: No CIGAR right now */
                                //hitsIt->cigar = startHit.cigar;
                          }else{
                               hitsIt->editD = startEditDistance;
                                hitsIt->toAlign = false ;
                          }
                }else{
                      auto extend = (overHangLeft > 0)?(readLen - overHangLeft):readLen ;
                      extend = (overHangRight > 0)?(extend-overHangRight):extend;

                      /* ROB: get rid of the copy */
                      //auto thisTxpSeq = concatText.substr(globalPos, extend);
                      const char* thisTxpSeq = concatText.data() + globalPos;
                      int thisTargetLen = extend;

                      //auto thisEditDistance = edit_distance(read, thisTxpSeq, 50) ;
                      /* ROB : slight interface change */
                      //EdlibAlignResult thisEdlibResult;
                      if(hitsIt->fwd){
                        ae_(read.c_str(), read.length(), thisTxpSeq, thisTargetLen, edlibNewAlignConfig(1+editThreshold*3, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));
                      }else{
                        auto revRead = rapmap::utils::reverseComplement(read);
                        ae_(revRead.c_str(), read.length(), thisTxpSeq, thisTargetLen, edlibNewAlignConfig(1+editThreshold*3, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));
                      }
                      auto& thisEdlibResult = ae_.result();
                      auto thisEditDistance = thisEdlibResult.editDistance ;


                      if(hitsIt->toAlign and thisEditDistance!=0){
                          std::cout<<"WTF\n";
                      }
                      if(thisEditDistance <= editThreshold){
                              //selectedHits.emplace_back(txpID,pos,startHit.fwd,hitsIt->readLen, thisEditDistance,"II");
                              hitsIt->editD = thisEditDistance;
                              hitsIt->toAlign = true;
                              /* ROB: No CIGAR right now */
                              /*
                              char* cigar_ = edlibAlignmentToCigar(thisEdlibResult.alignment, thisEdlibResult.alignmentLength, EDLIB_CIGAR_STANDARD);
                              std::string cigar(cigar_);
                              hitsIt->cigar = cigar;
                              */
                      } else {
                        hitsIt->editD = thisEditDistance;
                        hitsIt->toAlign = false;
                      }
                }

            }
        }

          //std::cout << "\n DONE \n" ;
          return true ;
   }

   template <typename ReadStructT>
       void operator()(ReadStructT& leftReadT,ReadStructT& rightReadT,
                  std::vector<rapmap::utils::QuasiAlignment>& hits,
                  int32_t editThreshold
               ) {
           auto leftHits = hits;
           auto rightHits = hits;
           //compute(leftReadT,hits,editThreshold);
           //for(auto hit: hits)
             //  std::cout<< hit.editD << "\n";
           for(auto& hit: rightHits){
                hit.pos = hit.matePos;
                hit.fwd = hit.mateIsFwd;
                hit.toAlign = hit.mateToAlign;
                hit.editD = hit.mateEditD;
           }

           compute(rightReadT,rightHits,editThreshold);
           compute(leftReadT,leftHits,editThreshold);

           bool edit = false; bool edit_r = true;
           for(int i=0; i<leftHits.size(); i++) {
               if(leftHits[i].editD==-1)
                    leftHits[i].editD = (editThreshold+1)*3;
                if(rightHits[i].editD==-1)
                    rightHits[i].editD = (editThreshold+1)*3;

                hits[i].editD = leftHits[i].editD + rightHits[i].editD;
                hits[i].toAlign =  hits[i].editD <= 2*editThreshold;

           }
       }

    template <typename ReadStructT>
	void operator()(ReadStructT& leftReadT,
			std::vector<rapmap::utils::QuasiAlignment>& hits,
			int32_t editThreshold) { 
			compute(leftReadT,hits,editThreshold);
	}


private:
  RapMapIndexT* rmi_ ;
  AlignerEngine ae_;

};

#endif
