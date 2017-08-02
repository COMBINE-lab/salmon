#include "HitManagerSA.hpp"
#include "BooMap.hpp"
#include "FastxParser.hpp"
#include "FrugalBooMap.hpp"
#include "HitManager.hpp"
#include "SelectiveAlignment.hpp"

#include "edlib.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <set>

#include <type_traits>

namespace rapmap {
namespace hit_manager_sa {

/*
template <typename RapMapIndexT, typename OffsetT>
int32_t alignRead(std::string& read,
        RapMapIndexT& rmi,
        OffsetT& pos){

    AlignerEngine ae_ ;
    auto& txpStarts = rmi_->txpOffsets ;
    auto& SA = rmi_->SA;
    auto& concatText = rmi_->seq;
    auto& txpLens = rmi_->txpLens ;
    auto& khash = rmi_->khash ;


}*/



template <typename RapMapIndexT>
bool mergeLeftRightMap(fastx_parser::ReadPair& rpair, SAHitMap& leftMap,
                       SAHitMap& rightMap,
                       std::vector<QuasiAlignment>& jointHits,
                       uint32_t editDistance, RapMapIndexT& rmi, uint32_t maxInsertSize_) {
  // using OffsetT = typename RapMapIndexT::IndexType;
  typedef struct covInfo {
  public:
    covInfo() {
      cov = 0;
      endPos = 0;
      rc = false;
      mmp = false;
    };
    // covInfo(SATxpQueryPos s,int c, int e,bool r) : saTxpQp(s), cov(c),
    // endPos(e) {rc = r;};
    covInfo(int c, int e, bool r, bool m) : cov(c), endPos(e), rc(r), mmp(m){};
    // SATxpQueryPos saTxpQp ;
    void operator=(const covInfo& CV) {
      cov = CV.cov;
      endPos = CV.endPos;
      rc = CV.rc;
      mmp = CV.mmp;
    }

    int cov;
    int endPos;
    bool rc;
    bool mmp;
  };

  constexpr const int32_t signedZero{0};
  //uint32_t maxInsertSize_{1000};
  bool foundHit{false};

  SECollector<RapMapIndexT> hitSECollector(&rmi);
  std::vector<QuasiAlignment> hits;

  std::set<uint32_t> commonTids;
  std::set<uint32_t> leftTids;
  for (auto& phLeft : leftMap) {
    auto leftTid = phLeft.first;
    leftTids.insert(leftTid);
  }
  for (auto& phRight : rightMap) {
    auto rightTid = phRight.first;
    if (leftTids.find(rightTid) != leftTids.end()) {
      commonTids.insert(rightTid);
    }
  }

  // This seems redundant --- if there are
  // no commonTids, then foundHit will remain the
  // default value (false) listed above.
  //if (commonTids.size() > 0) {

    // for each common tid check
    // if both the end tell us
    // something
    for (auto tid : commonTids) {

      if (leftMap[tid].active and rightMap[tid].active) {

        auto& lefttqvec = leftMap[tid].tqvec;
        auto& righttqvec = rightMap[tid].tqvec;

        std::sort(lefttqvec.begin(), lefttqvec.end(),
                  [](const SATxpQueryPos& a, const SATxpQueryPos& b) -> bool {
                    return a.pos < b.pos;
                  });

        std::sort(righttqvec.begin(), righttqvec.end(),
                  [](const SATxpQueryPos& a, const SATxpQueryPos& b) -> bool {
                    return a.pos < b.pos;
                  });

        // look at the things in pair
        // Make possible pairs
        // Discard pairs which are
        // impossible in the sense
        // that they are far(1000 bp) apart

        std::map<int32_t, covInfo> leftHitCov;
        {
          auto it = lefttqvec.begin();
          while (it != lefttqvec.end()) {
            int numOfCov = 0;
            int32_t hitKey = it->pos - it->queryPos;
            if (leftHitCov.count(hitKey) == 0) {
              covInfo cInfo(it->matchedLen, it->pos + it->matchedLen,
                            it->queryRC,
                            it->matchedLen == rpair.first.seq.length());
              leftHitCov[hitKey] = cInfo;
            } else {
              auto& curr = leftHitCov[hitKey];
              if (curr.endPos > it->pos) {
                curr.cov += it->pos + it->matchedLen - curr.endPos;
                curr.endPos = it->pos + it->matchedLen;
              } else {
                curr.cov += it->matchedLen;
                curr.endPos = it->pos + it->matchedLen;
              }
            }
            ++it;
          }
        }

        std::map<int32_t, covInfo> rightHitCov;
        {
          auto it = righttqvec.begin();
          while (it != righttqvec.end()) {
            int numOfCov = 0;

            int32_t hitKey = it->pos - it->queryPos;
            if (rightHitCov.count(hitKey) == 0) {
              covInfo cInfo(it->matchedLen, it->pos + it->matchedLen,
                            it->queryRC,
                            it->matchedLen == rpair.second.seq.length());
              rightHitCov[hitKey] = cInfo;
            } else {
              auto& curr = rightHitCov[hitKey];
              if (curr.endPos > it->pos) {
                curr.cov += (it->pos + it->matchedLen - curr.endPos);
                curr.endPos = it->pos + it->matchedLen;
              } else {
                curr.cov += it->matchedLen;
                curr.endPos = it->pos + it->matchedLen;
              }
            }
            ++it;
          }
        }

        int maxCovScore = 0;
        std::vector<std::pair<int32_t, int32_t>> validPairs;
        bool leftFirst =
            leftHitCov.cbegin()->first - rightHitCov.cbegin()->first < 0;
        if (leftFirst) {
          for (auto& leftTQPos : leftHitCov) {
            int32_t hitPos1 = leftTQPos.first;  // pos - leftTQPos.queryPos ;
            if (std::abs(hitPos1 - rightHitCov.cbegin()->first) >
                maxInsertSize_) {
              continue;
            }
            for (auto& rightTQPos : rightHitCov) {

              // Check for the distance
              // go with it otherwise
              // go with something else
              int32_t hitPos2 = rightTQPos.first; // pos - rightTQPos.queryPos ;

              if (std::abs(hitPos1 - hitPos2) < maxInsertSize_) {
                validPairs.push_back(std::make_pair(hitPos1, hitPos2));
                auto score = leftTQPos.second.cov + rightTQPos.second.cov;
                if (maxCovScore < score) { maxCovScore = score; }
              } else {
                break;
              }
            }
          } // end for*/

        } else {

          for (auto& rightTQPos : rightHitCov) {
            int32_t hitPos2 = rightTQPos.first; // pos - rightTQPos.queryPos ;
            if (std::abs(hitPos2 - leftHitCov.cbegin()->first) >
                maxInsertSize_) {
              continue;
            }
            for (auto& leftTQPos : leftHitCov) {
              // Check for the distance
              // go with it otherwise
              // go with something else
              int32_t hitPos1 = leftTQPos.first;  // pos - leftTQPos.queryPos ;

              if (std::abs(hitPos1 - hitPos2) < maxInsertSize_) {
                validPairs.push_back(std::make_pair(hitPos1, hitPos2));
                auto score = leftTQPos.second.cov + rightTQPos.second.cov;
                if (maxCovScore < score) { maxCovScore = score; }
              } else {
                break;
              }
            }
          } // end for*/
        }

        // Now I have valid pairs so I can directly make QuasiAlignment vectors
        //
        if (validPairs.size() > 0) {
          foundHit = true;
        }

        /*
        int maxCovScore = 0;
        for(auto& vp : validPairs){

          auto score = leftHitCov[(vp.first)].cov + rightHitCov[(vp.second)].cov
        ;
          if(maxCovScore < score){
            maxCovScore = score;
          }
        }
        */
        hits.clear();
        for (auto& vp : validPairs) {
          // I will try to make QuasiAlignment objects
          // real alignments and edit distance are not present
          // yet

          int32_t hitPos1 = vp.first;  // vp.first.pos - vp.first.queryPos ;
          int32_t hitPos2 = vp.second; // vp.second.pos - vp.second.queryPos ;

          if (leftHitCov[hitPos1].cov + rightHitCov[hitPos2].cov ==
              maxCovScore) {
            int32_t startRead1 = std::max(hitPos1, signedZero);
            int32_t startRead2 = std::max(hitPos2, signedZero);

            bool read1First{(startRead1 < startRead2)};
            int32_t fragStartPos = read1First ? startRead1 : startRead2;
            int32_t fragEndPos = read1First
                                     ? (startRead2 + rpair.second.seq.length())
                                     : (startRead1 + rpair.first.seq.length());
            uint32_t fragLen = fragEndPos - fragStartPos;
            jointHits.emplace_back(tid, hitPos1, !leftHitCov[vp.first].rc,
                                   rpair.first.seq.length(),
                                   0, // leftHitCov[vp.first].saTxpQp.lcpLength,
                                   fragLen, true);
            auto& qaln = jointHits.back();
            // qaln.fwd = !leftHitCov[vp.first].rc;
            qaln.mateLen = rpair.second.seq.length();
            qaln.matePos = hitPos2;
            qaln.mateIsFwd =
                !(rightHitCov[vp.second].rc); //! vp.second.queryRC ;
            qaln.mateStatus = rapmap::utils::MateStatus::PAIRED_END_PAIRED;
            if (leftHitCov[hitPos1].mmp) {
              qaln.toAlign = true;
              qaln.editD = 0;
            }
            if (rightHitCov[hitPos2].mmp) {
              qaln.mateToAlign = true;
              qaln.mateEditD = 0;
            }
            break;
          }
        }

      } // end if
    }   // end for

    // This seems redundant
    //} else {
    //foundHit = false;
    //}

  // since we use an ordered set, these should be sorted.
  /*
  if(jointHits.size()>0){
    std::sort(jointHits.begin(), jointHits.end(),
              [](const QuasiAlignment& a, const QuasiAlignment& b)-> bool {
                return a.tid < b.tid;
              } );
  }
  */
  return foundHit;
}


//template <typename RapMapIndexT>
bool mergeMap(fastx_parser::ReadSeq& rp, SAHitMap& leftMap,
                       std::vector<QuasiAlignment>& jointHits) {

	for (auto& ph : leftMap) {
		auto tid = ph.first;
		//sort transcripts with respect to
		//positions
		std::sort(ph.second.tqvec.begin(),
			ph.second.tqvec.end(),
			[](const SATxpQueryPos& a, const SATxpQueryPos& b)-> bool {
				return a.pos < b.pos;
			} );

		std::map<int32_t , int> hitCov ;
		std::map<int32_t , int> hitEndPos ;
		auto it = ph.second.tqvec.begin();
		while(it != ph.second.tqvec.end()){
		    int numOfCov = 0;
		    int32_t hitKey = it->pos - it->queryPos ;
		    if(hitCov.count(hitKey) == 0){
			hitCov[hitKey] = it->matchedLen;
			hitEndPos[hitKey] = it->pos + it->matchedLen;
		    }
		    else{
			if(hitEndPos[hitKey] > it->pos){
			    hitCov[hitKey] += ( it->pos + it->matchedLen - hitEndPos[hitKey] );

			}else{
			    hitCov[hitKey] += it->matchedLen;
			    hitEndPos[hitKey] = it->pos + it->matchedLen ;
			}
		    }
		    ++it ;
		}
		//std::cout << "\n outside while \n";

		auto minPosIt = std::min_element(ph.second.tqvec.begin(),
				ph.second.tqvec.end(),
				[](const SATxpQueryPos& a, const SATxpQueryPos& b) -> bool {
				    return a.pos < b.pos;
				});

		using pair_type = decltype(hitCov)::value_type;
		auto maxCov = std::max_element(hitCov.begin(),hitCov.end(),
			[](const pair_type& p1, const pair_type& p2) -> bool {
				return p1.second < p2.second ;
			});

		bool hitRC = ph.second.tqvec[0].queryRC;
		int32_t hitPos = maxCov->first;
		bool isFwd = !hitRC;
		jointHits.emplace_back(tid, hitPos, isFwd, rp.seq.length(), 0,0,true);
	
                auto& qaln = jointHits.back();
                // qaln.fwd = !leftHitCov[vp.first].rc;
                qaln.mateStatus = rapmap::utils::MateStatus::SINGLE_END;
                if (maxCov->second== rp.seq.length()) {
                  qaln.toAlign = true;
                  qaln.editD = 0;
                }
	}
	if(jointHits.size()>0)
		return true;
	else
		return false;
}



template <typename RapMapIndexT>
bool mergeLeftRightSAInts(
    fastx_parser::ReadPair& rpair, bool lhp, bool rhp,
    std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& leftFwdSAInts,
    std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& leftRcSAInts,
    std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>&
        rightFwdSAInts,
    std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& rightRcSAInts,
    std::vector<QuasiAlignment>& jointHits, RapMapIndexT& rmi, bool maxNumHits,
    bool consistentHits, rapmap::utils::HitCounters& hctr,
    uint32_t editDistance, uint32_t maxInsertSize_) {

  using OffsetT = typename RapMapIndexT::IndexType;
  AlignerEngine ae_;

  auto& SA = rmi.SA;
  auto& txpStarts = rmi.txpOffsets;
  auto& txpIDs = rmi.positionIDs;

  //uint32_t maxInsertSize_{1000};

  bool foundHit{false};

  SAHitMap leftFwdMap;
  SAHitMap leftRcMap;
  SAHitMap rightFwdMap;
  SAHitMap rightRcMap;

  if (leftFwdSAInts.size() > 0 or leftRcSAInts.size() > 0) {
    if (leftFwdSAInts.size() > 1) {
      leftFwdMap = rapmap::hit_manager::unionSAHits(
          leftFwdSAInts, rmi, rpair.first.seq.length(), consistentHits);
    } else if (leftFwdSAInts.size() > 0) {
      auto inHit = leftFwdSAInts.front();
      for (auto i = inHit.begin; i < inHit.end; ++i) {
        auto globalPos = SA[i];
        auto tid = rmi.transcriptAtPosition(globalPos);
        auto txpPos = globalPos - txpStarts[tid];
        leftFwdMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos,
                                           inHit.lcpLength, inHit.queryRC);
        leftFwdMap[tid].active = true;
        leftFwdMap[tid].tqvec.back().matchedLen = inHit.len;
      }
    }

    if (leftRcSAInts.size() > 1) {
      leftRcMap = rapmap::hit_manager::unionSAHits(
          leftRcSAInts, rmi, rpair.first.seq.length(), consistentHits);
    } else if (leftRcSAInts.size() > 0) {
      auto inHit = leftRcSAInts.front();
      for (auto i = inHit.begin; i < inHit.end; ++i) {
        auto globalPos = SA[i];
        auto tid = rmi.transcriptAtPosition(globalPos);
        auto txpPos = globalPos - txpStarts[tid];
        leftRcMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos,
                                          inHit.lcpLength, inHit.queryRC);
        leftRcMap[tid].active = true;
        leftRcMap[tid].tqvec.back().matchedLen = inHit.len;
      }
    }
  }

  if (rightFwdSAInts.size() > 0 or rightRcSAInts.size() > 0) {
    if (rightFwdSAInts.size() > 1) {
      rightFwdMap = rapmap::hit_manager::unionSAHits(
          rightFwdSAInts, rmi, rpair.second.seq.length(), consistentHits);
    } else if (rightFwdSAInts.size() > 0) {
      auto inHit = rightFwdSAInts.front();
      for (auto i = inHit.begin; i < inHit.end; ++i) {
        auto globalPos = SA[i];
        auto tid = rmi.transcriptAtPosition(globalPos);
        auto txpPos = globalPos - txpStarts[tid];
        rightFwdMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos,
                                            inHit.lcpLength, inHit.queryRC);
        rightFwdMap[tid].active = true;
        rightFwdMap[tid].tqvec.back().matchedLen = inHit.len;
      }
    }
    if (rightRcSAInts.size() > 1) {
      rightRcMap = rapmap::hit_manager::unionSAHits(
          rightRcSAInts, rmi, rpair.second.seq.length(), consistentHits);
    } else if (rightRcSAInts.size() > 0) {
      auto inHit = rightRcSAInts.front();
      for (auto i = inHit.begin; i < inHit.end; ++i) {
        auto globalPos = SA[i];
        auto tid = rmi.transcriptAtPosition(globalPos);
        auto txpPos = globalPos - txpStarts[tid];
        rightRcMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos,
                                           inHit.lcpLength, inHit.queryRC);
        rightRcMap[tid].active = true;
        rightRcMap[tid].tqvec.back().matchedLen = inHit.len;
      }
    }
  }

  //jointHits in two different orientations
  //std::vector<QuasiAlignment> fwdRcHits;
  //std::vector<QuasiAlignment> rcFwdHits;

  // there are two possibilities
  bool fwdRc{false};
  bool rcFwd{false};
  // 1. The left is from forward and right is from reverse
  if (leftFwdMap.size() > 0 and rightRcMap.size() > 0) {
    // only consider transcripts that are common between both
    fwdRc = mergeLeftRightMap(rpair, leftFwdMap, rightRcMap, jointHits,//fwdRcHits,
                              editDistance, rmi,maxInsertSize_);
    foundHit = true;
  }
  size_t fwdRcOffset = jointHits.size();

  if (leftRcMap.size() > 0 and rightFwdMap.size() > 0) {
    rcFwd = mergeLeftRightMap(rpair, leftRcMap, rightFwdMap, jointHits,//rcFwdHits,
                              editDistance, rmi, maxInsertSize_);
    foundHit = true;
  }

  // merge two sets if jointHits in different orientations instead of sort
  // use std::merge for this
  std::inplace_merge(jointHits.begin(), jointHits.begin() + fwdRcOffset,
                     jointHits.end(),
                     [](const QuasiAlignment& l, const QuasiAlignment& r)->bool {
                       return l.tid < r.tid;
                     });
  /*
  uint32_t i = 0;
  uint32_t j = 0;
  while(i<fwdRcHits.size() and j<rcFwdHits.size()){
    if(fwdRcHits[i].tid < rcFwdHits[j].tid){
      jointHits.push_back(fwdRcHits[i]);
      i++;
    } else {
      jointHits.push_back(rcFwdHits[j]);
      j++;
    }
  }
  while (i<fwdRcHits.size()) {
    jointHits.push_back(fwdRcHits[i]);
    i++;
  }
  while (j<rcFwdHits.size()) {
    jointHits.push_back(rcFwdHits[j]);
    j++;
  }
  */

  /*if(fwdRc and rcFwd and jointHits.size()>1){
    std::sort(jointHits.begin(), jointHits.end(),
              [](const QuasiAlignment& a, const QuasiAlignment& b)-> bool {
                return a.tid < b.tid;
              } );
  }*/
  return foundHit;
}


template <typename RapMapIndexT>
bool mergeSAInts(
    fastx_parser::ReadSeq& rp, bool lhp,
    std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& leftFwdSAInts,
    std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& leftRcSAInts,
    std::vector<QuasiAlignment>& jointHits, RapMapIndexT& rmi, bool maxNumHits,
    bool consistentHits, rapmap::utils::HitCounters& hctr) {

  //if(leftFwdSAInts.size()==0 and leftRcSAInts.size()==0)
  //	std::cout<<"err\n";

  using OffsetT = typename RapMapIndexT::IndexType;
  AlignerEngine ae_;

  auto& SA = rmi.SA;
  auto& txpStarts = rmi.txpOffsets;
  auto& txpIDs = rmi.positionIDs;

  //uint32_t maxInsertSize_{1000};

  bool foundHit{false};

  SAHitMap leftFwdMap;
  SAHitMap leftRcMap;

  if (leftFwdSAInts.size() > 0 or leftRcSAInts.size() > 0) {
    if (leftFwdSAInts.size() > 1) {
      leftFwdMap = rapmap::hit_manager::unionSAHits(
          leftFwdSAInts, rmi, rp.seq.length(), consistentHits);
    } else if (leftFwdSAInts.size() > 0) {
      auto inHit = leftFwdSAInts.front();
      for (auto i = inHit.begin; i < inHit.end; ++i) {
        auto globalPos = SA[i];
        auto tid = rmi.transcriptAtPosition(globalPos);
        auto txpPos = globalPos - txpStarts[tid];
        leftFwdMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos,
                                           inHit.lcpLength, inHit.queryRC);
        leftFwdMap[tid].active = true;
        leftFwdMap[tid].tqvec.back().matchedLen = inHit.len;
      }
    }

    if (leftRcSAInts.size() > 1) {
      leftRcMap = rapmap::hit_manager::unionSAHits(
          leftRcSAInts, rmi, rp.seq.length(), consistentHits);
    } else if (leftRcSAInts.size() > 0) {
      auto inHit = leftRcSAInts.front();
      for (auto i = inHit.begin; i < inHit.end; ++i) {
        auto globalPos = SA[i];
        auto tid = rmi.transcriptAtPosition(globalPos);
        auto txpPos = globalPos - txpStarts[tid];
        leftRcMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos,
                                          inHit.lcpLength, inHit.queryRC);
        leftRcMap[tid].active = true;
        leftRcMap[tid].tqvec.back().matchedLen = inHit.len;
      }
    }
  }

  bool fwdRc{false};
  bool rcFwd{false};

  if (leftFwdMap.size() > 0) {
    // only consider transcripts that are common between both
    fwdRc = mergeMap(rp, leftFwdMap, jointHits);
    foundHit = true;
  }

  size_t fwdRcOffset = jointHits.size();

  if (leftRcMap.size() > 0) {
    rcFwd = mergeMap(rp, leftRcMap, jointHits);
    foundHit = true;
  }

  // merge two sets if jointHits in different orientations instead of sort
  // use std::merge for this
  std::inplace_merge(jointHits.begin(), jointHits.begin() + fwdRcOffset,
                     jointHits.end(),
                     [](const QuasiAlignment& l, const QuasiAlignment& r)->bool {
                       return l.tid < r.tid;
                     });

  return foundHit;

}



using SAIndex32BitDense =
    RapMapSAIndex<int32_t, RegHashT<uint64_t, rapmap::utils::kmerVal<int32_t>,
                                    rapmap::utils::KmerKeyHasher>>;
using SAIndex64BitDense =
    RapMapSAIndex<int64_t, RegHashT<uint64_t, rapmap::utils::kmerVal<int64_t>,
                                    rapmap::utils::KmerKeyHasher>>;
using SAIndex32BitPerfect =
    RapMapSAIndex<int32_t,
                  PerfectHashT<uint64_t, rapmap::utils::kmerVal<int32_t>>>;
using SAIndex64BitPerfect =
    RapMapSAIndex<int64_t,
                  PerfectHashT<uint64_t, rapmap::utils::kmerVal<int64_t>>>;

template bool mergeLeftRightSAInts<SAIndex32BitDense>(
    fastx_parser::ReadPair& rpair, bool lhp, bool rhp,
    std::vector<SAIntervalHit<int32_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int32_t>>& leftRcSAInts,
    std::vector<SAIntervalHit<int32_t>>& rightFwdSAInts,
    std::vector<SAIntervalHit<int32_t>>& rightRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex32BitDense& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr,
    uint32_t editDistance, uint32_t maxInsertSize_);

template bool mergeLeftRightSAInts<SAIndex64BitDense>(
    fastx_parser::ReadPair& rpair, bool lhp, bool rhp,
    std::vector<SAIntervalHit<int64_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int64_t>>& leftRcSAInts,
    std::vector<SAIntervalHit<int64_t>>& rightFwdSAInts,
    std::vector<SAIntervalHit<int64_t>>& rightRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex64BitDense& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr,
    uint32_t editDistance, uint32_t maxInsertSize_);

template bool mergeLeftRightSAInts<SAIndex32BitPerfect>(
    fastx_parser::ReadPair& rpair, bool lhp, bool rhp,
    std::vector<SAIntervalHit<int32_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int32_t>>& leftRcSAInts,
    std::vector<SAIntervalHit<int32_t>>& rightFwdSAInts,
    std::vector<SAIntervalHit<int32_t>>& rightRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex32BitPerfect& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr,
    uint32_t editDistance, uint32_t maxInsertSize_);

template bool mergeLeftRightSAInts<SAIndex64BitPerfect>(
    fastx_parser::ReadPair& rpair, bool lhp, bool rhp,
    std::vector<SAIntervalHit<int64_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int64_t>>& leftRcSAInts,
    std::vector<SAIntervalHit<int64_t>>& rightFwdSAInts,
    std::vector<SAIntervalHit<int64_t>>& rightRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex64BitPerfect& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr,
    uint32_t editDistance, uint32_t maxInsertSize_);



template bool mergeSAInts<SAIndex32BitDense>(
    fastx_parser::ReadSeq& rp, bool lhp,
    std::vector<SAIntervalHit<int32_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int32_t>>& leftRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex32BitDense& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr);

template bool mergeSAInts<SAIndex64BitDense>(
    fastx_parser::ReadSeq& rp, bool lhp,
    std::vector<SAIntervalHit<int64_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int64_t>>& leftRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex64BitDense& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr);

template bool mergeSAInts<SAIndex32BitPerfect>(
    fastx_parser::ReadSeq& rp, bool lhp,
    std::vector<SAIntervalHit<int32_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int32_t>>& leftRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex32BitPerfect& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr);

template bool mergeSAInts<SAIndex64BitPerfect>(
    fastx_parser::ReadSeq& rp, bool lhp,
    std::vector<SAIntervalHit<int64_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int64_t>>& leftRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex64BitPerfect& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr);







template bool mergeLeftRightMap<SAIndex32BitDense>(
    fastx_parser::ReadPair& rpair, SAHitMap& leftMap, SAHitMap& rightMap,
    std::vector<QuasiAlignment>& jointHits, uint32_t editDistance,
    SAIndex32BitDense& rmi, uint32_t maxInsertSize_);

template bool mergeLeftRightMap<SAIndex64BitDense>(
    fastx_parser::ReadPair& rpair, SAHitMap& leftMap, SAHitMap& rightMap,
    std::vector<QuasiAlignment>& jointHits, uint32_t editDistance,
    SAIndex64BitDense& rmi, uint32_t maxInsertSize_);

template bool mergeLeftRightMap<SAIndex32BitPerfect>(
    fastx_parser::ReadPair& rpair, SAHitMap& leftMap, SAHitMap& rightMap,
    std::vector<QuasiAlignment>& jointHits, uint32_t editDistance,
    SAIndex32BitPerfect& rmi, uint32_t maxInsertSize_);

template bool mergeLeftRightMap<SAIndex64BitPerfect>(
    fastx_parser::ReadPair& rpair, SAHitMap& leftMap, SAHitMap& rightMap,
    std::vector<QuasiAlignment>& jointHits, uint32_t editDistance,
    SAIndex64BitPerfect& rmi, uint32_t maxInsertSize_);
}
}
