#include "HitManager.hpp"
#include "HitManagerSA.hpp"
#include "BooMap.hpp"
#include "FrugalBooMap.hpp"
#include "FastxParser.hpp"

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

    bool mergeLeftRightMap(fastx_parser::ReadPair& rpair,
    		              SAHitMap& leftMap,
						  SAHitMap& rightMap,
						  std::vector<QuasiAlignment>& jointHits
    		              ){
    	//using OffsetT = typename RapMapIndexT::IndexType;
        typedef struct covInfo{
		  public:
              covInfo() {cov=0;endPos=0;rc=false;};
              //covInfo(SATxpQueryPos s,int c, int e,bool r) : saTxpQp(s), cov(c), endPos(e) {rc = r;};
              covInfo(int c, int e,bool r) : cov(c), endPos(e), rc(r) {};
			//SATxpQueryPos saTxpQp ;
            void operator = (const covInfo &CV ) {
                cov = CV.cov;
                endPos = CV.endPos;
                rc = CV.rc;
            }

			int cov ;
			int endPos ;
            bool rc ;
		};

    	constexpr const int32_t signedZero{0};
    	uint32_t maxInsertSize_{1000} ;
    	bool foundHit{false} ;


    		std::set<int> commonTids ;
    		for(auto& phLeft : leftMap){
    			auto leftTid = phLeft.first ;
    			for(auto& phRight: rightMap){
    				auto rightTid = phRight.first ;
    				if(leftTid == rightTid){
    					commonTids.insert(leftTid);
    				}
    			}
    		}
    		if(commonTids.size() > 0){
    			//for each common tid check
    			//if both the end tell us
    			//something
    			for(auto tid:commonTids){
    				if(leftMap[tid].active and rightMap[tid].active){
						auto& lefttqvec = leftMap[tid].tqvec ;
						auto& righttqvec = rightMap[tid].tqvec ;

						std::sort(lefttqvec.begin(),
								  lefttqvec.end(),
								  [](const SATxpQueryPos& a, const SATxpQueryPos& b)-> bool {
															return a.pos < b.pos;
													} );

						std::sort(righttqvec.begin(),
								  righttqvec.end(),
								  [](const SATxpQueryPos& a, const SATxpQueryPos& b)-> bool {
															return a.pos < b.pos;
													} );

						//look at the things in pair
						//Make possible pairs
						//Discard pairs which are
						//impossible in the sense
						//that they are far(1000 bp) apart

						std::map<int32_t , covInfo> leftHitCov ;
                        {
							auto it = lefttqvec.begin();
							while(it != lefttqvec.end()){
								int numOfCov = 0;
                                int32_t hitKey = it->pos - it->queryPos ;
								if(leftHitCov.count(hitKey) == 0){
                                    covInfo cInfo(it->matchedLen,it->pos + it->matchedLen,it->queryRC);
                                    leftHitCov[hitKey] = cInfo;
								}
								else {
									auto& curr = leftHitCov[hitKey];
									if(curr.endPos > it->pos) {
										curr.cov +=  it->pos + it->matchedLen - curr.endPos;
										curr.endPos = it->pos + it->matchedLen;
									} else {
										curr.cov += it->matchedLen;
										curr.endPos = it->pos + it->matchedLen;
									}
								}
								++it ;
							}
                        }

                        std::map<int32_t , covInfo> rightHitCov ;
                        {
							auto it = righttqvec.begin();
							while(it != righttqvec.end()){
								int numOfCov = 0;
								int32_t hitKey = it->pos - it->queryPos ;
								if(rightHitCov.count(hitKey) == 0){
									covInfo cInfo(it->matchedLen,it->pos + it->matchedLen,it->queryRC);
									rightHitCov[hitKey] = cInfo;
								}
								else{
									auto& curr = rightHitCov[hitKey];
									if(curr.endPos > it->pos){
										curr.cov += ( it->pos + it->matchedLen - curr.endPos );
										curr.endPos = it->pos + it->matchedLen ;
									}else{
										curr.cov += it->matchedLen;
										curr.endPos = it->pos + it->matchedLen ;
									}
								}
								++it ;
                            }
                        }


						std::vector<std::pair<int32_t, int32_t>> validPairs ;

			            for(auto& leftTQPos : leftHitCov){
							for(auto& rightTQPos: rightHitCov){

								// Check for the distance
								// go with it otherwise
								// go with something else
								int32_t hitPos1 = leftTQPos.first; // pos - leftTQPos.queryPos ;
								int32_t hitPos2 = rightTQPos.first; // pos - rightTQPos.queryPos ;

								if(std::abs(hitPos1 - hitPos2) < maxInsertSize_){
									validPairs.push_back(std::make_pair(hitPos1, hitPos2)) ;
								}
							}
						}//end for

						//Now I have valid pairs so I can directly make QuasiAlignment vectors
						//
						if(validPairs.size() > 0){
							foundHit = true;
						}

						/*
						auto maxCovPair = std::max_element(validPairs.begin(), validPairs.end(),
														[leftHitCov,rightHitCov](const std::pair<SATxpQueryPos, SATxpQueryPos>& a,
																const std::pair<SATxpQueryPos, SATxpQueryPos>& b){
																int score1 = leftHitCov[static_cast<int32_t>(a.first.pos - a.first.queryPos)]
																						+ (rightHitCov[static_cast<int32_t>(a.second.pos - a.second.queryPos)] ;
																int score2 = leftHitCov[static_cast<int32_t>(b.first.pos - b.first.queryPos)]
																						+ (rightHitCov[static_cast<int32_t>(b.second.pos - b.second.queryPos)] ;
																return (score1 < score2);
						                                 });
						*/
						int maxCovScore = 0;
						for(auto& vp : validPairs){

							auto score = leftHitCov[(vp.first)].cov + rightHitCov[(vp.second)].cov ;
                            if(maxCovScore < score){
								maxCovScore = score;
							}
						}

						for(auto& vp: validPairs){
							//I will try to make QuasiAlignment objects
							//real alignments and edit distance are not present
							//yet

							int32_t hitPos1 = vp.first;  // vp.first.pos - vp.first.queryPos ;
							int32_t hitPos2 = vp.second; // vp.second.pos - vp.second.queryPos ;

							if(leftHitCov[hitPos1].cov + rightHitCov[hitPos2].cov == maxCovScore){
								int32_t startRead1 = std::max(hitPos1, signedZero);
								int32_t startRead2 = std::max(hitPos2, signedZero);

								bool read1First{(startRead1 < startRead2)};
								int32_t fragStartPos = read1First ? startRead1 : startRead2;
								int32_t fragEndPos = read1First ?
                                    (startRead1 + rpair.first.seq.length()) : (startRead2 + rpair.second.seq.length());
								uint32_t fragLen = fragEndPos - fragStartPos;
								//if(leftHitCov[vp.first]->saTxpQp.queryRC != (rightHitCov[vp.second]->saTxpQp).queryRC) {
									jointHits.emplace_back(
										tid,
                                        hitPos1,
										!leftHitCov[vp.first].rc,
										rpair.first.seq.length(),
										0,//leftHitCov[vp.first].saTxpQp.lcpLength,
										fragLen,
										true
									);
									auto& qaln = jointHits.back();
                                    //qaln.fwd = !leftHitCov[vp.first].rc;
									qaln.mateLen = rpair.second.seq.length();
									qaln.matePos = hitPos2 ;
                                    qaln.mateIsFwd = !(rightHitCov[vp.second].rc);  //!vp.second.queryRC ;
									qaln.mateStatus = rapmap::utils::MateStatus::PAIRED_END_PAIRED;
                                    if(leftHitCov[vp.first].cov == rpair.first.seq.length()){
                                        qaln.editD = 0;
                                        qaln.toAlign = true;
                                    }
                                    if(rightHitCov[vp.second].cov == rpair.second.seq.length()){
                                        qaln.mateEditD = 0;
                                        qaln.mateToAlign = true;
                                    }
                                //}
							}
						}
    				}//end if
    			}//end for
    		}else{
    			foundHit = false;
    		}

        if(jointHits.size()>0){
        	std::sort(jointHits.begin(), jointHits.end(),
		[](const QuasiAlignment& a, const QuasiAlignment& b)-> bool {
                     return a.tid < b.tid;
		} );
        }
    	return foundHit ;
    }

    //template <typename readT>
    /*bool mergeLeftRightMap(fastx_parser::ReadPair& rpair,
    		              SAHitMap& leftMap,
						  SAHitMap& rightMap,
						  std::vector<QuasiAlignment>& jointHits
    		               ){
    	//using OffsetT = typename RapMapIndexT::IndexType;

    	constexpr const int32_t signedZero{0};
    	uint32_t maxInsertSize_{1000} ;
    	 bool foundHit{false} ;


    		std::set<int> commonTids ;
    		for(auto& phLeft : leftMap){
    			auto leftTid = phLeft.first ;
    			for(auto& phRight: rightMap){
    				auto rightTid = phRight.first ;
    				if(leftTid == rightTid){
    					commonTids.insert(leftTid);
    				}
    			}
    		}
            //std::cout<<commonTids.size()<<"\n";
    		if(commonTids.size() > 0){
    			//for each common tid check
    			//if both the end tell us
    			//something
    			for(auto tid:commonTids){
    				if(leftMap[tid].active and rightMap[tid].active){
						auto& lefttqvec = leftMap[tid].tqvec ;
						auto& righttqvec = rightMap[tid].tqvec ;
						std::sort(lefttqvec.begin(),
								  lefttqvec.end(),
								  [](const SATxpQueryPos& a, const SATxpQueryPos& b)-> bool {
															return a.pos < b.pos;
													} );

						std::sort(righttqvec.begin(),
								  righttqvec.end(),
								  [](const SATxpQueryPos& a, const SATxpQueryPos& b)-> bool {
															return a.pos < b.pos;
													} );

						//look at the things in pair
						//Make possible pairs
						//Discard pairs which are
						//impossible in the sense
						//that they are far(1000 bp) apart
						std::map<int32_t , int> leftHitCov ;
						{
							std::map<int32_t , int> hitEndPos ;
							auto it = lefttqvec.begin();
							while(it != lefttqvec.end()){
								int numOfCov = 0;
								int32_t hitKey = it->pos - it->queryPos ;
								if(leftHitCov.count(hitKey) == 0){
									leftHitCov[hitKey] = it->matchedLen;
									hitEndPos[hitKey] = it->pos + it->matchedLen;
								}
								else{
									if(hitEndPos[hitKey] > it->pos){
										leftHitCov[hitKey] += ( it->pos + it->matchedLen - hitEndPos[hitKey] );
										hitEndPos[hitKey] = it->pos + it->matchedLen ;
									}else{
										leftHitCov[hitKey] += it->matchedLen;
										hitEndPos[hitKey] = it->pos + it->matchedLen ;
									}
								}
								++it ;
							}
						}

						std::map<int32_t , int> rightHitCov ;
						{
							std::map<int32_t , int> hitEndPos ;
							auto it = righttqvec.begin();
							while(it != righttqvec.end()){
								int numOfCov = 0;
								int32_t hitKey = it->pos - it->queryPos ;
								if(rightHitCov.count(hitKey) == 0){
									rightHitCov[hitKey] = it->matchedLen;
									hitEndPos[hitKey] = it->pos + it->matchedLen;
								}
								else{
									if(hitEndPos[hitKey] > it->pos){
										rightHitCov[hitKey] += ( it->pos + it->matchedLen - hitEndPos[hitKey] );
										hitEndPos[hitKey] = it->pos + it->matchedLen ;
									}else{
										rightHitCov[hitKey] += it->matchedLen;
										hitEndPos[hitKey] = it->pos + it->matchedLen ;
									}
								}
								++it ;
							}
						}


						//std::vector<std::pair<SATxpQueryPos, SATxpQueryPos>> validPairs ;

                        std::vector<std::pair<int32_t, int32_t>> validPairs ;


                        for(auto& leftTQPos : leftHitCov){
                            for(auto& rightTQPos: rightHitCov){

						//for(auto& leftTQPos : lefttqvec){
							//for(auto& rightTQPos: righttqvec){
								// Check for the distance
								// go with it otherwise
								// go with something else
                                int32_t hitPos1 = leftTQPos.first;    //leftTQPos.pos - leftTQPos.queryPos ;
								int32_t hitPos2 = rightTQPos.first;   //rightTQPos.pos - rightTQPos.queryPos ;

								if(std::abs(hitPos1 - hitPos2) < maxInsertSize_){
									validPairs.push_back(std::make_pair(leftTQPos.first, rightTQPos.first)) ;
								}
							}
						}//end for
                        //std::cout<<"validPairs \t" << validPairs.size() <<"\t" <<lefttqvec.size() <<"\t" <<righttqvec.size() << "\n";
						//Now I have valid pairs so I can directly make QuasiAlignment vectors
						//
						if(validPairs.size() > 0){
							foundHit = true;
						}


						int maxCovScore = 0;
						for(auto& vp : validPairs){
								maxCovScore = score ;
							}
						}

						for(auto& vp: validPairs){
							//I will try to make QuasiAlignment objects
							//real alignments and edit distance are not present
							//yet

							int32_t hitPos1 = vp.first;  //  vp.first.pos - vp.first.queryPos ;
                            int32_t hitPos2 = vp.second; //  vp.second.pos - vp.second.queryPos ;

							if(leftHitCov[hitPos1] + rightHitCov[hitPos2] == maxCovScore){
								int32_t startRead1 = std::max(hitPos1, signedZero);
								int32_t startRead2 = std::max(hitPos2, signedZero);

								bool read1First{(startRead1 < startRead2)};
								int32_t fragStartPos = read1First ? startRead1 : startRead2;
								int32_t fragEndPos = read1First ?
									   (startRead2 + rpair.second.seq.length()) : (startRead1 + rpair.first.seq.length());
								uint32_t fragLen = fragEndPos - fragStartPos;

								jointHits.emplace_back(
										tid,
										hitPos1,
										!vp.first,//.queryRC,
										rpair.first.seq.length(),
										0,//vp.first.lcpLength,
										fragLen,
										true
										);
								auto& qaln = jointHits.back();

								// qaln.mateLen = rpair.second.seq.length();
								// qaln.matePos = hitPos2 ;
								// qaln.mateIsFwd = false; //!vp.second.queryRC ;
							}
						}
    				}//end if
    			}//end for
    		}else{
    			foundHit = false;
    		}
    if(jointHits.size()>0){
        std::sort(jointHits.begin(), jointHits.end(),
								  [](const QuasiAlignment& a, const QuasiAlignment& b)-> bool {
                                  return a.tid < b.tid;
								} );
    }
    	return foundHit ;
    }*/


    template <typename RapMapIndexT>
	bool mergeLeftRightSAInts(
						fastx_parser::ReadPair& rpair,
						bool lhp,
						bool rhp,
						std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& leftFwdSAInts,
						std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& leftRcSAInts,
						std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& rightFwdSAInts,
						std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& rightRcSAInts,
						std::vector<QuasiAlignment>& jointHits,
						RapMapIndexT& rmi,
						bool maxNumHits,
						bool consistentHits,
						rapmap::utils::HitCounters& hctr
						){

    	using OffsetT = typename RapMapIndexT::IndexType;
    	AlignerEngine ae_ ;

    	auto& SA = rmi.SA;
    	auto& txpStarts = rmi.txpOffsets;
    	auto& txpIDs = rmi.positionIDs;

    	OffsetT maxInsertSize_{1000} ;

    	bool foundHit{false} ;




    	SAHitMap leftFwdMap ;
    	SAHitMap leftRcMap ;
    	SAHitMap rightFwdMap ;
    	SAHitMap rightRcMap ;

    	if(leftFwdSAInts.size() > 0 or leftRcSAInts.size() > 0){
    		if(leftFwdSAInts.size() > 1){
    			leftFwdMap = rapmap::hit_manager::unionSAHits(leftFwdSAInts, rmi, rpair.first.seq.length(),consistentHits);
    		}else if(leftFwdSAInts.size() > 0){
    			auto inHit = leftFwdSAInts.front() ;
    			for(auto i = inHit.begin ; i < inHit.end; ++i){
    				auto globalPos = SA[i] ;
    				auto tid = rmi.transcriptAtPosition(globalPos) ;
    				auto txpPos = globalPos - txpStarts[tid] ;
    				leftFwdMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos, inHit.lcpLength, inHit.queryRC);
    				leftFwdMap[tid].active = true;
    				leftFwdMap[tid].tqvec.back().matchedLen = inHit.len ;
    			}
    		}

    		if(leftRcSAInts.size() > 1){
    			leftRcMap = rapmap::hit_manager::unionSAHits(leftRcSAInts, rmi, rpair.first.seq.length(), consistentHits);
    		}else if(leftRcSAInts.size() > 0){
    			auto inHit = leftRcSAInts.front() ;
    			for(auto i = inHit.begin ; i < inHit.end; ++i){
    				auto globalPos = SA[i] ;
    				auto tid = rmi.transcriptAtPosition(globalPos) ;
    				auto txpPos = globalPos - txpStarts[tid] ;
    				leftRcMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos, inHit.lcpLength, inHit.queryRC);
    				leftRcMap[tid].active = true;
    				leftRcMap[tid].tqvec.back().matchedLen = inHit.len ;
    			}
    		}
    	}

    	if(rightFwdSAInts.size() > 0 or rightRcSAInts.size() > 0){
    		if(rightFwdSAInts.size() > 1){
    			rightFwdMap = rapmap::hit_manager::unionSAHits(rightFwdSAInts, rmi, rpair.second.seq.length(),consistentHits);
    		}else if(rightFwdSAInts.size() > 0){
    			auto inHit = rightFwdSAInts.front() ;
    			for(auto i = inHit.begin ; i < inHit.end; ++i){
    				auto globalPos = SA[i] ;
    				auto tid = rmi.transcriptAtPosition(globalPos) ;
    				auto txpPos = globalPos - txpStarts[tid] ;
    				rightFwdMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos, inHit.lcpLength, inHit.queryRC);
    				rightFwdMap[tid].active = true;
    				rightFwdMap[tid].tqvec.back().matchedLen = inHit.len ;
    			}
    		}
    		if(rightRcSAInts.size() > 1){
    			rightRcMap = rapmap::hit_manager::unionSAHits(rightRcSAInts, rmi, rpair.second.seq.length(), consistentHits);
    		}else if(rightRcSAInts.size() > 0){
    			auto inHit = rightRcSAInts.front() ;
    			for(auto i = inHit.begin ; i < inHit.end; ++i){
    				auto globalPos = SA[i] ;
    				auto tid = rmi.transcriptAtPosition(globalPos) ;
    				auto txpPos = globalPos - txpStarts[tid] ;
    				rightRcMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos, inHit.lcpLength, inHit.queryRC);
    				rightRcMap[tid].active = true;
    				rightRcMap[tid].tqvec.back().matchedLen = inHit.len ;
    			}
    		}
    	}

    	//there are two possibilities
    	bool fwdRc{false} ;
    	bool rcFwd{false} ;
    	//1. The left is from forward and right is from reverse
    	if(leftFwdMap.size() > 0 and rightRcMap.size() > 0){
    		// only consider transcripts that are common between both
    		fwdRc = mergeLeftRightMap(rpair, leftFwdMap, rightRcMap, jointHits) ;
    		foundHit = true ;
    	}

    	if(leftRcMap.size() > 0 and rightFwdMap.size() > 0){
            rcFwd = mergeLeftRightMap(rpair, leftRcMap, rightFwdMap, jointHits) ;
    		foundHit = true ;
    	}


    	return foundHit ;
    }

    using SAIndex32BitDense = RapMapSAIndex<int32_t, RegHashT<uint64_t, rapmap::utils::kmerVal<int32_t>,
                   					     rapmap::utils::KmerKeyHasher>>;
    using SAIndex64BitDense = RapMapSAIndex<int64_t, RegHashT<uint64_t, rapmap::utils::kmerVal<int64_t>,
   									     rapmap::utils::KmerKeyHasher>>;
    using SAIndex32BitPerfect = RapMapSAIndex<int32_t, PerfectHashT<uint64_t, rapmap::utils::kmerVal<int32_t>>>;
    using SAIndex64BitPerfect = RapMapSAIndex<int64_t, PerfectHashT<uint64_t, rapmap::utils::kmerVal<int64_t>>>;

    template
    bool mergeLeftRightSAInts<SAIndex32BitDense>(
         						fastx_parser::ReadPair& rpair,
         						bool lhp,
         						bool rhp,
         						std::vector<SAIntervalHit<int32_t>>& leftFwdSAInts,
         						std::vector<SAIntervalHit<int32_t>>& leftRcSAInts,
         						std::vector<SAIntervalHit<int32_t>>& rightFwdSAInts,
         						std::vector<SAIntervalHit<int32_t>>& rightRcSAInts,
         						std::vector<QuasiAlignment>& jointHits,
         						SAIndex32BitDense& rmi,
         						bool maxNumHits,
         						bool consistentHits,
         						rapmap::utils::HitCounters& hctr
         						);
    template
    bool mergeLeftRightSAInts<SAIndex64BitDense>(
    							fastx_parser::ReadPair& rpair,
         						bool lhp,
         						bool rhp,
         						std::vector<SAIntervalHit<int64_t>>& leftFwdSAInts,
         						std::vector<SAIntervalHit<int64_t>>& leftRcSAInts,
         						std::vector<SAIntervalHit<int64_t>>& rightFwdSAInts,
         						std::vector<SAIntervalHit<int64_t>>& rightRcSAInts,
         						std::vector<QuasiAlignment>& jointHits,
         						SAIndex64BitDense& rmi,
         						bool maxNumHits,
         						bool consistentHits,
         						rapmap::utils::HitCounters& hctr
         						);
    template
    bool mergeLeftRightSAInts<SAIndex32BitPerfect>(
    							fastx_parser::ReadPair& rpair,
         						bool lhp,
         						bool rhp,
         						std::vector<SAIntervalHit<int32_t>>& leftFwdSAInts,
         						std::vector<SAIntervalHit<int32_t>>& leftRcSAInts,
         						std::vector<SAIntervalHit<int32_t>>& rightFwdSAInts,
         						std::vector<SAIntervalHit<int32_t>>& rightRcSAInts,
         						std::vector<QuasiAlignment>& jointHits,
         						SAIndex32BitPerfect& rmi,
         						bool maxNumHits,
         						bool consistentHits,
         						rapmap::utils::HitCounters& hctr
         						);
    template
    bool mergeLeftRightSAInts<SAIndex64BitPerfect>(
    							fastx_parser::ReadPair& rpair,
         						bool lhp,
         						bool rhp,
         						std::vector<SAIntervalHit<int64_t>>& leftFwdSAInts,
         						std::vector<SAIntervalHit<int64_t>>& leftRcSAInts,
         						std::vector<SAIntervalHit<int64_t>>& rightFwdSAInts,
         						std::vector<SAIntervalHit<int64_t>>& rightRcSAInts,
         						std::vector<QuasiAlignment>& jointHits,
         						SAIndex64BitPerfect& rmi,
         						bool maxNumHits,
         						bool consistentHits,
         						rapmap::utils::HitCounters& hctr
         						);



    }
}
