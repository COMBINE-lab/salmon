#include "HitManager.hpp"
#include "HitManagerSA.hpp"
#include "BooMap.hpp"
#include "FrugalBooMap.hpp"
#include "FastxParser.hpp"
#include "SelectiveAlignment.hpp"
#include "SalmonUtils.hpp"

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
    bool mergeLeftRightMap(fastx_parser::ReadPair& rpair,
    		              SAHitMap& leftMap,
						  SAHitMap& rightMap,
						  std::vector<QuasiAlignment>& jointHits,
                          uint32_t editDistance,
                          RapMapIndexT& rmi){
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

        SECollector<RapMapIndexT> hitSECollector(&rmi);
        std::vector<QuasiAlignment> hits;

        if(salmon::utils::GroundTruth::foundTxps.count(rpair.first.name)==0){
            std::vector<std::string> txps;
            salmon::utils::GroundTruth::foundTxps[rpair.first.name] = txps;
        }


    		std::set<int> commonTids ;
    		for(auto& phLeft : leftMap){
    			auto leftTid = phLeft.first ;
    			for(auto& phRight: rightMap){
    				auto rightTid = phRight.first ;
                     //if(rpair.second.name=="13915359_0_6053_768_153/2"){
                    //    std::cout<<rpair.second.name << " " << rmi.txpNames[rightTid] << " " << rpair.first.name << " " << rmi.txpNames[leftTid] <<"\n";
                    //}

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

                //if(rpair.second.name=="13915359_0_6053_768_153/2"){
                //    std::cout<<rpair.second.name << "\n " << rmi.txpNames[tid] <<"\n\n";
                //}

   				if(leftMap[tid].active and rightMap[tid].active){

                    //if(rpair.second.name=="13915359_0_6053_768_153/2"){
                    //    std::cout<<rpair.second.name << "\n\n after active " << rmi.txpNames[tid] <<"\n\n";
                   // }


                    //salmon::utils::GroundTruth::foundTxps[rpair.first.name].push_back(rmi.txpNames[tid]);
                    //if(rpair.first.name=="13915359_0_6053_768_153/2")
                     //   std::cout<<"\n\n added to the vector: "<<salmon::utils::GroundTruth::foundTxps[rpair.first.name].front() << "\n\n";

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

                        hits.clear();
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
                                break;
                            }
						}
                        /*hitSECollector(rpair.first,rpair.second, hits, editDistance);

                        auto minDist = 200;
                        std::for_each(hits.begin(), hits.end(),
                            [&minDist](QuasiAlignment& a) {
                            if (a.editD < minDist and a.editD != -1) { minDist = a.editD; }
                        });



                        hits.erase(std::remove_if(hits.begin(),hits.end(),
                            [&minDist](QuasiAlignment& a) {
                                return (a.editD > minDist);
                        }), hits.end());*/

                        //if(hits.size()>0)
                        //    jointHits.push_back(hits.back());

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
						rapmap::utils::HitCounters& hctr,
						uint32_t editDistance){

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
        /*if(rpair.second.seq == "AAAACTAGCCAGGCGTGGTGGTGGGCACCTGTAATCCCAGCTACGTGGGAGGCTGAGGCAAGAGAATTGCTTGAA" or rpair.second.seq=="TTCAAGCAATTCTCTTGCCTCAGCCTCCCACGTAGCTGGGATTACAGGTGCCCACCACCACGCCTGGCTAGTTTT"){
            std::cout<<"\n\n" <<leftRcMap.size() << " " << rightFwdMap.size() << "\n\n" <<leftFwdMap.size() << " " << rightRcMap.size() << "\n\n";

        }*/


    	//there are two possibilities
    	bool fwdRc{false} ;
    	bool rcFwd{false} ;
    	//1. The left is from forward and right is from reverse
    	if(leftFwdMap.size() > 0 and rightRcMap.size() > 0){
    		// only consider transcripts that are common between both
            /*if(rpair.second.seq == "AAAACTAGCCAGGCGTGGTGGTGGGCACCTGTAATCCCAGCTACGTGGGAGGCTGAGGCAAGAGAATTGCTTGAA" or rpair.second.seq=="TTCAAGCAATTCTCTTGCCTCAGCCTCCCACGTAGCTGGGATTACAGGTGCCCACCACCACGCCTGGCTAGTTTT"){
                std::cout<<"\n\n" <<leftFwdMap.size() << " " << rightRcMap.size() << "\n\n";

            }*/
            fwdRc = mergeLeftRightMap(rpair, leftFwdMap, rightRcMap, jointHits,editDistance,rmi) ;
    		foundHit = true ;
    	}

    	if(leftRcMap.size() > 0 and rightFwdMap.size() > 0){
            /*if(rpair.second.seq == "AAAACTAGCCAGGCGTGGTGGTGGGCACCTGTAATCCCAGCTACGTGGGAGGCTGAGGCAAGAGAATTGCTTGAA" or rpair.second.seq=="TTCAAGCAATTCTCTTGCCTCAGCCTCCCACGTAGCTGGGATTACAGGTGCCCACCACCACGCCTGGCTAGTTTT"){
                std::cout<<"\n\n" <<leftRcMap.size() << " " << rightFwdMap.size() << "\n\n";

            }*/

            rcFwd = mergeLeftRightMap(rpair, leftRcMap, rightFwdMap, jointHits,editDistance,rmi) ;
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
         						rapmap::utils::HitCounters& hctr,
         						uint32_t editDistance);
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
         						rapmap::utils::HitCounters& hctr,
                                uint32_t editDistance);
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
         						rapmap::utils::HitCounters& hctr,
         						uint32_t editDistance);
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
         						rapmap::utils::HitCounters& hctr,
         						uint32_t editDistance);


    template
    bool mergeLeftRightMap<SAIndex32BitDense>(
                          fastx_parser::ReadPair& rpair,
    		              SAHitMap& leftMap,
						  SAHitMap& rightMap,
						  std::vector<QuasiAlignment>& jointHits,
                          uint32_t editDistance,
                          SAIndex32BitDense& rmi);

    template
    bool mergeLeftRightMap<SAIndex64BitDense>(
                          fastx_parser::ReadPair& rpair,
    		              SAHitMap& leftMap,
						  SAHitMap& rightMap,
						  std::vector<QuasiAlignment>& jointHits,
                          uint32_t editDistance,
                          SAIndex64BitDense& rmi);

    template
    bool mergeLeftRightMap<SAIndex32BitPerfect>(
                          fastx_parser::ReadPair& rpair,
    		              SAHitMap& leftMap,
						  SAHitMap& rightMap,
						  std::vector<QuasiAlignment>& jointHits,
                          uint32_t editDistance,
                          SAIndex32BitPerfect& rmi);

    template
    bool mergeLeftRightMap<SAIndex64BitPerfect>(
                          fastx_parser::ReadPair& rpair,
    		              SAHitMap& leftMap,
						  SAHitMap& rightMap,
						  std::vector<QuasiAlignment>& jointHits,
                          uint32_t editDistance,
                          SAIndex64BitPerfect& rmi);


    }
}
