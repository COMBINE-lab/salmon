#include "Util.hpp"
#include <cstring> 

#define ALLOW_VERBOSE 0

namespace pufferfish {
    namespace util {

void joinReadsAndFilterSingle( pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>>& leftMemClusters,
                               //phmap::flat_hash_map<size_t, std::vector<pufferfish::util::MemCluster>> &leftMemClusters,
                              std::vector<pufferfish::util::JointMems> &jointMemsList,
                              uint32_t perfectCoverage,
                              double coverageRatio) {

    uint32_t maxCoverage{0};
    for (auto &leftClustItr : leftMemClusters) {
        // reference id
        size_t tid = leftClustItr.first;
        // left mem clusters
        auto &lClusts = *(leftClustItr.second);
        // Compare the left clusters to the right clusters to filter by positional constraints
        for (auto lclust = lClusts.begin(); lclust != lClusts.end(); lclust++) {
            auto totalCoverage = lclust->coverage;
            if (totalCoverage >= coverageRatio * maxCoverage or totalCoverage == perfectCoverage) {
                jointMemsList.emplace_back(tid, lclust, lclust, 0, MateStatus::PAIRED_END_LEFT);
                uint32_t currCoverage = jointMemsList.back().coverage();
                if (maxCoverage < currCoverage) {
                    maxCoverage = currCoverage;
                }
            }
        }
    }
}


pufferfish::util::MergeResult joinReadsAndFilter(
                                                 pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>>& leftMemClusters,
                                                 pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>>& rightMemClusters,
                                                 //phmap::flat_hash_map<size_t, std::vector<pufferfish::util::MemCluster>> &leftMemClusters,
                                                 //phmap::flat_hash_map<size_t, std::vector<pufferfish::util::MemCluster>> &rightMemClusters,
                        std::vector<pufferfish::util::JointMems> &jointMemsList,
                        uint32_t maxFragmentLength,
                        uint32_t perfectCoverage,
                        double coverageRatio,
                                                 uint64_t firstDecoyIndex,
                                                 const pufferfish::util::MappingConstraintPolicy& mpol,
                                                 pufferfish::util::HitCounters& hctr) {

  // NOTE : We will fill in `jointMemsList` with iterators to MemClusters from the left and right read.
  // multiple JointMems can share the same iterator (i.e., multiple JointMems can point to the same MemCluster).

    using pufferfish::util::MergeResult;
    MergeResult mergeRes{MergeResult::HAD_NONE};
    
    bool noOrphans = mpol.noOrphans;
    bool noDiscordant = mpol.noDiscordant;
    bool noDovetail = mpol.noDovetail;

    // for filtering concordant chain pairs *on the same target*
    const double thresh = mpol.post_merge_chain_sub_thresh;
    const double ithresh = mpol.inv_post_merge_chain_sub_thresh;
    // for filtering orphan chains with respect to the *best chain for each read*
    const double orphan_chain_sub_thresh = mpol.orphan_chain_sub_thresh;

#if ALLOW_VERBOSE
    std::cerr << "\n[JOINREADSANDFILTER]\n";
#endif // ALLOW_VERBOSE

    // If we will allow orphans in the output, then we compute here the
    // maximum coverage of any optimal chain for the left and right read.
    // Later, we will report orphans, in addition to concordant alignments,
    // if they have sufficiently high coverage with respect to the maximum.
    uint64_t maxLeft{0}, maxRight{0}, maxLeftCnt{0}, maxRightCnt{0};
    bool isMaxLeftAndRight = false;
    if (!noOrphans) {
        for (auto &kv : leftMemClusters) {
#if ALLOW_VERBOSE
          std::cerr << "\ntid:" << kv.first << "\n"; 
#endif // ALLOW_VERBOSE
            auto &lClusts = *(kv.second);
            for (auto clust = lClusts.begin(); clust != lClusts.end(); clust++) {
                if (maxLeft == clust->coverage) {
                    maxLeftCnt += 1;
                } else if (maxLeft < clust->coverage) {
                    maxLeft = clust->coverage;
                    maxLeftCnt = 1;
                }
            }
        } // leftMemClusters
        for (auto &kv : rightMemClusters) {
#if ALLOW_VERBOSE
          std::cerr << "\ntid:" << kv.first << "\n"; 
#endif // ALLOW_VERBOSE
          auto &rClusts = *(kv.second);
          for (auto clust = rClusts.begin(); clust != rClusts.end(); clust++) {
            if (maxRight == clust->coverage) {
              maxRightCnt += 1;
            } else if (maxRight < clust->coverage) {
              maxRight = clust->coverage;
              maxRightCnt = 1;
            }
          }
        } // rightMemClusters
    } // !noOrphans


    // The maximum coverage of a mem cluster for the left or right read
    auto maxLeftOrRight = maxLeft > maxRight ? maxLeft : maxRight;

    //orphan reads should be taken care of maybe with a flag!
    uint32_t maxCoverage{0};
    uint8_t round{0};
    int32_t sameTxpCount{0};
    int32_t numConcordant{0};
    int32_t numDiscordant{0};
    size_t index_for_current_transcript{0};
    bool hadDovetail{false};
    //phmap::parallel_hash_set<uint32_t> refsWithJointMems;
    while (round == 0 or (round == 1 and !jointMemsList.size() and !noDiscordant)) {
      bool concordantSearch = (round == 0);
      for (auto &leftClustItr : leftMemClusters) {
            // reference id
            size_t tid = leftClustItr.first;
            // left mem clusters
            auto &lClusts = *(leftClustItr.second);
            // right mem clusters for the same reference id
            auto &rClusts = rightMemClusters[tid];
            using score_type = decltype(lClusts.begin()->coverage);
            score_type best_pair_in_target{std::numeric_limits<score_type>::min()};
            index_for_current_transcript = jointMemsList.size();
	    // if we are allowing orphans, then don't 
	    // report orphans to any reference that had joint mem hits
	    // regardless of whether they are concordant or discordant.
            /*
	    if (!noOrphans and !rClusts.empty()) { 
		refsWithJointMems.insert(tid);
	    }
            */

            // Compare the left clusters to the right clusters to filter by positional constraints
            for (auto lclust = lClusts.begin(); lclust != lClusts.end(); lclust++) {

                for (auto rclust = rClusts.begin(); rclust != rClusts.end(); rclust++) {
                    // if both the left and right clusters are oriented in the same direction, skip this pair
                    // NOTE: This should be optional as some libraries could allow this.
	            bool satisfiesOri = lclust->isFw != rclust->isFw;
                    if (concordantSearch and !satisfiesOri) { // if priority 0, ends should be concordant
                        continue;
                    }
		    
		    bool isDovetail{false};
		    if (satisfiesOri) {
			    isDovetail = lclust->isFw ? (lclust->approxReadStartPos() > rclust->approxReadStartPos()) :
				    (rclust->approxReadStartPos() > lclust->approxReadStartPos());
			    if (isDovetail and (static_cast<uint64_t>(tid) < firstDecoyIndex)) { hadDovetail = true; }
		    }
		    // if noDovetail is set, then dovetail mappings are considered discordant
                    // otherwise we consider then concordant.
		    if (isDovetail and noDovetail and concordantSearch) {
			continue;
		    }

                    // FILTER 1
                    // filter read pairs based on the fragment length which is approximated by the distance between the left most start and right most hit end
                    int32_t fragmentLen = rclust->lastRefPos() + rclust->lastMemLen() - lclust->firstRefPos();
                    if (lclust->firstRefPos() > rclust->firstRefPos()) {
                        fragmentLen = lclust->lastRefPos() + lclust->lastMemLen() - rclust->firstRefPos();
                    }
                    if (fragmentLen < 0) { // @fatemeh : should we even be checking for this?
                        std::cerr << "Fragment length cannot be smaller than zero!\n";
                        exit(1);
                    }

                    // FILTERING fragments with size smaller than maxFragmentLength
                    // FILTER just in case of priority 0 (round 0)
                    if ((static_cast<uint32_t>(fragmentLen) < maxFragmentLength) or (round > 0)) {
                        // This will add a new potential mapping. Coverage of a mapping for read pairs is left->coverage + right->coverage
                        // If we found a perfect coverage, we would only add those mappings that have the same perfect coverage
                        auto totalCoverage = lclust->coverage + rclust->coverage;
                        if ( (totalCoverage >= coverageRatio * maxCoverage) or
                              (totalCoverage == perfectCoverage) ) {

                            if (totalCoverage >= best_pair_in_target * thresh){
                              ++sameTxpCount;
                              numConcordant += concordantSearch ? 1 : 0;                           
    
                              //if ((totalCoverage > best_pair_in_target * ithresh) and (jointMemsList.size() > index_for_current_transcript)) {
                              //  jointMemsList.erase(jointMemsList.begin()+index_for_current_transcript, jointMemsList.end());
                              //}
                              best_pair_in_target = std::max(best_pair_in_target, totalCoverage);

                              jointMemsList.emplace_back(tid, lclust, rclust, fragmentLen);
                              uint32_t currCoverage = jointMemsList.back().coverage();
                              if (maxCoverage < currCoverage) {
                                  maxCoverage = currCoverage;
                                  if ( (lclust->coverage < maxLeft) or (rclust->coverage < maxRight)) {
                                      isMaxLeftAndRight = false;
                                  } else {
                                      isMaxLeftAndRight = true;
                                  }
                              }
                          } //if (totalCoverage >= best_pair_in_target) 

                        }
                    }
                }
            } // for (auto rclust = rClusts.begin(); rclust != rClusts.end(); rclust++) 

            //auto start_iter = jointMemsList.begin() + index_for_current_transcript;
            //auto nhit = std::distance(start_iter, jointMemsList.end());
            //if (nhit > 1) {
              jointMemsList.erase(
                  std::remove_if(
                    jointMemsList.begin() + index_for_current_transcript,
                    jointMemsList.end(),
                    [best_pair_in_target, thresh](const decltype(*jointMemsList.begin())& jm) -> bool {
                                   return jm.coverage() < thresh * best_pair_in_target;
                    }),
                    jointMemsList.end());
            //}

        } // @fatemeh : this nesting just seems too many levels deep.  Can we re-work the logic here to make things simpler?
        round++;
    }
    numDiscordant = sameTxpCount - numConcordant;
    (void) numDiscordant;

    hctr.numDovetails += hadDovetail ? 1 : 0;
#if ALLOW_VERBOSE
    // If we couldn't find any pair and we are allowed to add orphans
        std::cerr << "isMaxLeftAndRight:" << isMaxLeftAndRight << "\n";
#endif // ALLOW_VERBOSE
   
    // if we've collected any mappings the same transcript, either concordant or discordant (if we are allowing it)
    // then don't consider orphans.
    bool noPairedMappings = (sameTxpCount == 0);
    bool leftOrphan = false; bool rightOrphan = false;
    if (!noOrphans and noPairedMappings and (!jointMemsList.size() or !isMaxLeftAndRight or maxLeftCnt > 1 or maxRightCnt > 1)) {
        auto orphanFiller = [&jointMemsList, &maxCoverage, &coverageRatio, &maxLeftOrRight, &leftOrphan, &rightOrphan, thresh, ithresh, orphan_chain_sub_thresh]
        (pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>> &memClusters,
                 bool isLeft) {
            
            /**
             * NOTE: we handle orphans a little bit differently than we handle (properly paired)
             * chains above.  Since orphans are likely to be from lower-quality fragments and/or 
             * to be more ambiguous than chain pairs, we want to be a bit more aggressive in 
             * how we filter them.  First, rather than use the standard coverageRatio that 
             * we use for paired-end chains, we will instead define the orphan coverage ratio 
             * that comes from the `orphan_chain_sub_thresh`.
             * Second, we will only consider the orphan chains on a target that pass
             * this global threshold, so we don't worry about allowing sub-optimality here.
             **/
            decltype(coverageRatio) orphanCoverageRatio = orphan_chain_sub_thresh;
            constexpr const decltype(thresh) orphanThresh = 1.0;
            constexpr const decltype(thresh) invOrphanThresh = 1.0;

            // fragmentLen is set to 0
            for (auto &clustItr : memClusters) {
              double best_pair_in_target{std::numeric_limits<double>::min()};
              size_t index_for_current_transcript = jointMemsList.size();
              
              // reference id
              size_t tid = clustItr.first;
              // left mem clusters
              auto& Clusts = *(clustItr.second);
              for (auto clust = Clusts.begin(); clust != Clusts.end();
                   clust++) {
                
                if (clust->coverage >= orphanCoverageRatio * maxLeftOrRight) {
                   
                   if (clust->coverage >= best_pair_in_target * orphanThresh) {

                    if ((clust->coverage > best_pair_in_target * invOrphanThresh) and
                        (jointMemsList.size() > index_for_current_transcript)) {
                      jointMemsList.erase(jointMemsList.begin() + index_for_current_transcript,
                                          jointMemsList.end());
                    }
                    
                    best_pair_in_target = std::max(best_pair_in_target, clust->coverage);
                    
                    if (isLeft) {
                      leftOrphan = true;
                      jointMemsList.emplace_back(tid, clust, /*dummy*/ clust, 0,
                                                 MateStatus::PAIRED_END_LEFT);
                    } else {
                      rightOrphan = true;
                      jointMemsList.emplace_back(tid, /*dummy*/ clust, clust, 0,
                                                 MateStatus::PAIRED_END_RIGHT);
                    }
                    uint32_t currCoverage = jointMemsList.back().coverage();
                    if (maxCoverage < currCoverage) {
                      maxCoverage = currCoverage;
                    }

                  } // if (clust->coverage >= best_pair_in_target * orphanThresh) 
                } 

                }

                // don't need this since we always keep only the best orphan per-target
                /*
                jointMemsList.erase(
                    std::remove_if(
                        jointMemsList.begin() + index_for_current_transcript,
                        jointMemsList.end(),
                        [best_pair_in_target, orphanThresh](const decltype(*jointMemsList.begin())& jm) -> bool {
                          return jm.coverage() < orphanThresh * best_pair_in_target;
                        }),
                    jointMemsList.end());
                */
            } // for (auto &clustItr : memClusters) [ over all transcripts]
        };
        orphanFiller(leftMemClusters, true);
        orphanFiller(rightMemClusters, false);
    }
    if (sameTxpCount == 0) {
      if (leftOrphan and !rightOrphan) {
        mergeRes = MergeResult::HAD_ONLY_LEFT;
      } else if (!leftOrphan and rightOrphan) {
        mergeRes = MergeResult::HAD_ONLY_RIGHT;
      } else if (leftOrphan and rightOrphan) {
        mergeRes = MergeResult::HAD_EMPTY_INTERSECTION;
      } else {
        mergeRes = MergeResult::HAD_NONE;
      }
    } else {
      // round is always incremented, so if it's value is 1, we found a concordant
      // mapping and incremented round only one time.
      mergeRes = (round == 1) ? MergeResult::HAD_CONCORDANT : MergeResult::HAD_DISCORDANT;
    }

#if ALLOW_VERBOSE
        std::cerr << "\nBefore filter " << jointMemsList.size() << " maxCov:" << maxCoverage << "\n";
#endif // ALLOW_VERBOSE

    jointMemsList.erase(std::remove_if(jointMemsList.begin(), jointMemsList.end(),
                                       [&maxCoverage, coverageRatio](pufferfish::util::JointMems &pairedReadMems) -> bool {
                                           return pairedReadMems.coverage() < coverageRatio * maxCoverage;
                                       }),
                        jointMemsList.end());

#if ALLOW_VERBOSE
        std::cerr << "\nAfter:" << jointMemsList.size() << " maxCov:" << maxCoverage << "\n";
        std::cerr << "\n[END OF JOINREADSANDFILTER]\n";
#endif // ALLOW_VERBOSE

    return mergeRes;
}



      char * getRefSeqOwned(compact::vector<uint64_t, 2> &refseq, uint64_t refAccPos, uint32_t refLen) {
        if (refLen == 0) return nullptr;
        char* seq = new char[refLen];
        std::memset(seq, 0, refLen);
        uint64_t c = 0;
        uint64_t bucket_offset = (refAccPos) * 2;
        auto len_on_vector = refLen * 2;
        int32_t toFetch = len_on_vector;
        while (toFetch > 0) {
          uint32_t len = (toFetch >= 64) ? 64 : toFetch;
          toFetch -= len;
          uint64_t word = refseq.get_int(bucket_offset, len);
          for (uint32_t i = 0; i < len; i += 2) {
            uint8_t next_bits = ((word >> i) & 0x03);
            seq[c++] = "ACGT"[next_bits];
          }
          bucket_offset += len;
        }
        return seq;
      }


        char complement(char &c) {
          switch (c) {
            case 'A':
              c = 'T';
                  return c;
            case 'T':
              c = 'A';
                  return c;
            case 'C':
              c = 'G';
                  return c;
            case 'G':
              c = 'C';
                  return c;
          }
          return 'N';
        }

        std::string revcomp(std::string s) {
          int n = s.size();
          int halfLength = s.size() / 2;
          for (int i = 0; i < halfLength; i++) {
            char temp = complement(s[i]);
            s[i] = complement(s[n - 1 - i]);
            s[n - 1 - i] = temp;
          }
          if (s.size() % 2 != 0) {
            s[halfLength] = complement(s[halfLength]);
          }
          return s;
        }

        bool isRevcomp(std::string s) {
          int n = s.size();
          int halfLength = n / 2;
          for (int i = 0; i < halfLength; i++) {
            char temp = complement(s[n - 1 - i]);
            if (temp != s[i])
              return false;
          }
          return true;
        }

        std::vector<std::pair<uint64_t, bool>> explode(const stx::string_view str, const char &ch) {
          std::string next;
          std::vector<std::pair<uint64_t, bool>> result;
          // For each character in the string
          for (auto it = str.begin(); it != str.end(); it++) {
            // If we've hit the terminal character
            if (*it == '+' or *it == '-') {
              bool orientation = true;
              // If we have some characters accumulated
              // Add them to the result vector
              if (!next.empty()) {
                if (*it == '-') {
                  orientation = false;
                }
                try {
                  uint64_t nid = std::stoll(next);
                  result.emplace_back(nid, orientation);
                } catch (std::exception &e) {
                  // not a numeric contig id
                  std::cerr << "tried to convert " << next << " into a long long\n";
                  std::exit(1);
                }
                next.clear();
              }
            } else if (*it != ch) {
              // Accumulate the next character into the sequence
              next += *it;
            }
          }
          if (!next.empty()) {
            std::cerr << "impossible is the opposite of possible " << next << "\n";
            std::cerr << "The line is " << str << "\n";
            result.emplace_back(std::stoll(next),
                                true); // this case shouldn't even happen
          }
          return result;
        }


        std::vector<extension> getExts(uint8_t edgeVec) {
          std::vector<extension> ext;
          uint8_t mask = 1;
          std::vector<char> nuclmap = {'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A'};
          for (uint8_t i = 0; i < 8; i++) {
            if (edgeVec & (mask << i)) {
              if (i < 4)
                ext.push_back({nuclmap[i], Direction::FORWARD});
              else
                ext.push_back({nuclmap[i], Direction::BACKWORD});
            }
          }
          return ext;

        }

/*std::vector<std::pair<std::string, bool>> explode(const stx::string_view str,
                                                  const char& ch) {
  std::string next;
  std::vector<std::pair<std::string, bool>> result;
  // For each character in the string
  for (auto it = str.begin(); it != str.end(); it++) {
    // If we've hit the terminal character
    if (*it == '+' or *it == '-') {
      bool orientation = true;
      // If we have some characters accumulated
      // Add them to the result vector
      if (!next.empty()) {
        if (*it == '-') {
          orientation = false;
        }
        result.emplace_back(next, orientation);
        next.clear();
      }
    } else if (*it != ch) {
      // Accumulate the next character into the sequence
      next += *it;
    }
  }
  if (!next.empty())
    result.emplace_back(next, true); // this case shouldn't even happen
  return result;
}
*/
        bool is_number(const std::string &s) {
          return !s.empty() && std::find_if(s.begin(), s.end(), [](char c) {
              return !std::isdigit(c);
          }) == s.end();
        }

// tokenize the file names
// later TODO: replace string streams with string_view 
        std::vector<std::string> tokenize(const std::string &s, char delim) {
          std::stringstream ss(s);
          std::string item;
          std::vector<std::string> elems;
          while (std::getline(ss, item, delim)) {
            elems.push_back(item);
          }
          return elems;
        }


// Avoiding un-necessary stream creation + replacing strings with string view
// is a bit > than a 2x win!
// implementation from : https://marcoarena.wordpress.com/tag/string_view/

        std::vector<stx::string_view> split(stx::string_view str, char delims) {
          std::vector<stx::string_view> ret;

          stx::string_view::size_type start = 0;
          auto pos = str.find_first_of(delims, start);
          while (pos != stx::string_view::npos) {
            if (pos != start) {
              ret.push_back(str.substr(start, pos - start));
            }
            start = pos + 1;
            pos = str.find_first_of(delims, start);
          }
          if (start < str.length()) {
            ret.push_back(str.substr(start, str.length() - start));
          }
          return ret;
        }
    }
}
