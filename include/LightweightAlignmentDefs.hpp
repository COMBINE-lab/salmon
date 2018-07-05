#ifndef LIGHTWEIGHT_ALIGNMENT_DEFS_HPP
#define LIGHTWEIGHT_ALIGNMENT_DEFS_HPP

#include "BWAMemStaticFuncs.hpp"
#include "EffectiveLengthStats.hpp"
#include "RapMapUtils.hpp"


using BulkReadExpT = ReadExperiment<EquivalenceClassBuilder<TGValue>>;

class SMEMAlignment {
public:
  SMEMAlignment()
      : pos(0), fwd(false), mateIsFwd(false),
        transcriptID_(std::numeric_limits<TranscriptID>::max()),
        format_(LibraryFormat::formatFromID(0)), score_(0.0), fragLength_(0),
        logProb(salmon::math::LOG_0), logBias(salmon::math::LOG_0) {}

  SMEMAlignment(TranscriptID transcriptIDIn, LibraryFormat format,
                double scoreIn = 0.0, int32_t hitPosIn = 0,
                uint32_t fragLengthIn = 0,
                double logProbIn = salmon::math::LOG_0)
      : pos(hitPosIn), fwd(false), mateIsFwd(false),
        transcriptID_(transcriptIDIn), format_(format), score_(scoreIn),
        fragLength_(fragLengthIn), fragLen(fragLengthIn), logProb(logProbIn) {}

  SMEMAlignment(const SMEMAlignment& o) = default;
  SMEMAlignment(SMEMAlignment&& o) = default;
  SMEMAlignment& operator=(SMEMAlignment& o) = default;
  SMEMAlignment& operator=(SMEMAlignment&& o) = default;

  inline TranscriptID transcriptID() const { return transcriptID_; }
  inline uint32_t fragLength() const { return fragLength_; }
  inline uint32_t fragLengthPedantic(uint32_t txpLen) const {
    return fragLength_;
  }
  inline LibraryFormat libFormat() const { return format_; }
  inline double score() const { return score_; }
  inline int32_t hitPos() const { return pos; }
  // inline double coverage() {  return static_cast<double>(kmerCount) /
  // fragLength_; };
  uint32_t kmerCount;
  double logProb;
  double logBias;
  template <typename Archive> void save(Archive& archive) const {
    archive(transcriptID_, format_.formatID(), score_, pos, fragLength_);
  }

  template <typename Archive> void load(Archive& archive) {
    uint8_t formatID;
    archive(transcriptID_, formatID, score_, pos, fragLength_);
    format_ = LibraryFormat::formatFromID(formatID);
  }

  rapmap::utils::MateStatus mateStatus;
  int32_t pos;
  int32_t matePos; // JUST FOR COMPATIBILITY WITH QUASI!
  bool fwd;
  bool mateIsFwd;
  uint32_t readLen;
  uint32_t mateLen;
  uint32_t fragLen;

private:
  TranscriptID transcriptID_;
  LibraryFormat format_;
  double score_;
  uint32_t fragLength_;
};

uint32_t basesCovered(std::vector<uint32_t>& kmerHits) {
  std::sort(kmerHits.begin(), kmerHits.end());
  uint32_t covered{0};
  uint32_t lastHit{0};
  uint32_t kl{20};
  for (auto h : kmerHits) {
    covered += std::min(h - lastHit, kl);
    lastHit = h;
  }
  return covered;
}

uint32_t basesCovered(std::vector<uint32_t>& posLeft,
                      std::vector<uint32_t>& posRight) {
  return basesCovered(posLeft) + basesCovered(posRight);
}

class KmerVote {
public:
  KmerVote(int32_t vp, uint32_t rp, uint32_t vl)
      : votePos(vp), readPos(rp), voteLen(vl) {}
  int32_t votePos{0};
  uint32_t readPos{0};
  uint32_t voteLen{0};
  /*
  std::string str(){
      return "<" + votePos  + ", "  + readPos  + ", "  + voteLen + ">";
  }
  */
};
class MatchFragment {
public:
  MatchFragment(uint32_t refStart_, uint32_t queryStart_, uint32_t length_)
      : refStart(refStart_), queryStart(queryStart_), length(length_) {}

  uint32_t refStart, queryStart, length;
  uint32_t weight;
  double score;
};

bool precedes(const MatchFragment& a, const MatchFragment& b) {
  return (a.refStart + a.length) < b.refStart and
         (a.queryStart + a.length) < b.queryStart;
}

class TranscriptHitList {
public:
  int32_t bestHitPos{0};
  uint32_t bestHitCount{0};
  double bestHitScore{0.0};

  std::vector<KmerVote> votes;
  std::vector<KmerVote> rcVotes;

  uint32_t targetID;
  uint32_t fwdCov{0};
  uint32_t revCov{0};

  bool isForward_{true};

  void addFragMatch(uint32_t tpos, uint32_t readPos, uint32_t voteLen) {
    int32_t votePos =
        static_cast<int32_t>(tpos) - static_cast<int32_t>(readPos);
    votes.emplace_back(votePos, readPos, voteLen);
    fwdCov += voteLen;
  }

  void addFragMatchRC(uint32_t tpos, uint32_t readPos, uint32_t voteLen,
                      uint32_t readLen) {
    // int32_t votePos = static_cast<int32_t>(tpos) - (readPos) + voteLen;
    int32_t votePos = static_cast<int32_t>(tpos) - (readLen - readPos);
    rcVotes.emplace_back(votePos, readPos, voteLen);
    revCov += voteLen;
  }

  uint32_t totalNumHits() { return std::max(votes.size(), rcVotes.size()); }

  bool computeBestLocFast_(std::vector<KmerVote>& sVotes,
                           Transcript& transcript, std::string& read, bool isRC,
                           int32_t& maxClusterPos, uint32_t& maxClusterCount,
                           double& maxClusterScore) {
    bool updatedMaxScore{true};
    if (sVotes.size() == 0) {
      return updatedMaxScore;
    }
    uint32_t readLen = read.length();
    uint32_t votePos = sVotes.front().votePos;

    uint32_t cov = isRC ? revCov : fwdCov;
    if (cov > maxClusterCount) {
      maxClusterCount = cov;
      maxClusterPos = votePos;
      maxClusterScore = maxClusterCount / static_cast<double>(readLen);
      updatedMaxScore = true;
    }
    return updatedMaxScore;
  }

  bool computeBestLoc_(std::vector<KmerVote>& sVotes, Transcript& transcript,
                       std::string& read, bool isRC, int32_t& maxClusterPos,
                       uint32_t& maxClusterCount, double& maxClusterScore) {
    // Did we update the highest-scoring cluster? This will be set to
    // true iff we have a cluster of a higher score than the score
    // currently given in maxClusterCount.
    bool updatedMaxScore{false};

    if (sVotes.size() == 0) {
      return updatedMaxScore;
    }

    struct VoteInfo {
      uint32_t coverage = 0;
      int32_t rightmostBase = 0;
    };

    uint32_t readLen = read.length();

    boost::container::flat_map<uint32_t, VoteInfo> hitMap;
    int32_t currClust{static_cast<int32_t>(sVotes.front().votePos)};

    for (size_t j = 0; j < sVotes.size(); ++j) {

      int32_t votePos = sVotes[j].votePos;
      uint32_t readPos = sVotes[j].readPos;
      uint32_t voteLen = sVotes[j].voteLen;

      if (votePos >= currClust) {
        if (votePos - currClust > 10) {
          currClust = votePos;
        }
        auto& hmEntry = hitMap[currClust];

        hmEntry.coverage += std::min(voteLen, (votePos + readPos + voteLen) -
                                                  hmEntry.rightmostBase);
        hmEntry.rightmostBase = votePos + readPos + voteLen;
      } else if (votePos < currClust) {
        std::cerr << "Should not have votePos = " << votePos
                  << " <  currClust = " << currClust << "\n";
        std::exit(1);
      }

      if (hitMap[currClust].coverage > maxClusterCount) {
        maxClusterCount = hitMap[currClust].coverage;
        maxClusterPos = currClust;
        maxClusterScore = maxClusterCount / static_cast<double>(readLen);
        updatedMaxScore = true;
      }
    }
    return updatedMaxScore;
  }

  bool computeBestLoc2_(std::vector<KmerVote>& sVotes, uint32_t tlen,
                        int32_t& maxClusterPos, uint32_t& maxClusterCount,
                        double& maxClusterScore) {

    bool updatedMaxScore{false};

    if (sVotes.size() == 0) {
      return updatedMaxScore;
    }

    double weights[] = {1.0,
                        0.983471453822,
                        0.935506985032,
                        0.860707976425,
                        0.765928338365,
                        0.6592406302,
                        0.548811636094,
                        0.441902209585,
                        0.344153786865,
                        0.259240260646,
                        0.188875602838};

    int32_t maxGap = 4;
    uint32_t leftmost = (sVotes.front().votePos > maxGap)
                            ? (sVotes.front().votePos - maxGap)
                            : 0;
    uint32_t rightmost = static_cast<uint32_t>(
                             std::min(sVotes.back().votePos + maxGap,
                                      static_cast<int32_t>(tlen)));

    uint32_t span = (rightmost - leftmost);
    std::vector<double> probAln(span, 0.0);
    double kwidth = 1.0 / (2.0 * maxGap);

    size_t nvotes = sVotes.size();
    for (size_t j = 0; j < nvotes; ++j) {
      uint32_t votePos = sVotes[j].votePos;
      uint32_t voteLen = sVotes[j].voteLen;

      auto x = j + 1;
      while (x < nvotes and static_cast<uint32_t>(sVotes[x].votePos) == votePos) {
        voteLen += sVotes[x].voteLen;
        j += 1;
        x += 1;
      }

      uint32_t dist{0};
      size_t start = (votePos >= static_cast<uint32_t>(maxGap)) ? (votePos - maxGap - leftmost)
                                         : (votePos - leftmost);
      size_t mid = votePos - leftmost;
      size_t end = std::min(votePos + maxGap - leftmost, rightmost - leftmost);
      for (size_t k = start; k < end; k += 1) {
        dist = (mid > k) ? mid - k : k - mid;
        probAln[k] += weights[dist] * voteLen;
        if (probAln[k] > maxClusterScore) {
          maxClusterScore = probAln[k];
          maxClusterPos = k + leftmost;
          updatedMaxScore = true;
        }
      }
    }

    return updatedMaxScore;
  }

  inline uint32_t numSampledHits_(Transcript& transcript, std::string& readIn,
                                  int32_t votePos, int32_t posInRead,
                                  int32_t voteLen, bool isRC,
                                  uint32_t numTries) {

    // The read starts at this position in the transcript (may be negative!)
    int32_t readStart = votePos;
    // The (uncorrected) length of the read
    int32_t readLen = readIn.length();
    // Pointer to the sequence of the read
    const char* read = readIn.c_str();
    // Don't mess around with unsigned arithmetic here
    int32_t tlen = transcript.RefLength;

    // If the read starts before the first base of the transcript,
    // trim off the initial overhang  and correct the other variables
    if (readStart < 0) {
      if (isRC) {
        uint32_t correction = -readStart;
        // std::cerr << "readLen = " << readLen << ", posInRead = " << posInRead
        // << ", voteLen = " << voteLen << ", correction = " << correction <<
        // "\n";  std::cerr << "tlen = " << tlen << ", votePos = " << votePos <<
        // "\n";
        read += correction;
        readLen -= correction;
        posInRead -= correction;
        readStart = 0;
      } else {
        uint32_t correction = -readStart;
        read += correction;
        readLen -= correction;
        posInRead -= correction;
        readStart = 0;
      }
    }
    // If the read hangs off the end of the transcript,
    // shorten its effective length.
    if (readStart + readLen >= tlen) {
      if (isRC) {
        uint32_t correction = (readStart + readLen) - transcript.RefLength + 1;
        // std::cerr << "Trimming RC hit: correction = " << correction << "\n";
        // std::cerr << "untrimmed read : "  << read << "\n";
        read += correction;
        readLen -= correction;
        if (voteLen > readLen) {
          voteLen = readLen;
        }
        posInRead = 0;
      } else {
        readLen = tlen - (readStart + 1);
        voteLen = std::max(voteLen, readLen - (posInRead + voteLen));
      }
    }
    // Finally, clip any reverse complement reads starting at 0
    if (isRC) {

      if (voteLen > readStart) {
        readLen -= (readLen - (posInRead + voteLen));
      }
    }

    // If the read is too short, it's not useful
    if (readLen <= 15) {
      return 0;
    }
    // The step between sample centers (given the number of samples we're going
    // to take)
    double step = (readLen - 1) / static_cast<double>(numTries - 1);
    // The strand of the transcript from which we'll extract sequence
    auto dir = (isRC) ? salmon::stringtools::strand::reverse
                      : salmon::stringtools::strand::forward;

    bool superVerbose{false};

    if (superVerbose) {
      std::stringstream ss;
      ss << "Supposed hit " << (isRC ? "RC" : "") << "\n";
      ss << "info: votePos = " << votePos << ", posInRead = " << posInRead
         << ", voteLen = " << voteLen << ", readLen = " << readLen
         << ", tran len = " << tlen << ", step = " << step << "\n";
      if (readStart + readLen > tlen) {
        ss << "ERROR!!!\n";
        std::cerr << "[[" << ss.str() << "]]";
        std::exit(1);
      }
      ss << "Transcript name = " << transcript.RefName << "\n";
      ss << "T : ";
      try {
        for (int32_t j = 0; j < readLen; ++j) {
          if (isRC) {
            if (j == posInRead) {
              char red[] = "\x1b[30m";
              red[3] = '0' + static_cast<char>(fmt::RED);
              ss << red;
            }

            if (j == posInRead + voteLen) {
              const char RESET_COLOR[] = "\x1b[0m";
              ss << RESET_COLOR;
            }
            ss << transcript.charBaseAt(readStart + readLen - j, dir);
          } else {
            if (j == posInRead) {
              char red[] = "\x1b[30m";
              red[3] = '0' + static_cast<char>(fmt::RED);
              ss << red;
            }

            if (j == posInRead + voteLen) {
              const char RESET_COLOR[] = "\x1b[0m";
              ss << RESET_COLOR;
            }

            ss << transcript.charBaseAt(readStart + j);
          }
        }
        ss << "\n";
        char red[] = "\x1b[30m";
        red[3] = '0' + static_cast<char>(fmt::RED);
        const char RESET_COLOR[] = "\x1b[0m";

        ss << "R : " << std::string(read, posInRead) << red
           << std::string(read + posInRead, voteLen) << RESET_COLOR;
        if (readLen > posInRead + voteLen) {
          ss << std::string(read + posInRead + voteLen);
        }
        ss << "\n\n";
      } catch (std::exception& e) {
        std::cerr << "EXCEPTION !!!!!! " << e.what() << "\n";
      }
      std::cerr << ss.str() << "\n";
      ss.clear();
    }

    // The index of the current sample within the read
    int32_t readIndex = 0;

    // The number of loci in the subvotes and their
    // offset patternns
    size_t lpos = 3;
    int leftPattern[] = {-4, -2, 0};
    int rightPattern[] = {0, 2, 4};
    int centerPattern[] = {-4, 0, 4};

    // The number of subvote hits we've had
    uint32_t numHits = 0;
    // Take the samples
    for (size_t i = 0; i < numTries; ++i) {
      // The sample will be centered around this point
      readIndex =
          static_cast<uint32_t>(std::round(readStart + i * step)) - readStart;

      // The number of successful sub-ovtes we have
      uint32_t subHit = 0;
      // Select the center sub-vote pattern, unless we're near the end of a read
      int* pattern = &centerPattern[0];
      if (readIndex + pattern[0] < 0) {
        pattern = &rightPattern[0];
      } else if (readIndex + pattern[lpos - 1] >= readLen) {
        pattern = &leftPattern[0];
      }

      // collect the subvotes
      for (size_t j = 0; j < lpos; ++j) {
        // the pattern offset
        int offset = pattern[j];
        // and sample position it implies within the read
        int readPos = readIndex + offset;

        if (readStart + readPos >= tlen) {
          std::cerr << "offset = " << offset << ", readPos = " << readPos
                    << ", readStart = " << readStart
                    << ", readStart + readPos = " << readStart + readPos
                    << ", tlen = " << transcript.RefLength << "\n";
        }

        subHit +=
            (isRC)
                ? (transcript.charBaseAt(readStart + readLen - readPos, dir) ==
                   salmon::stringtools::charCanon[static_cast<uint8_t>(read[readPos])])
                : (transcript.charBaseAt(readStart + readPos) ==
                   salmon::stringtools::charCanon[static_cast<uint8_t>(read[readPos])]);
      }
      // if the entire subvote was successful, this is a hit
      numHits += (subHit == lpos);
    }
    // return the number of hits we had
    return numHits;
  }

  bool computeBestLoc3_(std::vector<KmerVote>& sVotes, Transcript& transcript,
                        std::string& read, bool isRC, int32_t& maxClusterPos,
                        uint32_t& maxClusterCount, double& maxClusterScore) {

    bool updatedMaxScore{false};

    if (sVotes.size() == 0) {
      return updatedMaxScore;
    }

    struct LocHitCount {
      int32_t loc;
      uint32_t nhits;
    };

    uint32_t numSamp = 15;
    std::vector<LocHitCount> hitCounts;
    size_t nvotes = sVotes.size();
    int32_t prevPos = -std::numeric_limits<int32_t>::max();
    for (size_t j = 0; j < nvotes; ++j) {
      int32_t votePos = sVotes[j].votePos;
      int32_t posInRead = sVotes[j].readPos;
      int32_t voteLen = sVotes[j].voteLen;
      if (prevPos == votePos) {
        continue;
      }
      auto numHits = numSampledHits_(transcript, read, votePos, posInRead,
                                     voteLen, isRC, numSamp);
      hitCounts.push_back({votePos, numHits});
      prevPos = votePos;
    }

    uint32_t maxGap = 8;
    uint32_t hitIdx = 0;
    uint32_t accumHits = 0;
    int32_t hitLoc = hitCounts[hitIdx].loc;
    while (hitIdx < hitCounts.size()) {
      uint32_t idx2 = hitIdx;
      while (idx2 < hitCounts.size() and
             std::abs(hitCounts[idx2].loc - hitLoc) <= maxGap) {
        accumHits += hitCounts[idx2].nhits;
        ++idx2;
      }

      double score = static_cast<double>(accumHits) / numSamp;
      if (score > maxClusterScore) {
        maxClusterCount = accumHits;
        maxClusterScore = score;
        maxClusterPos = hitCounts[hitIdx].loc;
        updatedMaxScore = true;
      }
      accumHits = 0;
      ++hitIdx;
      hitLoc = hitCounts[hitIdx].loc;
    }

    return updatedMaxScore;
  }

  bool computeBestChain(Transcript& transcript, std::string& read) {
    std::sort(votes.begin(), votes.end(),
              [](const KmerVote& v1, const KmerVote& v2) -> bool {
                if (v1.votePos == v2.votePos) {
                  return v1.readPos < v2.readPos;
                }
                return v1.votePos < v2.votePos;
              });

    std::sort(rcVotes.begin(), rcVotes.end(),
              [](const KmerVote& v1, const KmerVote& v2) -> bool {
                if (v1.votePos == v2.votePos) {
                  return v1.readPos < v2.readPos;
                }
                return v1.votePos < v2.votePos;
              });

    int32_t maxClusterPos{0};
    uint32_t maxClusterCount{0};
    double maxClusterScore{0.0};

    // we don't need the return value from the first call
    static_cast<void>(computeBestLoc_(votes, transcript, read, false,
                                      maxClusterPos, maxClusterCount,
                                      maxClusterScore));
    bool revIsBest =
        computeBestLoc_(rcVotes, transcript, read, true, maxClusterPos,
                        maxClusterCount, maxClusterScore);
    isForward_ = not revIsBest;

    bestHitPos = maxClusterPos;
    bestHitCount = maxClusterCount;
    bestHitScore = maxClusterScore;
    return true;
  }

  bool isForward() { return isForward_; }
};

template <typename AlnT>
void processMiniBatch(BulkReadExpT& readExp, ForgettingMassCalculator& fmCalc,
                      uint64_t firstTimestepOfRound, ReadLibrary& readLib,
                      const SalmonOpts& salmonOpts,
                      AlnGroupVecRange<AlnT> batchHits,
                      std::vector<Transcript>& transcripts,
                      ClusterForest& clusterForest,
                      FragmentLengthDistribution& fragLengthDist,
                      BiasParams& observedGCParams,
                      /**
                       * NOTE : test new el model in future
                       * EffectiveLengthStats& effLengthStats,
                       **/
                      std::atomic<uint64_t>& numAssignedFragments,
                      std::default_random_engine& randEng, bool initialRound,
                      std::atomic<bool>& burnedIn, double& maxZeroFrac);

template <typename CoverageCalculator>
inline void collectHitsForRead(SalmonIndex* sidx, const bwtintv_v* a,
                               smem_aux_t* auxHits, mem_opt_t* memOptions,
                               const SalmonOpts& salmonOpts,
                               const uint8_t* read, uint32_t readLen,
                               std::vector<CoverageCalculator>& hits) {
  // std::unordered_map<uint64_t, CoverageCalculator>& hits) {

  bwaidx_t* idx = sidx->bwaIndex();
  mem_collect_intv(salmonOpts, memOptions, sidx, readLen, read, auxHits);

  // For each MEM
  int firstSeedLen{-1};
  for (decltype(auxHits->mem.n) i = 0; i < auxHits->mem.n; ++i) {
    // A pointer to the interval of the MEMs occurences
    bwtintv_t* p = &auxHits->mem.a[i];
    // The start and end positions in the query string (i.e. read) of the MEM
    int qstart = p->info >> 32;
    uint32_t qend = static_cast<uint32_t>(p->info);
    int step, count, slen = (qend - qstart); // seed length

    /*
    if (firstSeedLen > -1) {
        if (slen < firstSeedLen) { return; }
    } else {
        firstSeedLen = slen;
    }
    */

    int64_t k;
    step = (p->x[2] > static_cast<bwtint_t>(memOptions->max_occ)) ? p->x[2] / memOptions->max_occ : 1;
    // For every occurrence of the MEM
    for (k = count = 0;
         k < static_cast<decltype(k)>(p->x[2]) && count < static_cast<decltype(count)>(memOptions->max_occ);
         k += step, ++count) {
      //bwtint_t pos;
      bwtint_t startPos, endPos;
      //int len, isRev, isRevStart, isRevEnd, refID, refIDStart, refIDEnd;
      int isRev, isRevStart, isRevEnd, refID, refIDStart, refIDEnd;
      int queryStart = qstart;
      //len = slen;
      uint32_t rlen = readLen;

      // Get the position in the reference index of this MEM occurrence
      int64_t refStart = bwt_sa(idx->bwt, p->x[0] + k);

      //pos = startPos = bns_depos(idx->bns, refStart, &isRevStart);
      startPos = bns_depos(idx->bns, refStart, &isRevStart);
      endPos = bns_depos(idx->bns, refStart + slen - 1, &isRevEnd);
      // If we span the forward/reverse boundary, discard the hit
      if (isRevStart != isRevEnd) {
        continue;
      }
      // Otherwise, isRevStart = isRevEnd so just assign isRev = isRevStart
      isRev = isRevStart;

      // If the hit is reversed --- swap the start and end
      if (isRev) {
        if (endPos > startPos) {
          salmonOpts.jointLog->warn("Hit is supposedly reversed, "
                                    "but startPos = {} < endPos = {}",
                                    startPos, endPos);
        }
        auto temp = startPos;
        startPos = endPos;
        endPos = temp;
      }
      // Get the ID of the reference sequence in which it occurs
      refID = refIDStart = bns_pos2rid(idx->bns, startPos);
      refIDEnd = bns_pos2rid(idx->bns, endPos);

      if (refID < 0) {
        continue;
      } // bridging multiple reference sequences or the forward-reverse
        // boundary;

      auto tlen = idx->bns->anns[refID].len;

      // The refence sequence-relative (e.g. transcript-relative) position of
      // the MEM
      long hitLoc = static_cast<long>(isRev ? endPos : startPos) -
                    idx->bns->anns[refID].offset;

      if ((refIDStart != refIDEnd)) {
        // If a seed spans two transcripts

        // If we're not considering splitting such seeds, then
        // just discard this seed and continue.
        if (not salmonOpts.splitSpanningSeeds) {
          continue;
        }

        // std::cerr << "Seed spans two transcripts! --- attempting to split:
        // \n";
        if (!isRev) {
          // If it's going forward, we have a situation like this
          // packed transcripts: t1 ===========|t2|==========>
          // hit:                          |==========>

          // length of hit in t1
          auto len1 = tlen - hitLoc;
          // length of hit in t2
          auto len2 = slen - len1;
          if (std::max(len1, len2) < memOptions->min_seed_len) {
            continue;
          }

          /** Keeping this here for now in case I need to debug splitting seeds
          again std::cerr << "\t hit is in the forward direction: "; std::cerr
          << "t1 part has length " << len1 << ", t2 part has length " << len2 <<
          "\n";
          */

          // If the part in t1 is larger then just cut off the rest
          if (len1 >= len2) {
            slen = len1;
            int32_t votePos = static_cast<int32_t>(hitLoc) - queryStart;
            // std::cerr << "\t\t t1 (of length " << tlen << ") has larger hit
            // --- new hit length = " << len1 << "; starts at pos " <<
            // queryStart
            // << " in the read (votePos will be " << votePos << ")\n";
          } else {
            // Otherwise, make the hit be in t2.
            // Because the hit spans the boundary where t2 begins,
            // the new seed begins matching at position 0 of
            // transcript t2
            hitLoc = 0;
            slen = len2;
            // The seed originally started at position q, now it starts  len1
            // characters to the  right of that
            queryStart += len1;
            refID = refIDEnd;
            int32_t votePos = static_cast<int32_t>(hitLoc) - queryStart;
            tlen = idx->bns->anns[refID].len;
            // std::cerr << "\t\t t2 (of length " << tlen << ") has larger hit
            // --- new hit length = " << len2 << "; starts at pos " <<
            // queryStart
            // << " in the read (votePos will be " << votePos << ")\n";
          }
        } else {

          // If it's going in the reverse direction, we have a situation like
          // this packed transcripts: t1 <===========|t2|<========== hit:
          // X======Y>======Z> Which means we have packed transcripts: t1
          // <===========|t2|<========== hit:
          // <Z=====Y<======X length of hit in t1

          auto len2 = endPos - idx->bns->anns[refIDEnd].offset;
          auto len1 = slen - len2;
          if (std::max(len1, len2) < static_cast<decltype(len1)>(memOptions->min_seed_len)) {
            continue;
          }

          /** Keeping this here for now in case I need to debug splitting seeds
          again std::cerr << "\t hit is in the reverse direction: "; std::cerr
          << "\n\n"; std::cerr << "startPos = " << startPos << ", endPos = " <<
          endPos << ", offset[refIDStart] = "
                    <<  idx->bns->anns[refIDStart].offset << ", offset[refIDEnd]
          = " << idx->bns->anns[refIDEnd].offset << "\n"; std::cerr << "\n\n";
          std::cerr << "t1 part has length " << len1 << ", t2 part has length "
          << len2 << "\n\n";
          */

          if (len1 >= len2) {
            slen = len1;
            hitLoc = tlen - len2;
            queryStart += len2;
            rlen -= len2;
            int32_t votePos =
                static_cast<int32_t>(hitLoc) - (rlen - queryStart);
            // std::cerr << "\t\t t1 (hitLoc: " << hitLoc << ") (of length " <<
            // tlen << ") has larger hit --- new hit length = " << len1 << ";
            // starts at pos " << queryStart << " in the read (votePos will be "
            // << votePos << ")\n";
          } else {
            slen = len2;
            refID = bns_pos2rid(idx->bns, endPos);
            tlen = idx->bns->anns[refID].len;
            hitLoc = len2;
            rlen = hitLoc + queryStart;
            int32_t votePos =
                static_cast<int32_t>(hitLoc) - (rlen - queryStart);
            // std::cerr << "\t\t t2 (of length " << tlen << ") (hitLoc: " <<
            // hitLoc << ") has larger hit --- new hit length = " << len2 << ";
            // starts at pos " << queryStart << " in the read (votePos will be "
            // << votePos << ")\n";
          }
        }
      }

      auto hitIt = std::find_if(hits.begin(), hits.end(),
                                [refID](CoverageCalculator& c) -> bool {
                                  return (c.targetID) == static_cast<decltype(c.targetID)>(refID);
                                });
      if (isRev) {
        if (hitIt == hits.end()) {
          CoverageCalculator hit;
          hit.targetID = refID;
          hit.addFragMatchRC(hitLoc, queryStart, slen, rlen);
          hits.emplace_back(hit);
        } else {
          hitIt->addFragMatchRC(hitLoc, queryStart, slen, rlen);
          // hits[refID].addFragMatchRC(hitLoc, queryStart , slen, rlen);
        }
      } else {
        if (hitIt == hits.end()) {
          CoverageCalculator hit;
          hit.targetID = refID;
          hit.addFragMatch(hitLoc, queryStart, slen);
          hits.emplace_back(hit);
        } else {
          hitIt->addFragMatch(hitLoc, queryStart, slen);
          // hits[refID].addFragMatch(hitLoc, queryStart, slen);
        }
      }
    } // for k
  }
}

/*
constexpr inline bool consistentNames(header_sequence_qual& r) { return true; }

constexpr inline bool consistentNames(
    std::pair<header_sequence_qual, header_sequence_qual>& rp) {
  auto l1 = rp.first.header.length();
  auto l2 = rp.second.header.length();
  char* sptr = static_cast<char*>(memchr(&rp.first.header[0], ' ', l1));

  bool compat = false;
  // If we didn't find a space in the name of read1
  if (sptr == NULL) {
    if (l1 > 1) {
      compat = (l1 == l2);
      compat = compat and
               (memcmp(&rp.first.header[0], &rp.second.header[0], l1 - 1) == 0);
      compat =
          compat and ((rp.first.header[l1 - 1] == '1' and
                       rp.second.header[l2 - 1] == '2') or
                      (rp.first.header[l1 - 1] == rp.second.header[l2 - 1]));
    } else {
      compat = (l1 == l2);
      compat = compat and (rp.first.header[0] == rp.second.header[0]);
    }
  } else {
    size_t offset = sptr - (&rp.first.header[0]);

    // If read2 matches read1 up to and including the space
    if (offset + 1 < l2) {
      compat = memcmp(&rp.first.header[0], &rp.second.header[0], offset) == 0;
      // and after the space, read1 and read2 have an identical character or
      // read1 has a '1' and read2 has a '2', then this is a consistent pair.
      compat = compat and
               ((rp.first.header[offset + 1] == rp.second.header[offset + 1]) or
                (rp.first.header[offset + 1] == '1' and
                 rp.second.header[offset + 1] == '2'));
    } else {
      compat = false;
    }
  }
  return compat;
}
*/

/**
 *  Returns true if the @hit is within @cutoff bases of the end of
 *  transcript @txp and false otherwise.
 */
template <typename CoverageCalculator>
inline bool
nearEndOfTranscript(CoverageCalculator& hit, Transcript& txp,
                    int32_t cutoff = std::numeric_limits<int32_t>::max()) {
  // check if hit appears close to the end of the given transcript
  bool isForward = hit.isForward();
  int32_t hitPos = static_cast<int32_t>(hit.bestHitPos);
  return (hitPos <= cutoff or
          std::abs(static_cast<int32_t>(txp.RefLength) - hitPos) <= cutoff);
}

template <typename CoverageCalculator>
inline void getHitsForFragment(
    fastx_parser::ReadPair& frag,
    // std::pair<header_sequence_qual, header_sequence_qual>& frag,
    SalmonIndex* sidx, smem_i* itr, const bwtintv_v* a, smem_aux_t* auxHits,
    mem_opt_t* memOptions, BulkReadExpT& readExp,
    const SalmonOpts& salmonOpts, double coverageThresh,
    uint64_t& upperBoundHits, AlignmentGroup<SMEMAlignment>& hitList,
    uint64_t& hitListCount, std::vector<Transcript>& transcripts) {

  // std::unordered_map<uint64_t, CoverageCalculator> leftHits;
  // std::unordered_map<uint64_t, CoverageCalculator> rightHits;

  std::vector<CoverageCalculator> leftHits;
  std::vector<CoverageCalculator> rightHits;

  uint32_t leftReadLength{0};
  uint32_t rightReadLength{0};

  auto& eqBuilder = readExp.equivalenceClassBuilder();
  bool allowOrphans{salmonOpts.allowOrphans};

  /**
   * As soon as we can decide on an acceptable way to validate read names,
   * we'll inform the user and quit if we see something inconsistent.  However,
   * we first need a reasonable way to verify potential naming formats from
   * many different sources.
   */
  /*
  if (!consistentNames(frag)) {
      fmt::MemoryWriter errstream;

      errstream << "Inconsistent paired-end reads!\n";
      errstream << "mate1 : " << frag.first.header << "\n";
      errstream << "mate2 : " << frag.second.header << "\n";
      errstream << "Paired-end reads should appear consistently in their
  respective files.\n"; errstream << "Please fix the paire-end input before
  quantifying with salmon; exiting.\n";

      std::cerr << errstream.str();
      std::exit(-1);
  }
  */

  //---------- End 1 ----------------------//
  {
    std::string readStr = frag.first.seq;
    uint32_t readLen = readStr.size();

    leftReadLength = readLen;

    for (decltype(readLen) p = 0; p < readLen; ++p) {
      readStr[p] = nst_nt4_table[static_cast<int>(readStr[p])];
    }

    collectHitsForRead(sidx, a, auxHits, memOptions, salmonOpts,
                       reinterpret_cast<const uint8_t*>(readStr.c_str()),
                       readLen, leftHits);
  }

  //---------- End 2 ----------------------//
  {
    std::string readStr = frag.second.seq;
    uint32_t readLen = readStr.size();

    rightReadLength = readLen;

    for (decltype(readLen) p = 0; p < readLen; ++p) {
      readStr[p] = nst_nt4_table[static_cast<int>(readStr[p])];
    }

    collectHitsForRead(sidx, a, auxHits, memOptions, salmonOpts,
                       reinterpret_cast<const uint8_t*>(readStr.c_str()),
                       readLen, rightHits);
  } // end right

  size_t numTrivialHits = (leftHits.size() + rightHits.size() > 0) ? 1 : 0;
  upperBoundHits += (leftHits.size() + rightHits.size() > 0) ? 1 : 0;
  size_t readHits{0};
  auto& alnList = hitList.alignments();
  hitList.isUniquelyMapped() = true;
  alnList.clear();
  // nothing more to do
  if (numTrivialHits == 0) {
    return;
  }

  double cutoffLeft{coverageThresh};  //* leftReadLength};
  double cutoffRight{coverageThresh}; //* rightReadLength};

  uint64_t leftHitCount{0};

  // Fraction of the optimal coverage that a lightweight alignment
  // must obtain in order to be retained.
  float fOpt{0.95};

  // First, see if there are transcripts where both ends of the
  // fragments map
  auto& minHitList =
      (leftHits.size() < rightHits.size()) ? leftHits : rightHits;
  auto& maxHitList =
      (leftHits.size() < rightHits.size()) ? rightHits : leftHits;

  struct JointHitPtr {
    uint32_t transcriptID;
    size_t leftIndex;
    size_t rightIndex;
  };

  std::vector<JointHitPtr> jointHits; // haha (variable name)!
  jointHits.reserve(minHitList.size());

  // vector-based code
  // Sort the left and right hits
  std::sort(
      leftHits.begin(), leftHits.end(),
      [](const CoverageCalculator& c1, const CoverageCalculator& c2) -> bool {
        return c1.targetID < c2.targetID;
      });
  std::sort(
      rightHits.begin(), rightHits.end(),
      [](const CoverageCalculator& c1, const CoverageCalculator& c2) -> bool {
        return c1.targetID < c2.targetID;
      });
  // Take the intersection of these two hit lists
  // Adopted from : http://en.cppreference.com/w/cpp/algorithm/set_intersection
  {
    auto leftIt = leftHits.begin();
    auto leftEnd = leftHits.end();
    auto rightIt = rightHits.begin();
    auto rightEnd = rightHits.end();
    while (leftIt != leftEnd && rightIt != rightEnd) {
      if (leftIt->targetID < rightIt->targetID) {
        ++leftIt;
      } else {
        if (!(rightIt->targetID < leftIt->targetID)) {
          jointHits.push_back(
              {leftIt->targetID,
               static_cast<size_t>(std::distance(leftHits.begin(), leftIt)),
               static_cast<size_t>(std::distance(rightHits.begin(), rightIt))});
          ++leftIt;
        }
        ++rightIt;
      }
    }
  }
  // End vector-based code

  /* map based code
  {
      auto notFound = maxHitList.end();
      for (auto& kv : minHitList) {
          uint64_t refID = kv.first;
          if (maxHitList.find(refID) != notFound) {
              jointHits.emplace_back(refID);
          }
      }
  }
  */

  // Check if the fragment generated orphaned
  // lightweight alignments.
  bool isOrphan = (jointHits.size() == 0);

  uint32_t firstTranscriptID = std::numeric_limits<uint32_t>::max();
  double bestScore = -std::numeric_limits<double>::max();
  bool sortedByTranscript = true;
  int32_t lastTranscriptId = std::numeric_limits<int32_t>::min();

  if (BOOST_UNLIKELY(isOrphan and allowOrphans)) {
    // std::vector<CoverageCalculator> allHits;
    // allHits.reserve(totalHits);

    // search for a hit on the left
    for (auto& tHitList : leftHits) {
      auto transcriptID = tHitList.targetID;
      auto& covChain = tHitList;
      Transcript& t = transcripts[transcriptID];
      if (!t.hasAnchorFragment()) {
        continue;
      }

      covChain.computeBestChain(t, frag.first.seq);
      double score = covChain.bestHitScore;

      // make sure orphaned fragment is near the end of the transcript
      // if (!nearEndOfTranscript(covChain, t, 1000)) { continue; }

      if (score >= fOpt * bestScore and score >= cutoffLeft) {

        if (score > bestScore) {
          bestScore = score;
        }
        bool isForward = covChain.isForward();
        int32_t hitPos = covChain.bestHitPos;
        auto fmt = salmon::utils::hitType(hitPos, isForward);

        if (leftHitCount == 0) {
          firstTranscriptID = transcriptID;
        } else if (hitList.isUniquelyMapped() and
                   transcriptID != firstTranscriptID) {
          hitList.isUniquelyMapped() = false;
        }

        if (static_cast<decltype(lastTranscriptId)>(transcriptID) < lastTranscriptId) {
          sortedByTranscript = false;
        }

        alnList.emplace_back(transcriptID, fmt, score, hitPos);
        alnList.back().fwd = isForward;
        alnList.back().mateStatus = rapmap::utils::MateStatus::PAIRED_END_LEFT;
        readHits += score;
        ++hitListCount;
        ++leftHitCount;
      }
    }

    // search for a hit on the right
    for (auto& tHitList : rightHits) {
      // Prior
      // auto transcriptID = tHitList.first;
      auto transcriptID = tHitList.targetID;
      auto& covChain = tHitList;
      Transcript& t = transcripts[transcriptID];
      if (!t.hasAnchorFragment()) {
        continue;
      }

      covChain.computeBestChain(t, frag.second.seq);
      double score = covChain.bestHitScore;

      // make sure orphaned fragment is near the end of the transcript
      // if (!nearEndOfTranscript(covChain, t, 1000)) { continue; }

      if (score >= fOpt * bestScore and score >= cutoffRight) {
        if (score > bestScore) {
          bestScore = score;
        }
        bool isForward = covChain.isForward();
        int32_t hitPos = covChain.bestHitPos;
        auto fmt = salmon::utils::hitType(hitPos, isForward);
        if (leftHitCount == 0) {
          firstTranscriptID = transcriptID;
        } else if (hitList.isUniquelyMapped() and
                   transcriptID != firstTranscriptID) {
          hitList.isUniquelyMapped() = false;
        }

        alnList.emplace_back(transcriptID, fmt, score, hitPos);
        alnList.back().fwd = isForward;
        alnList.back().mateStatus = rapmap::utils::MateStatus::PAIRED_END_RIGHT;
        readHits += score;
        ++hitListCount;
        ++leftHitCount;
      }
    }

    if (alnList.size() > 0) {
      auto newEnd =
          std::stable_partition(alnList.begin(), alnList.end(),
                                [bestScore, fOpt](SMEMAlignment& aln) -> bool {
                                  return aln.score() >= fOpt * bestScore;
                                });
      alnList.resize(std::distance(alnList.begin(), newEnd));
      if (!sortedByTranscript) {
        std::sort(alnList.begin(), alnList.end(),
                  [](const SMEMAlignment& x, const SMEMAlignment& y) -> bool {
                    return x.transcriptID() < y.transcriptID();
                  });
      }
    } else {
      return;
      /*
      // If we didn't have any *significant* hits --- add any *trivial* orphan
      hits size_t totalHits = leftHits.size() + rightHits.size();
      std::vector<uint32_t> txpIDs;
      txpIDs.reserve(totalHits);
      std::vector<double> auxProbs;
      auxProbs.reserve(totalHits);

      size_t txpIDsHash{0};
      std::vector<CoverageCalculator> allHits;
      allHits.reserve(totalHits);
      std::merge(leftHits.begin(), leftHits.end(),
                 rightHits.begin(), rightHits.end(),
                 std::back_inserter(allHits),
                 [](CoverageCalculator& c1, CoverageCalculator& c2) -> bool {
                  return c1.targetID < c2.targetID;
                 });
      double totProb{0.0};
      for (auto& h : allHits) {
          boost::hash_combine(txpIDsHash, h.targetID);
          txpIDs.push_back(h.targetID);
          double refLen =  std::max(1.0,
      static_cast<double>(transcripts[h.targetID].RefLength)); double startProb
      = 1.0 / refLen; auxProbs.push_back(startProb); totProb += startProb;
      }
      if (totProb > 0.0) {
          double norm = 1.0 / totProb;
          for (auto& p : auxProbs) { p *= norm; }

          TranscriptGroup tg(txpIDs, txpIDsHash);
          eqBuilder.addGroup(std::move(tg), auxProbs);
      } else {
          salmonOpts.jointLog->warn("Unexpected empty hit group [orphaned]");
      }
      */
    }
  } else { // Not an orphan
    for (auto jhp : jointHits) {
      auto& jointHitPtr = jhp;
      auto transcriptID = jhp.transcriptID;
      Transcript& t = transcripts[transcriptID];
      auto& leftHitList = leftHits[jhp.leftIndex];
      leftHitList.computeBestChain(t, frag.first.seq);
      if (leftHitList.bestHitScore >= cutoffLeft) {
        auto& rightHitList = rightHits[jhp.rightIndex];

        rightHitList.computeBestChain(t, frag.second.seq);
        if (rightHitList.bestHitScore < cutoffRight) {
          continue;
        }

        auto end1Start = leftHitList.bestHitPos;
        auto end2Start = rightHitList.bestHitPos;

        double score =
            (leftHitList.bestHitScore + rightHitList.bestHitScore) * 0.5;
        if (score < fOpt * bestScore) {
          continue;
        }

        if (score > bestScore) {
          bestScore = score;
        }

        uint32_t fragLength = std::abs(static_cast<int32_t>(end1Start) -
                                       static_cast<int32_t>(end2Start)) +
                              rightReadLength;

        bool end1IsForward = leftHitList.isForward();
        bool end2IsForward = rightHitList.isForward();

        uint32_t end1Pos = (end1IsForward)
                               ? leftHitList.bestHitPos
                               : leftHitList.bestHitPos + leftReadLength;
        uint32_t end2Pos = (end2IsForward)
                               ? rightHitList.bestHitPos
                               : rightHitList.bestHitPos + rightReadLength;
        bool canDovetail = false;
        auto fmt = salmon::utils::hitType(
            end1Pos, end1IsForward, leftReadLength, end2Pos, end2IsForward,
            rightReadLength, canDovetail);

        if (readHits == 0) {
          firstTranscriptID = transcriptID;
        } else if (hitList.isUniquelyMapped() and
                   transcriptID != firstTranscriptID) {
          hitList.isUniquelyMapped() = false;
        }

        int32_t minHitPos = std::min(end1Pos, end2Pos);
        if (static_cast<decltype(lastTranscriptId)>(transcriptID) < lastTranscriptId) {
          sortedByTranscript = false;
        }
        // ANCHOR TEST
        t.setAnchorFragment();
        alnList.emplace_back(transcriptID, fmt, score, minHitPos, fragLength);
        alnList.back().fwd = end1IsForward;
        alnList.back().mateIsFwd = end2IsForward;
        alnList.back().mateStatus =
            rapmap::utils::MateStatus::PAIRED_END_PAIRED;
        ++readHits;
        ++hitListCount;
      }
    } // end for jointHits
    if (alnList.size() > 0) {
      auto newEnd =
          std::stable_partition(alnList.begin(), alnList.end(),
                                [bestScore, fOpt](SMEMAlignment& aln) -> bool {
                                  return aln.score() >= fOpt * bestScore;
                                });
      alnList.resize(std::distance(alnList.begin(), newEnd));
      if (!sortedByTranscript) {
        std::sort(alnList.begin(), alnList.end(),
                  [](const SMEMAlignment& x, const SMEMAlignment& y) -> bool {
                    return x.transcriptID() < y.transcriptID();
                  });
      }
    } else {
      // If we didn't have any *significant* hits --- add any *trivial* joint
      // hits
      return;
      /*
      std::vector<uint32_t> txpIDs;
      txpIDs.reserve(jointHits.size());
      std::vector<double> auxProbs;
      auxProbs.reserve(jointHits.size());

      size_t txpIDsHash{0};
      double totProb{0.0};
      for (auto& h : jointHits) {
          boost::hash_combine(txpIDsHash, h.transcriptID);
          txpIDs.push_back(h.transcriptID);
          double refLen =  std::max(1.0,
      static_cast<double>(transcripts[h.transcriptID].RefLength)); double
      startProb = 1.0 / refLen; auxProbs.push_back(startProb); totProb +=
      startProb;
      }
      if (totProb > 0.0) {
      double norm = 1.0 / totProb;
      for (auto& p : auxProbs) { p *= norm; }

      TranscriptGroup tg(txpIDs, txpIDsHash);
      eqBuilder.addGroup(std::move(tg), auxProbs);
      } else {
          salmonOpts.jointLog->warn("Unexpected empty hit group [paired]");
      }
      */
    }

  } // end else
}

/**
 *   Get hits for single-end fragment
 *
 *
 */
template <typename CoverageCalculator>
inline void getHitsForFragment(fastx_parser::ReadSeq& frag,
                               // jellyfish::header_sequence_qual& frag,
                               SalmonIndex* sidx, smem_i* itr,
                               const bwtintv_v* a, smem_aux_t* auxHits,
                               mem_opt_t* memOptions, BulkReadExpT& readExp,
                               const SalmonOpts& salmonOpts,
                               double coverageThresh, uint64_t& upperBoundHits,
                               AlignmentGroup<SMEMAlignment>& hitList,
                               uint64_t& hitListCount,
                               std::vector<Transcript>& transcripts) {

  uint64_t leftHitCount{0};

  // std::unordered_map<uint64_t, CoverageCalculator> hits;
  std::vector<CoverageCalculator> hits;

  auto& eqBuilder = readExp.equivalenceClassBuilder();

  //uint32_t readLength{0};

  //---------- get hits ----------------------//
  {
    std::string readStr = frag.seq;
    uint32_t readLen = frag.seq.size();

    //readLength = readLen;

    for (int p = 0; p < static_cast<int>(readLen); ++p) {
      readStr[p] = nst_nt4_table[static_cast<int>(readStr[p])];
    }

    char* readPtr = const_cast<char*>(readStr.c_str());

    collectHitsForRead(sidx, a, auxHits, memOptions, salmonOpts,
                       reinterpret_cast<const uint8_t*>(readStr.c_str()),
                       readLen, hits);
  }

  upperBoundHits += (hits.size() > 0) ? 1 : 0;

  int32_t lastTranscriptId = std::numeric_limits<int32_t>::min();
  bool sortedByTranscript{true};
  double fOpt{0.95};
  double bestScore = -std::numeric_limits<double>::max();

  size_t readHits{0};
  auto& alnList = hitList.alignments();
  hitList.isUniquelyMapped() = true;
  alnList.clear();

  uint32_t firstTranscriptID = std::numeric_limits<uint32_t>::max();
  double cutoff{coverageThresh}; //* readLength};
  for (auto& tHitList : hits) {
    // Prior
    // auto hitID = tHitList.first;
    // auto& covVec = tHitList.second;
    auto hitID = tHitList.targetID;
    auto& covVec = tHitList;

    // Coverage score
    Transcript& t = transcripts[hitID];
    covVec.computeBestChain(t, frag.seq);
    double score = covVec.bestHitScore;
    if (score >= fOpt * bestScore and covVec.bestHitScore >= cutoff) {

      bool isForward = covVec.isForward();
      if (score < fOpt * bestScore) {
        continue;
      }

      if (score > bestScore) {
        bestScore = score;
      }

      auto hitPos = covVec.bestHitPos;
      auto fmt = salmon::utils::hitType(hitPos, isForward);

      if (leftHitCount == 0) {
        firstTranscriptID = hitID;
      } else if (hitList.isUniquelyMapped() and hitID != firstTranscriptID) {
        hitList.isUniquelyMapped() = false;
      }

      auto transcriptID = hitID;

      if (static_cast<int32_t>(transcriptID) < lastTranscriptId) {
        sortedByTranscript = false;
      }

      alnList.emplace_back(transcriptID, fmt, score, hitPos);
      alnList.back().fwd = isForward;
      alnList.back().mateStatus = rapmap::utils::MateStatus::SINGLE_END;
      readHits += score;
      ++hitListCount;
      ++leftHitCount;
    }
  }
  if (alnList.size() > 0) {
    auto newEnd =
        std::stable_partition(alnList.begin(), alnList.end(),
                              [bestScore, fOpt](SMEMAlignment& aln) -> bool {
                                return aln.score() >= fOpt * bestScore;
                              });
    alnList.resize(std::distance(alnList.begin(), newEnd));
    if (!sortedByTranscript) {
      std::sort(alnList.begin(), alnList.end(),
                [](const SMEMAlignment& x, const SMEMAlignment& y) -> bool {
                  return x.transcriptID() < y.transcriptID();
                });
    }
  } else {
    // If we didn't have any *significant* hits --- add any *trivial* joint hits
    return;
    /*
    std::vector<uint32_t> txpIDs;
    txpIDs.reserve(hits.size());
    double uniProb = 1.0 / hits.size();
    std::vector<double> auxProbs(hits.size(), uniProb);

    size_t txpIDsHash{0};
    for (auto& h : hits) {
        boost::hash_combine(txpIDsHash, h.targetID);
        txpIDs.push_back(h.targetID);
    }

    TranscriptGroup tg(txpIDs, txpIDsHash);
    eqBuilder.addGroup(std::move(tg), auxProbs);
    */
  }
}

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename ParserT, typename CoverageCalculator>
void processReadsMEM(
    ParserT* parser, BulkReadExpT& readExp, ReadLibrary& rl,
    AlnGroupVec<QuasiAlignment>& structureVec,
    std::atomic<uint64_t>& numObservedFragments,
    std::atomic<uint64_t>& numAssignedFragments,
    std::atomic<uint64_t>& validHits, std::atomic<uint64_t>& upperBoundHits,
    SalmonIndex* sidx, std::vector<Transcript>& transcripts,
    ForgettingMassCalculator& fmCalc, ClusterForest& clusterForest,
    FragmentLengthDistribution& fragLengthDist, BiasParams& observedGCParams,
    /**
     * NOTE : test new el model in future
     * EffectiveLengthStats& obsEffLength,
     **/
    mem_opt_t* memOptions, const SalmonOpts& salmonOpts, double coverageThresh,
    std::mutex& iomutex, bool initialRound, std::atomic<bool>& burnedIn,
    volatile bool& writeToCache) {
  // ERROR
  salmonOpts.jointLog->error("Quasimapping cannot be used with the FMD index "
                             "--- please report this bug on GitHub");
  std::exit(1);
}

template <typename ParserT, typename CoverageCalculator>
void processReadsMEM(
    ParserT* parser, BulkReadExpT& readExp, ReadLibrary& rl,
    AlnGroupVec<SMEMAlignment>& structureVec,
    std::atomic<uint64_t>& numObservedFragments,
    std::atomic<uint64_t>& numAssignedFragments,
    std::atomic<uint64_t>& validHits, std::atomic<uint64_t>& upperBoundHits,
    SalmonIndex* sidx, std::vector<Transcript>& transcripts,
    ForgettingMassCalculator& fmCalc, ClusterForest& clusterForest,
    FragmentLengthDistribution& fragLengthDist, BiasParams& observedGCParams,
    /**
     * NOTE : test new el model in future
     * EffectiveLengthStats& obsEffLengths,
     **/
    mem_opt_t* memOptions, const SalmonOpts& salmonOpts, double coverageThresh,
    std::mutex& iomutex, bool initialRound, std::atomic<bool>& burnedIn,
    volatile bool& writeToCache) {
  uint64_t count_fwd = 0, count_bwd = 0;
  // Seed with a real random value, if available
  std::random_device rd;

  // Create a random uniform distribution
  std::default_random_engine eng(rd());

  uint64_t prevObservedFrags{1};
  uint64_t leftHitCount{0};
  uint64_t hitListCount{0};

  // Super-MEM iterator
  smem_i* itr = smem_itr_init(sidx->bwaIndex()->bwt);
  const bwtintv_v* a = nullptr;
  smem_aux_t* auxHits = smem_aux_init();

  //auto expectedLibType = rl.format();

  uint64_t firstTimestepOfRound = fmCalc.getCurrentTimestep();

  size_t locRead{0};
  uint64_t localUpperBoundHits{0};
  size_t rangeSize{0};
  double maxZeroFrac{0.0};
  auto rg = parser->getReadGroup();
  while (parser->refill(rg)) {
    rangeSize = rg.size();
    if (rangeSize > structureVec.size()) {
      salmonOpts.jointLog->error("rangeSize = {}, but structureVec.size() = {} "
                                 "--- this shouldn't happen.\n"
                                 "Please report this bug on GitHub",
                                 rangeSize, structureVec.size());
      std::exit(1);
    }

    for (size_t i = 0; i < rangeSize; ++i) { // For all the read in this batch
      localUpperBoundHits = 0;

      auto& hitList = structureVec[i];
      getHitsForFragment<CoverageCalculator>(
          rg[i],
          // j->data[i],
          sidx, itr, a, auxHits, memOptions, readExp, salmonOpts,
          coverageThresh, localUpperBoundHits, hitList, hitListCount,
          transcripts);
      if (initialRound) {
        upperBoundHits += localUpperBoundHits;
      }

      // If the read mapped to > maxReadOccs places, discard it
      if (hitList.size() > salmonOpts.maxReadOccs) {
        hitList.alignments().clear();
      }
      validHits += hitList.size();
      locRead++;
      ++numObservedFragments;
      if (numObservedFragments % 50000 == 0) {
        iomutex.lock();
        const char RESET_COLOR[] = "\x1b[0m";
        char green[] = "\x1b[30m";
        green[3] = '0' + static_cast<char>(fmt::GREEN);
        char red[] = "\x1b[30m";
        red[3] = '0' + static_cast<char>(fmt::RED);
        if (initialRound) {
          fmt::print(stderr, "\033[A\r\r{}processed{} {} {}fragments{}\n",
                     green, red, numObservedFragments, green, RESET_COLOR);
          fmt::print(stderr, "hits: {}; hits per frag:  {}", validHits,
                     validHits / static_cast<float>(prevObservedFrags));
        } else {
          fmt::print(stderr, "\r\r{}processed{} {} {}fragments{}", green, red,
                     numObservedFragments, green, RESET_COLOR);
        }
        iomutex.unlock();
      }

    } // end for i < j->nb_filled

    prevObservedFrags = numObservedFragments;
    AlnGroupVecRange<SMEMAlignment> hitLists = {structureVec.begin(), structureVec.begin()+rangeSize};
      /*boost::make_iterator_range(
        structureVec.begin(), structureVec.begin() + rangeSize);*/
    processMiniBatch<SMEMAlignment>(
        readExp, fmCalc, firstTimestepOfRound, rl, salmonOpts, hitLists,
        transcripts, clusterForest, fragLengthDist, observedGCParams,
        /**
         * NOTE : test new el model in future
         * obsEffLengths,
         **/
        numAssignedFragments, eng, initialRound, burnedIn, maxZeroFrac);
  }

  if (maxZeroFrac > 0.0) {
    salmonOpts.jointLog->info("Thread saw mini-batch with a maximum of "
                              "{0:.2f}\% zero probability fragments",
                              maxZeroFrac);
  }

  smem_aux_destroy(auxHits);
  smem_itr_destroy(itr);
}

#endif // LIGHTWEIGHT_ALIGNMENT_DEFS_HPP
