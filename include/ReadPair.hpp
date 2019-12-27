#ifndef READ_PAIR
#define READ_PAIR

#include "LibraryFormat.hpp"
//#include "RapMapUtils.hpp"
#include "pufferfish/Util.hpp"
#include "SalmonMath.hpp"
#include "SalmonUtils.hpp"
#include "StadenUtils.hpp"

#include "spdlog/fmt/fmt.h"

struct ReadPair {
  bam_seq_t* read1 = nullptr;
  bam_seq_t* read2 = nullptr;
  salmon::utils::OrphanStatus orphanStatus;
  double logProb;
  LibraryFormat libFmt{ReadType::PAIRED_END, ReadOrientation::NONE,
                       ReadStrandedness::U};

  ReadPair()
      : read1(staden::utils::bam_init()), read2(staden::utils::bam_init()),
        orphanStatus(salmon::utils::OrphanStatus::Paired),
        logProb(salmon::math::LOG_0) {}

  ReadPair(bam_seq_t* r1, bam_seq_t* r2, salmon::utils::OrphanStatus os,
           double lp, LibraryFormat lf)
      : read1(r1), read2(r2), orphanStatus(os), logProb(lp), libFmt(lf) {}

  ReadPair(ReadPair&& other) {
    orphanStatus = other.orphanStatus;
    logProb = other.logProb;
    std::swap(read1, other.read1);
    std::swap(read2, other.read2);
    libFmt = other.libFmt;
  }

  ReadPair& operator=(ReadPair&& other) {
    orphanStatus = other.orphanStatus;
    logProb = other.logProb;
    std::swap(read1, other.read1);
    std::swap(read2, other.read2);
    libFmt = other.libFmt;
    return *this;
  }

  ReadPair(ReadPair& other) = default;

  ReadPair& operator=(ReadPair& other) = default;

  ReadPair* clone() {
    return new ReadPair(bam_dup(read1), bam_dup(read2), orphanStatus, logProb,
                        libFmt);
  }

  ~ReadPair() {
    staden::utils::bam_destroy(read1);
    staden::utils::bam_destroy(read2);
  }

  inline bam_seq_t* getRead1() { return read1; }
  inline bam_seq_t* getRead2() { return read2; }
  inline int32_t readLen() { return bam_seq_len(getRead1()); }
  inline int32_t mateLen() { return bam_seq_len(getRead2()); }

  inline LibraryFormat& libFormat() { return libFmt; }
  inline bool isPaired() const {
    return (orphanStatus == salmon::utils::OrphanStatus::Paired);
  }
  inline bool isLeftOrphan() const {
    return (orphanStatus == salmon::utils::OrphanStatus::LeftOrphan);
  }
  inline bool isRightOrphan() const {
    return (orphanStatus == salmon::utils::OrphanStatus::RightOrphan);
  }
  inline bam_seq_t* get5PrimeRead() {
    return (isPaired() or isLeftOrphan()) ? read1 : nullptr;
  }

  inline pufferfish::util::MateStatus mateStatus() const {
    if (isPaired()) {
      return pufferfish::util::MateStatus::PAIRED_END_PAIRED;
    } else if (isLeftOrphan()) {
      return pufferfish::util::MateStatus::PAIRED_END_LEFT;
    } else if (isRightOrphan()) {
      return pufferfish::util::MateStatus::PAIRED_END_RIGHT;
    }

    std::cerr << "ReadPair.hpp : mateStatus() --- should not get here ";
    std::cerr << "this may be a bug.  Please report it\n";

    return pufferfish::util::MateStatus::PAIRED_END_PAIRED;
  }

  inline bool haveASTag() const {
    uint8_t* tp = reinterpret_cast<uint8_t*>(bam_aux_find(read1, "AS"));
    return !(tp == NULL);
  }

  inline int32_t getAS() const {
    int32_t s{0};
    auto bf1 = bam_flag(read1);
    auto bf2 = bam_flag(read2);
    // if the read is mapped
    if (!(bf1 & BAM_FUNMAP)) {
      uint8_t* tl = reinterpret_cast<uint8_t*>(bam_aux_find(read1, "AS"));
      s += (tl == NULL) ? 0 : bam_aux_i(tl); 
    }
    if (!(bf2 & BAM_FUNMAP)) {
      uint8_t* tr = reinterpret_cast<uint8_t*>(bam_aux_find(read2, "AS"));
      s += (tr == NULL) ? 0 : bam_aux_i(tr);
    }
    return s;
  }

  inline int32_t pos() const { return left(); }
  inline bool fwd() const { return !bam_strand(read1); }
  inline bool isInward() const {
    bool fw1 = !bam_strand(read1);
    bool fw2 = !bam_strand(read2);
    return (fw1 != fw2);
  }

  /**
   * returns 0 on success, -1 on failure.
   */
  int writeToFile(scram_fd* fp) {
    int r1 = scram_put_seq(fp, read1);
    if (r1 == 0 and isPaired()) {
      return scram_put_seq(fp, read2);
    } else {
      return r1;
    }
  }

  inline char* getName() const { return bam_name(read1); }

  inline uint32_t getNameLength() {
    uint32_t l = bam_name_len(read1);
    char* r = getName();
    if (l > 2 and r[l - 2] == '/') {
      return l - 2;
    }
    return l; // bam_name_len(read1);
  }

  // from the leftmost end of the 5' read to the rightmost
  // end of the 3' read (can be less than the length of a single read)
  inline uint32_t fragLengthPedantic(uint32_t txpLen) const {
    int32_t txpLenSigned = static_cast<int32_t>(txpLen);
    if (!isPaired()) {
      return 0;
    }
    bool fw1 = !bam_strand(read1);
    bool fw2 = !bam_strand(read2);
    if (fw1 != fw2) {
      int32_t p1 = fw1 ? bam_pos(read1) : bam_pos(read2);
      p1 = (p1 < 0) ? 0 : p1;
      p1 = (p1 > txpLenSigned) ? txpLenSigned : p1;
      int32_t p2 = fw1 ? bam_pos(read2) + bam_seq_len(read2)
                       : bam_pos(read1) + bam_seq_len(read1);
      p2 = (p2 < 0) ? 0 : p2;
      p2 = (p2 > txpLenSigned) ? txpLenSigned : p2;
      return (p1 > p2) ? p1 - p2 : p2 - p1;
    }

    return 0;
  }

  inline uint32_t fragLen() const {
    if (!isPaired()) {
      return 0;
    }
    auto leftmost1 = bam_pos(read1);
    auto leftmost2 = bam_pos(read2);

    // The length of the mapped read that is "rightmost" w.r.t. the forward
    // strand.
    auto rightmostLen =
        (leftmost1 < leftmost2) ? bam_seq_len(read2) : bam_seq_len(read1);
    return std::abs(leftmost1 - leftmost2) + rightmostLen;

    // return std::abs(read1->core.isize) + std::abs(read1->core.l_qseq) +
    // std::abs(read2->core.l_qseq);
  }

  inline bool isRight() const {
    return isPaired() ? false : (isRightOrphan() ? true : false);
  }
  inline bool isLeft() const {
    return isPaired() ? false : (isLeftOrphan() ? true : false);
  }

  inline int32_t left() const {
    if (isPaired()) {
      return std::min(bam_pos(read1), bam_pos(read2));
    } else {
      return bam_pos(read1);
    }
  }

  inline int32_t right() const {
    if (isPaired()) {
      return std::max(bam_pos(read1) + bam_seq_len(read1),
                      bam_pos(read2) + bam_seq_len(read2));
    } else {
      return bam_pos(read1) + bam_seq_len(read1);
    }
  }

  inline ReadType fragType() const { return ReadType::PAIRED_END; }
  inline int32_t transcriptID() const { return bam_ref(read1); }

  inline double logQualProb() {
    return salmon::math::LOG_1;
    /*
    double q1 = bam_map_qual(read1);
    double q2 = bam_map_qual(read2);
    double logP1 = (q1 == 255) ? salmon::math::LOG_1 : std::log(std::pow(10.0,
    -q1 * 0.1)); double logP2 = (q2 == 255) ? salmon::math::LOG_1 :
    std::log(std::pow(10.0, -q2 * 0.1)); return logP1 + logP2;
    */
  }
};

#endif // READ_PAIR
