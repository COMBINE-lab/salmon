#ifndef UNPAIRED_READ
#define UNPAIRED_READ

//extern "C" {
//#ifdef HAVE_CONFIG_H
//#undef HAVE_CONFIG_H
//#endif
//
//#include "io_lib/os.h"
//#include "io_lib/scram.h"
//
////#define HAVE_CONFIG_H
//}

#include "LibraryFormat.hpp"
#include "RapMapUtils.hpp"
#include "SalmonMath.hpp"
//#include "StadenUtils.hpp"
#include "SamTypes.hpp"

struct UnpairedRead {
  SamRecord* read = nullptr;
  double logProb;
  LibraryFormat libFmt{ReadType::PAIRED_END, ReadOrientation::NONE,
                       ReadStrandedness::U};

  UnpairedRead()
      : read(combinelab::samutils::bam_init()), logProb(salmon::math::LOG_0) {}

  UnpairedRead(SamRecord* r, double lp, LibraryFormat lf)
      : read(r), logProb(lp), libFmt(lf) {}

  UnpairedRead(UnpairedRead&& other) {
    logProb = other.logProb;
    std::swap(read, other.read);
    libFmt = other.libFmt;
  }

  UnpairedRead& operator=(UnpairedRead&& other) {
    logProb = other.logProb;
    std::swap(read, other.read);
    libFmt = other.libFmt;
    return *this;
  }

  UnpairedRead(UnpairedRead& other) = default;

  UnpairedRead& operator=(UnpairedRead& other) = default;

  UnpairedRead* clone() {
    return new UnpairedRead(combinelab::samutils::bam_dup(read), logProb, libFmt);
  }

  ~UnpairedRead() { combinelab::samutils::bam_destroy(read); }

  inline LibraryFormat& libFormat() { return libFmt; }
  inline bool isPaired() const { return false; }
  inline bool isLeftOrphan() const { return false; }
  inline bool isRightOrphan() const { return false; }
  inline SamRecord* get5PrimeRead() { return read; }

  inline rapmap::utils::MateStatus mateStatus() const {
    return rapmap::utils::MateStatus::SINGLE_END;
  }

  inline int32_t pos() const { return left(); }
  inline bool fwd() const { return !combinelab::samutils::bam_strand(read); }
  inline bool isInward() const { return false; }
  // return 0 on success, -1 on failure
  // libstaden
  // int writeToFile(scram_fd* fp) { return scram_put_seq(fp, read); }
  // samtools
  int writeToFile(SamFile* fp, SamHeader* hdr) { return sam_write1(fp, hdr, read); }

  inline char* getName() { return combinelab::samutils::bam_name(read); }

  inline uint32_t getNameLength() { return combinelab::samutils::bam_name_len(read); }

  inline bool isRight() const { return combinelab::samutils::bam_flag(read) & BAM_FREVERSE; }
  inline bool isLeft() const { return !isRight(); }
  inline int32_t left() const { return combinelab::samutils::bam_pos(read); }
  inline int32_t right() const { return left() + combinelab::samutils::bam_seq_len(read); }
  // will always be at least the length of a single read
  inline uint32_t fragLen() const { return 0; }
  // from the leftmost end of the 5' read to the rightmost
  // end of the 3' read (can be less than the length of a single read)
  inline uint32_t fragLengthPedantic(uint32_t txpLen) const { return 0; }
  inline ReadType fragType() const { return ReadType::SINGLE_END; }
  inline int32_t transcriptID() const { return combinelab::samutils::bam_ref(read); }

  inline double logQualProb() const {
    return salmon::math::LOG_1;
    /*
    int q = bam_map_qual(read);
    //double logP = (q == 255) ? calcQuality(read) : std::log(std::pow(10.0, -q
    * 0.1)); double logP = (q == 255) ? salmon::math::LOG_1 :
    std::log(std::pow(10.0, -q * 0.1)); return logP;
    */
  }
};

#endif // UNPAIRED_READ
