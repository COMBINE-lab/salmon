#ifndef UNPAIRED_READ
#define UNPAIRED_READ

extern "C" {
#ifdef HAVE_CONFIG_H
#undef HAVE_CONFIG_H
#endif

#include "io_lib/os.h"
#include "io_lib/scram.h"

//#define HAVE_CONFIG_H
}

#include "LibraryFormat.hpp"
//#include "RapMapUtils.hpp"
#include "pufferfish/Util.hpp"
#include "SalmonMath.hpp"
#include "StadenUtils.hpp"

struct UnpairedRead {
  bam_seq_t* read = nullptr;
  double logProb;
  LibraryFormat libFmt{ReadType::SINGLE_END, ReadOrientation::NONE,
                       ReadStrandedness::U};

  UnpairedRead()
      : read(staden::utils::bam_init()), logProb(salmon::math::LOG_0) {}

  UnpairedRead(bam_seq_t* r, double lp, LibraryFormat lf)
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
    return new UnpairedRead(bam_dup(read), logProb, libFmt);
  }

  ~UnpairedRead() { staden::utils::bam_destroy(read); }

  inline bam_seq_t* getRead1() { return read; }
  inline bam_seq_t* getRead2() { return read; }
  inline int32_t readLen() { return bam_seq_len(getRead1()); }
  inline int32_t mateLen() { return bam_seq_len(getRead2()); }

  inline LibraryFormat& libFormat() { return libFmt; }
  inline bool isPaired() const { return false; }
  inline bool isLeftOrphan() const { return false; }
  inline bool isRightOrphan() const { return false; }
  inline bam_seq_t* get5PrimeRead() { return read; }

  inline pufferfish::util::MateStatus mateStatus() const {
    return pufferfish::util::MateStatus::SINGLE_END;
  }

  inline bool haveASTag() const { 
    uint8_t* tp = reinterpret_cast<uint8_t*>(bam_aux_find(read, "AS"));
    return !(tp == NULL);
  }

  inline int32_t getAS() const {
    uint8_t* tp = reinterpret_cast<uint8_t*>(bam_aux_find(read, "AS"));
    return (tp == NULL) ? 0 : bam_aux_i(tp);
  }

  inline int32_t pos() const { return left(); }
  inline bool fwd() const { return !bam_strand(read); }
  inline bool isInward() const { return false; }
  // return 0 on success, -1 on failure
  int writeToFile(scram_fd* fp) { return scram_put_seq(fp, read); }

  inline char* getName() { return bam_name(read); }

  inline uint32_t getNameLength() { return bam_name_len(read); }

  inline bool isRight() const { return bam_flag(read) & BAM_FREVERSE; }
  inline bool isLeft() const { return !isRight(); }
  inline int32_t left() const { return bam_pos(read); }
  inline int32_t right() const { return left() + bam_seq_len(read); }
  // will always be at least the length of a single read
  inline uint32_t fragLen() const { return 0; }
  // from the leftmost end of the 5' read to the rightmost
  // end of the 3' read (can be less than the length of a single read)
  inline uint32_t fragLengthPedantic(uint32_t /*txpLen*/) const { return 0; }
  inline ReadType fragType() const { return ReadType::SINGLE_END; }
  inline int32_t transcriptID() const { return bam_ref(read); }

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
