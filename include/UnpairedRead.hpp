#ifndef UNPAIRED_READ
#define UNPAIRED_READ

extern "C" {
#ifdef HAVE_CONFIG_H
#undef HAVE_CONFIG_H
#endif

#include "io_lib/scram.h"
#include "io_lib/os.h"

//#define HAVE_CONFIG_H
}

#include "StadenUtils.hpp"
#include "SalmonMath.hpp"
#include "LibraryFormat.hpp"

struct UnpairedRead {
   bam_seq_t* read = nullptr;
   double logProb;
   LibraryFormat libFmt{ReadType::PAIRED_END, ReadOrientation::NONE, ReadStrandedness::U};

   UnpairedRead() : read(staden::utils::bam_init()), logProb(salmon::math::LOG_0) {}

   UnpairedRead(bam_seq_t* r, double lp, LibraryFormat lf) :
       read(r), logProb(lp), libFmt(lf) {}

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

   inline LibraryFormat& libFormat() { return libFmt; }
   inline bool isPaired() const { return false; }
   inline bool isLeftOrphan() const { return false; }
   inline bool isRightOrphan() const { return false; }

    // return 0 on success, -1 on failure
    int writeToFile(scram_fd* fp) {
        return scram_put_seq(fp, read);
    }

    inline char* getName() {
        return  bam_name(read);
    }

    inline uint32_t getNameLength() {
        return bam_name_len(read);
    }

   inline bool isRight() { return bam_flag(read) & BAM_FREVERSE; }
   inline bool isLeft()  { return !isRight(); }
   inline int32_t left() { return bam_pos(read); }
   inline int32_t right() { return left() + bam_seq_len(read); }
   inline uint32_t fragLen() { return 0; }
   inline ReadType fragType() { return ReadType::SINGLE_END; }
   inline int32_t transcriptID() const { return bam_ref(read); }

    inline double logQualProb() {
        return salmon::math::LOG_1;
        /*
        int q = bam_map_qual(read);
        //double logP = (q == 255) ? calcQuality(read) : std::log(std::pow(10.0, -q * 0.1));
        double logP = (q == 255) ? salmon::math::LOG_1 : std::log(std::pow(10.0, -q * 0.1));
        return logP;
        */
    }

};

#endif //UNPAIRED_READ
