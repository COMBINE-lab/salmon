#ifndef __SAMTYPES_HPP__
#define __SAMTYPES_HPP__

#include <memory>

extern "C" {
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/klist.h"
#include "htslib/thread_pool.h"
#include "htslib/bgzf.h"
}

#include "spdlog/spdlog.h"

using SamFile = samFile; // scram_fd
using SamParserThreadpool = htsThreadPool; // SAM_hdr
using SamHeader = bam_hdr_t;
using SamFlag = uint16_t;
using SamRecord = bam1_t;

namespace combinelab {
  namespace samutils {
    enum CIGAROp {
      OP_UNKNOWN=-1,
      OP_CMATCH=0,
      OP_CINS=1,
      OP_CDEL=2,
      OP_CREF_SKIP=3,
      OP_CSOFT_CLIP=4,
      OP_CHARD_CLIP=5,
      OP_CPAD=6,
      OP_CBASE_MATCH=7,
      OP_CBASE_MISMATCH=8//,
      //OP_CBASE_BACK=9
    };

    inline CIGAROp get_cigar_op(unsigned int opIn) {
      CIGAROp op = CIGAROp::OP_UNKNOWN;
      switch (opIn) {
      case 0:
        return CIGAROp::OP_CMATCH;
      case 1:
        return CIGAROp::OP_CINS;
      case 2:
        return CIGAROp::OP_CDEL;
      case 3:
        return CIGAROp::OP_CREF_SKIP;
      case 4:
        return CIGAROp::OP_CSOFT_CLIP;
      case 5:
        return CIGAROp::OP_CHARD_CLIP;
      case 6:
        return CIGAROp::OP_CPAD;
      case 7:
        return CIGAROp::OP_CBASE_MATCH;
      case 8:
        return CIGAROp::OP_CBASE_MISMATCH;
      default:
        return op;
      }
    }

    inline void closeOrDie(SamFile* fp, std::shared_ptr<spdlog::logger> logger) {
      int r = sam_close(fp);
      if ( r < 0 ) {
        logger->error("Error closing SAM/BAM file. Salmon will now exit. "
                      "Please report this bug on GitHub\n");
        std::exit(1);
      }
    }

    inline SamRecord* bam_dup(const SamRecord* r) { return ::bam_dup1(r); }
    inline void bam_destroy(SamRecord* r) { ::bam_destroy1(r); }
    inline SamRecord* bam_init() { return bam_init1(); }

    inline const constexpr char* header_get_ref_name(SamHeader* hdr, int32_t i) {
      return hdr->target_name[i];
    }

    inline const constexpr uint32_t header_get_ref_len(SamHeader* hdr, int32_t i) {
      return hdr->target_len[i];
    }

    inline const constexpr char* header_get_ref_name(SamHeader* hdr, uint32_t i) {
      return hdr->target_name[i];
    }

    inline const constexpr uint32_t header_get_ref_len(SamHeader* hdr, uint32_t i) {
      return hdr->target_len[i];
    }

    inline const constexpr bool bam_consume_seq(CIGAROp op) {
      return (0x193>>(op))&1;
    }

    inline const constexpr bool bam_consume_ref(CIGAROp op) {
      return (0x18d>>(op))&1;
    }

    inline const constexpr SamFlag bam_flag(SamRecord* bam) {
      return ((bam)->core.flag);
    }

    inline const constexpr int32_t bam_ref(SamRecord* bam) {
      return ((bam)->core.tid);
    }

    inline const constexpr int32_t bam_mate_ref(SamRecord* bam) {
      return ((bam)->core.mtid);
    }

    inline const constexpr uint8_t bam_name_len(SamRecord* bam) {
      return ((bam)->core.l_qname);
    }

    inline constexpr char* bam_name(SamRecord* bam) {
      return bam_get_qname(bam);
    }

    inline constexpr uint8_t* bam_seq(SamRecord* bam) {
      return bam_get_seq(bam);
    }

    inline constexpr auto bam_seq_len(SamRecord* bam) -> decltype((bam)->core.l_qseq) {
      return ((bam)->core.l_qseq);
    }


    inline constexpr uint8_t* bam_qual(SamRecord* bam) {
      return bam_get_qual(bam);
    }

    inline constexpr auto bam_qual_len(SamRecord* bam) -> decltype(bam_seq_len(bam)) {
      return bam_seq_len(bam);
    }


    inline constexpr uint32_t* bam_cigar(SamRecord* bam) {
      return bam_get_cigar(bam);
    }

    inline const constexpr uint32_t bam_cigar_len(SamRecord* bam) {
      return ((bam)->core.n_cigar);
    }

    inline constexpr auto BAMSeqI(SamRecord* bam, int32_t i) -> decltype(bam_seqi(bam_get_seq(bam), i)) {
      return bam_seqi(bam_get_seq(bam), i);
    }

    inline constexpr auto BAMSeqI(uint8_t* bamseq, int32_t i) -> decltype(bam_seqi(bamseq, i)) {
      return bam_seqi(bamseq, i);
    }

    inline const constexpr bool bam_strand(SamRecord* bam) {
      return bam_is_rev(bam);
    }

    inline const constexpr int32_t bam_pos(SamRecord* bam) {
      return ((bam)->core.pos);
    }

  }
}

#endif // __SAMTYPES_HPP__
