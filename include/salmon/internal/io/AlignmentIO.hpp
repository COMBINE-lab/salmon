#ifndef SALMON_INTERNAL_IO_ALIGNMENT_IO_HPP
#define SALMON_INTERNAL_IO_ALIGNMENT_IO_HPP

#include <cstdint>
#include <cstring>

extern "C" {
#include <htslib/hts.h>
#include <htslib/sam.h>
#undef max
#undef min
}

namespace salmon::io {

struct AlignmentHeaderTag {
  char* str{nullptr};
  int len{0};
  AlignmentHeaderTag* next{nullptr};
};

struct AlignmentHeaderRef {
  char* name{nullptr};
  uint32_t len{0};
  AlignmentHeaderTag* tag{nullptr};
};

struct AlignmentHeader {
  sam_hdr_t* raw{nullptr};
  int32_t nref{0};
  AlignmentHeaderRef* ref{nullptr};
  int32_t ref_count{1};
};

struct AlignmentFileHandle {
  samFile* raw{nullptr};
  AlignmentHeader* header{nullptr};
};

AlignmentHeader* wrapHeader(sam_hdr_t* raw);
void destroyHeader(AlignmentHeader* header);

} // namespace salmon::io

using bam_seq_t = bam1_t;
using AlignmentFileHandle = salmon::io::AlignmentFileHandle;
using AlignmentHeader = salmon::io::AlignmentHeader;
using AlignmentHeaderRef = salmon::io::AlignmentHeaderRef;
using AlignmentHeaderTag = salmon::io::AlignmentHeaderTag;
using cigar_op = int;

#ifndef BAM_CONSUME_SEQ
#define BAM_CONSUME_SEQ(op) (bam_cigar_type((op)) & 1)
#endif

#ifndef BAM_CONSUME_REF
#define BAM_CONSUME_REF(op) (bam_cigar_type((op)) & 2)
#endif

#ifndef BAM_CBASE_MATCH
#define BAM_CBASE_MATCH BAM_CEQUAL
#endif

#ifndef BAM_CBASE_MISMATCH
#define BAM_CBASE_MISMATCH BAM_CDIFF
#endif

constexpr int BAM_UNKNOWN = -1;
inline bam_seq_t* bam_dup(bam_seq_t* record) { return bam_dup1(record); }
namespace salmon::io {
inline bam_seq_t* bam_init() { return bam_init1(); }
inline void bam_destroy(bam_seq_t* record) {
  if (record != nullptr) {
    bam_destroy1(record);
  }
}
} // namespace salmon::io

inline bam_seq_t* bam_init() { return salmon::io::bam_init(); }
inline void bam_destroy(bam_seq_t* record) { salmon::io::bam_destroy(record); }
inline char* bam_name(bam_seq_t* record) { return bam_get_qname(record); }
inline uint32_t bam_name_len(bam_seq_t* record) {
  return static_cast<uint32_t>(std::strlen(bam_get_qname(record)));
}
inline int32_t bam_seq_len(bam_seq_t* record) { return record->core.l_qseq; }
inline int32_t bam_pos(bam_seq_t* record) { return record->core.pos; }
inline int32_t bam_ref(bam_seq_t* record) { return record->core.tid; }
inline int32_t bam_mate_ref(bam_seq_t* record) { return record->core.mtid; }
inline uint16_t bam_flag(bam_seq_t* record) { return record->core.flag; }
inline void bam_set_flag(bam_seq_t* record, uint16_t flag) {
  record->core.flag |= flag;
}
inline int32_t bam_map_qual(bam_seq_t* record) { return record->core.qual; }
inline int bam_strand(bam_seq_t* record) { return bam_is_rev(record); }
inline uint8_t* bam_seq(bam_seq_t* record) { return bam_get_seq(record); }
inline uint8_t* bam_qual(bam_seq_t* record) { return bam_get_qual(record); }
inline uint32_t* bam_cigar(bam_seq_t* record) { return bam_get_cigar(record); }
inline uint32_t bam_cigar_len(bam_seq_t* record) { return record->core.n_cigar; }
inline uint8_t* bam_aux_find(bam_seq_t* record, const char tag[2]) {
  return bam_aux_get(record, tag);
}
inline int32_t bam_aux_i(uint8_t* tag_value) { return bam_aux2i(tag_value); }

AlignmentFileHandle* openAlignmentFile(const char* path, const char* mode);
int closeAlignmentFile(AlignmentFileHandle* file);
inline AlignmentHeader* getAlignmentHeader(AlignmentFileHandle* file) {
  return (file == nullptr) ? nullptr : file->header;
}

inline int setAlignmentThreads(AlignmentFileHandle* file, int threads) {
  if (file == nullptr || file->raw == nullptr) {
    return -1;
  }
  return hts_set_threads(file->raw, threads);
}

inline int readAlignmentRecord(AlignmentFileHandle* file, bam_seq_t** read) {
  if (file == nullptr || file->raw == nullptr || file->header == nullptr) {
    return -1;
  }
  if (*read == nullptr) {
    *read = bam_init1();
  }
  return sam_read1(file->raw, file->header->raw, *read);
}

inline int writeAlignmentRecord(AlignmentFileHandle* file, bam_seq_t* read) {
  if (file == nullptr || file->raw == nullptr || file->header == nullptr) {
    return -1;
  }
  return (sam_write1(file->raw, file->header->raw, read) >= 0) ? 0 : -1;
}

inline void setAlignmentHeader(AlignmentFileHandle* file, AlignmentHeader* header) {
  if (file != nullptr) {
    file->header = header;
  }
}

inline int writeAlignmentHeader(AlignmentFileHandle* file) {
  if (file == nullptr || file->raw == nullptr || file->header == nullptr) {
    return -1;
  }
  return sam_hdr_write(file->raw, file->header->raw);
}

inline void alignmentHeaderIncrRef(AlignmentHeader* header) {
  if (header != nullptr) {
    ++header->ref_count;
  }
}

inline void alignmentHeaderDecrRef(AlignmentHeader* header) {
  if (header == nullptr) {
    return;
  }
  --header->ref_count;
  if (header->ref_count <= 0) {
    salmon::io::destroyHeader(header);
  }
}

#endif
