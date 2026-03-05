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

struct SAM_hdr_tag {
  char* str{nullptr};
  int len{0};
  SAM_hdr_tag* next{nullptr};
};

struct SAM_hdr_ref {
  char* name{nullptr};
  uint32_t len{0};
  SAM_hdr_tag* tag{nullptr};
};

struct SAM_hdr {
  sam_hdr_t* raw{nullptr};
  int32_t nref{0};
  SAM_hdr_ref* ref{nullptr};
  int32_t ref_count{1};
};

struct scram_fd {
  samFile* raw{nullptr};
  SAM_hdr* header{nullptr};
};

SAM_hdr* wrap_header(sam_hdr_t* raw);
void destroy_header(SAM_hdr* header);

} // namespace salmon::io

using bam_seq_t = bam1_t;
using scram_fd = salmon::io::scram_fd;
using SAM_hdr = salmon::io::SAM_hdr;
using SAM_hdr_ref = salmon::io::SAM_hdr_ref;
using SAM_hdr_tag = salmon::io::SAM_hdr_tag;
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

scram_fd* scram_open(const char* path, const char* mode);
int scram_close(scram_fd* file);
SAM_hdr* scram_get_header(scram_fd* file);
int scram_set_option(scram_fd* file, int option, int value);
int scram_get_seq(scram_fd* file, bam_seq_t** read);
int scram_put_seq(scram_fd* file, bam_seq_t* read);
void scram_set_header(scram_fd* file, SAM_hdr* header);
int scram_write_header(scram_fd* file);

inline void sam_hdr_incr_ref(SAM_hdr* header) {
  if (header != nullptr) {
    ++header->ref_count;
  }
}

inline void sam_hdr_decr_ref(SAM_hdr* header) {
  if (header == nullptr) {
    return;
  }
  --header->ref_count;
  if (header->ref_count <= 0) {
    salmon::io::destroy_header(header);
  }
}

#endif
