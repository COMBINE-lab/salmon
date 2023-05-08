/*
 * Copyright (c) 2013 Genome Research Ltd.
 * Author(s): James Bonfield
 * 
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 * 
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 * 
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *    Institute nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific
 *    prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _CRAM_STRUCTS_H_
#define _CRAM_STRUCTS_H_

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Defines in-memory structs for the basic file-format objects in the
 * CRAM format.
 *
 * The basic file format is:
 *     File-def SAM-hdr Container Container ...
 *
 * Container:
 *     Service-block data-block data-block ...
 *
 * Multiple blocks in a container are grouped together as slices,
 * also sometimes referred to as landmarks in the spec.
 */


#include <stdint.h>

#include "io_lib/hash_table.h"       // From io_lib aka staden-read
#include "io_lib/thread_pool.h"
#include "io_lib/mFILE.h"

#ifdef SAMTOOLS
// From within samtools/HTSlib
#  include "io_lib/string_alloc.h"
#else
// From within io_lib
#  include "io_lib/bam.h"              // For BAM header parsing
#endif

#define SEQS_PER_SLICE 10000
#define BASES_PER_SLICE (SEQS_PER_SLICE*500)
#define SLICE_PER_CNT  1

#define CRAM_SUBST_MATRIX "CGTNAGTNACTNACGNACGT"

// TN only in Cram v1
//#define TN_external

#define MAX_STAT_VAL 1024
//#define MAX_STAT_VAL 16
typedef struct {
    int freqs[MAX_STAT_VAL];
    HashTable *h;
    int nsamp; // total number of values added
    int nvals; // total number of unique values added
} cram_stats;

/* NB: matches java impl, not the spec */
enum cram_encoding {
    E_NULL               = 0,
    E_EXTERNAL           = 1,
    E_GOLOMB             = 2,
    E_HUFFMAN            = 3,
    E_BYTE_ARRAY_LEN     = 4,
    E_BYTE_ARRAY_STOP    = 5,
    E_BETA               = 6,
    E_SUBEXP             = 7,
    E_GOLOMB_RICE        = 8,
    E_GAMMA              = 9,
    E_NUM_CODECS         = 10, /* Number of codecs, not a real one. */
};

enum cram_external_type {
    E_INT                = 1,
    E_LONG               = 2,
    E_BYTE               = 3,
    E_BYTE_ARRAY         = 4,
    E_BYTE_ARRAY_BLOCK   = 5,
};

/* External IDs used by this implementation (only assumed during writing) */
enum cram_DS_ID {
    DS_CORE   = 0,
    DS_aux    = 1, // aux_blk
    DS_aux_OQ = 2,
    DS_aux_BQ = 3,
    DS_aux_BD = 4,
    DS_aux_BI = 5,
    DS_aux_FZ = 6, // also ZM:B
    DS_aux_oq = 7, // other qualities
    DS_aux_os = 8, // other sequences
    DS_aux_oz = 9, // other strings
    DS_ref,
    DS_RN, // name_blk
    DS_QS, // qual_blk
    DS_IN, // base_blk
    DS_SC, // soft_blk

    DS_BF, // start loop
    DS_CF,
    DS_AP,
    DS_RG,
    DS_MQ,
    DS_NS,
    DS_MF,
    DS_TS,
    DS_NP,
    DS_NF,
    DS_RL,
    DS_FN,
    DS_FC,
    DS_FP,
    DS_DL,
    DS_BA,
    DS_BS,
    DS_TL,
    DS_RI,
    DS_RS,
    DS_PD,
    DS_HC,
    DS_BB,
    DS_QQ,

    DS_TN, // end loop

    DS_RN_len,
    DS_SC_len,
    DS_BB_len,
    DS_QQ_len,

    DS_TC, // CRAM v1.0 tags
    DS_TM, // test
    DS_TV, // test
    
    DS_END,
};

/* "File Definition Structure" */
typedef struct {
    char    magic[4];
    uint8_t major_version;
    uint8_t minor_version;
    char    file_id[20];      // Filename or SHA1 checksum
} cram_file_def;

#define CRAM_MAJOR_VERS(v) ((v) >> 8)
#define CRAM_MINOR_VERS(v) ((v) & 0xff)
#define IS_CRAM_1_VERS(fd) (CRAM_MAJOR_VERS((fd)->version)==1)
#define IS_CRAM_2_VERS(fd) (CRAM_MAJOR_VERS((fd)->version)==2)
#define IS_CRAM_3_VERS(fd) (CRAM_MAJOR_VERS((fd)->version)==3)

struct cram_slice;

enum cram_block_method {
    BM_ERROR  = -1,
    RAW    = 0,
    GZIP   = 1,    // Z_FILTERED
    BZIP2  = 2,
    LZMA   = 3,
    RANS0  = 4,
    RANS1  = 10,   // Not externalised; stored as RANS (generic)
    GZIP_RLE = 11, // Z_RLE, NB: not externalised in CRAM
    GZIP_1 = 12,   // Z_DEFAULT_STRATEGY level 1, NB: not externalised in CRAM
};

enum cram_content_type {
    CT_ERROR           = -1,
    FILE_HEADER        = 0,
    COMPRESSION_HEADER = 1,
    MAPPED_SLICE       = 2,
    UNMAPPED_SLICE     = 3, // CRAM_1_VERS only
    EXTERNAL           = 4,
    CORE               = 5,
};

/* Compression metrics */
typedef struct {
    // number of trials and time to next trial
    int trial;
    int next_trial;

    // aggregate sizes during trials
    int sz_gz_rle;
    int sz_gz_def;
    int sz_gz_1;
    int sz_rans0;
    int sz_rans1;
    int sz_bzip2;
    int sz_lzma;

    // resultant method from trials
    int method;
    int strat;

    // Revisions of method, to allow culling of continually failing ones.
    int gz_rle_cnt;
    int gz_def_cnt;
    int gz_1_cnt;
    int rans0_cnt;
    int rans1_cnt;
    int bzip2_cnt;
    int lzma_cnt;
    int revised_method;

    double gz_rle_extra;
    double gz_def_extra;
    double gz_1_extra;
    double rans0_extra;
    double rans1_extra;
    double bzip2_extra;
    double lzma_extra;
} cram_metrics;

/* Block */
typedef struct {
    enum cram_block_method  method, orig_method;
    enum cram_content_type  content_type;
    int32_t  content_id;
    int32_t  comp_size;
    int32_t  uncomp_size;
    uint32_t crc32;
    int32_t  idx; /* offset into data */
    unsigned char    *data;

    // For bit I/O
    size_t alloc;
    size_t byte;
    int bit;

    // To aid compression
    cram_metrics *m; // used to track aux block compression only
} cram_block;

struct cram_codec; /* defined in cram_codecs.h */
struct cram_map;

#define CRAM_MAP_HASH 32
#define CRAM_MAP(a,b) (((a)*3+(b))&(CRAM_MAP_HASH-1))

/* Compression header block */
typedef struct {
    int32_t ref_seq_id;
    int32_t ref_seq_start;
    int32_t ref_seq_span;
    int32_t num_records;
    int32_t num_landmarks;
    int32_t *landmark;

    /* Flags from preservation map */
    int mapped_qs_included;
    int unmapped_qs_included;
    int unmapped_placed;
    int qs_included;
    int read_names_included;
    int AP_delta;
    // indexed by ref-base and subst. code
    char substitution_matrix[5][4];

    // TD Dictionary as a concatenated block
    cram_block *TD_blk;  // Tag Dictionary
    int nTL;		 // number of TL entries in TD
    unsigned char **TL;  // array of size nTL, pointer into TD_blk.
    HashTable *TD;       // for encoding, keyed on TD entries
    
    HashTable *preservation_map;
    struct cram_map *rec_encoding_map[CRAM_MAP_HASH];
    struct cram_map *tag_encoding_map[CRAM_MAP_HASH];

    struct cram_codec *codecs[DS_END];

    char *uncomp; // A single block of uncompressed data
    size_t uncomp_size, uncomp_alloc;

    unsigned int data_series; // See cram_fields enum below
} cram_block_compression_hdr;

typedef struct cram_map {
    int key;    /* 0xe0 + 3 bytes */
    enum cram_encoding encoding;
    int offset; /* Offset into a single block of memory */
    int size;   /* Size */
    struct cram_codec *codec;
    struct cram_map *next; // for noddy internal hash
} cram_map;

typedef struct {
    struct cram_codec *codec;
    cram_block *blk;
    cram_metrics *m;
} cram_tag_map;

/* Mapped or unmapped slice header block */
typedef struct {
    enum cram_content_type content_type;
    int32_t ref_seq_id;     /* if content_type == MAPPED_SLICE */
    int32_t ref_seq_start;  /* if content_type == MAPPED_SLICE */
    int32_t ref_seq_span;   /* if content_type == MAPPED_SLICE */
    int32_t num_records;
    int64_t record_counter;
    int32_t num_blocks;
    int32_t num_content_ids;
    int32_t *block_content_ids;
    int32_t ref_base_id;    /* if content_type == MAPPED_SLICE */
    unsigned char md5[16];
    HashTable *tags;        /* hash of optional tags */
    uint32_t BD_crc;        /* base call digest */
    uint32_t SD_crc;        /* quality score digest */
} cram_block_slice_hdr;

struct ref_entry;

/*
 * Container.
 *
 * Conceptually a container is split into slices, and slices into blocks.
 * However on disk it's just a list of blocks and we need to query the
 * block types to identify the start/end points of the slices.
 *
 * OR... are landmarks the start/end points of slices?
 */
typedef struct {
    int32_t  length;
    int32_t  ref_seq_id;
    int32_t  ref_seq_start;
    int32_t  ref_seq_span;
    int64_t  record_counter;
    int64_t  num_bases;
    int32_t  num_records;
    int32_t  num_blocks;
    int32_t  num_landmarks;
    int32_t *landmark;

    /* Size of container header above */
    size_t   offset;
    
    /* Compression header is always the first block? */
    cram_block_compression_hdr *comp_hdr;
    cram_block *comp_hdr_block;

    /* For construction purposes */
    int max_slice, curr_slice;   // maximum number of slices
    int max_rec, curr_rec;       // current and max recs per slice
    int max_c_rec, curr_c_rec;   // current and max recs per container
    int slice_rec;               // rec no. for start of this slice
    int curr_ref;                // current ref ID. -2 for no previous
    int last_pos;                // last record position
    struct cram_slice **slices, *slice;
    int pos_sorted;              // boolean, 1=>position sorted data
    int max_apos;                // maximum position, used if pos_sorted==0
    int last_slice;              // number of reads in last slice (0 for 1st)
    int multi_seq;               // true if packing multi seqs per cont/slice
    int unsorted;		 // true is AP_delta is 0.

    /* Copied from fd before encoding, to allow multi-threading */
    int ref_start, first_base, last_base, ref_id, ref_end;
    char *ref;
    //struct ref_entry *ref;

    /* For multi-threading */
    bam_seq_t **bams;

    /* Statistics for encoding */
    cram_stats *stats[DS_END];

    HashTable *tags_used; // cram_tag_map[], per tag types in use.

    int *refs_used;       // array of frequency of ref seq IDs

    // For experimental name delta
    char *last_name;

    uint32_t crc32;       // Raw container bytes CRC

    uint64_t s_num_bases; // number of bases in this slice
} cram_container;

/*
 * A single cram record
 */
typedef struct {
    struct cram_slice *s; // Filled out by cram_decode only

    int32_t ref_id;       // fixed for all recs in slice?
    int32_t flags;        // BF
    int32_t cram_flags;   // CF
    int32_t len;          // RL
    int32_t apos;         // AP
    int32_t rg;           // RG
    int32_t name;         // RN; idx to s->names_blk
    int32_t name_len;
    int32_t mate_line;    // index to another cram_record
    int32_t mate_ref_id;
    int32_t mate_pos;     // NP
    int32_t tlen;         // TS

    // Auxiliary data
    int32_t ntags;        // TC
    int32_t aux;          // idx to s->aux_blk
    int32_t aux_size;     // total size of packed ntags in aux_blk
#ifndef TN_external
    int32_t TN_idx;       // TN; idx to s->TN;
#else
    int32_t tn;           // idx to s->tn_blk
#endif
    int     TL;

    int32_t seq;          // idx to s->seqs_blk
    int32_t qual;         // idx to s->qual_blk
    int32_t cigar;        // idx to s->cigar
    int32_t ncigar;
    int32_t aend;         // alignment end
    int32_t mqual;        // MQ

    int32_t feature;      // idx to s->feature
    int32_t nfeature;     // number of features
    int32_t mate_flags;   // MF
} cram_record;

// Accessor macros as an analogue of the bam ones
#define cram_qname(c)    (&(c)->s->name_blk->data[(c)->name])
#define cram_seq(c)      (&(c)->s->seqs_blk->data[(c)->seq])
#define cram_qual(c)     (&(c)->s->qual_blk->data[(c)->qual])
#define cram_aux(c)      (&(c)->s->aux_blk->data[(c)->aux])
#define cram_seqi(c,i)   (cram_seq((c))[(i)])
#define cram_name_len(c) ((c)->name_len)
#define cram_strand(c)   (((c)->flags & BAM_FREVERSE) != 0)
#define cram_mstrand(c)  (((c)->flags & BAM_FMREVERSE) != 0)
#define cram_cigar(c)    (&((cr)->s->cigar)[(c)->cigar])

/*
 * A feature is a base difference, used for the sequence reference encoding.
 * (We generate these internally when writing CRAM.)
 */
typedef struct {
    union {
	struct {
	    int pos;
	    int code;
	    int base;    // substitution code
	} X;
	struct {
	    int pos;
	    int code;
	    int base;    // actual base & qual
	    int qual;
	} B;
	struct {
	    int pos;
	    int code;
	    int seq_idx; // index to s->seqs_blk
	    int len;
	} b;
	struct {
	    int pos;
	    int code;
	    int qual;
	} Q;
	struct {
	    int pos;
	    int code;
	    int len;
	    int seq_idx; // soft-clip multiple bases
	} S;
	struct {
	    int pos;
	    int code;
	    int len;
	    int seq_idx; // insertion multiple bases
	} I;
	struct {
	    int pos;
	    int code;
	    int base; // insertion single base
	} i;
	struct {
	    int pos;
	    int code;
	    int len;
	} D;
	struct {
	    int pos;
	    int code;
	    int len;
	} N;
	struct {
	    int pos;
	    int code;
	    int len;
	} P;
	struct {
	    int pos;
	    int code;
	    int len;
	} H;
    };
} cram_feature;

//// Turns [A-Z][A-Z] into an integer from 0 to 32*32
//#define ID(a) ((((a)[0]-'A')<<5)+(a)[1]-'A')

/*
 * A slice is really just a set of blocks, but it
 * is the logical unit for decoding a number of
 * sequences.
 */
typedef struct cram_slice {
    cram_block_slice_hdr *hdr;
    cram_block *hdr_block;
    cram_block **block;
    cram_block **block_by_id;

    /* State used during encoding/decoding */
    int last_apos, max_apos;

    /* Array of decoded cram records */
    cram_record *crecs;

    /* An dynamically growing buffers for data pointed
     * to by crecs[] array.
     */
    uint32_t  *cigar;
    uint32_t   cigar_alloc;
    uint32_t   ncigar;

    cram_feature *features;
    int           nfeatures;
    int           afeatures; // allocated size of features

#ifndef TN_external
    // TN field (Tag Name)
    uint32_t      *TN;
    int           nTN, aTN;  // used and allocated size for TN[]
#else
    cram_block *tn_blk;
    int tn_id;
#endif

    // For variable sized elements which are always external blocks.
    cram_block *name_blk;
    cram_block *seqs_blk;
    cram_block *qual_blk;
    cram_block *base_blk;
    cram_block *soft_blk;
    cram_block *aux_blk;  // BAM aux block, used when going from CRAM to BAM

    HashTable *pair[2];      // for identifying read-pairs in this slice.

    char *ref;               // slice of current reference
    int ref_start;           // start position of current reference;
    int ref_end;             // end position of current reference;
    int ref_id;

    uint32_t BD_crc;         // base call digest
    uint32_t SD_crc;         // quality score digest

    // For going from BAM to CRAM; an array of auxiliary blocks per type
    int naux_block;
    cram_block **aux_block;
} cram_slice;

/*-----------------------------------------------------------------------------
 * Consider moving reference handling to cram_refs.[ch]
 */
// from fa.fai / samtools faidx files
typedef struct ref_entry {
    char *name;
    char *fn;
    int64_t length;
    int64_t offset;
    int bases_per_line;
    int line_length;
    int64_t count;	   // for shared references so we know to dealloc seq
    char *seq;
    mFILE *mf;
} ref_entry;

// References structure.
typedef struct {
    string_alloc_t *pool;  // String pool for holding filenames and SN vals

    HashTable *h_meta;     // ref_entry*, index by name
    ref_entry **ref_id;    // ref_entry*, index by ID
    int nref;              // number of ref_entry

    char *fn;              // current file opened
    FILE *fp;              // and the FILE* to go with it.

    int count;             // how many cram_fd sharing this refs struct

    pthread_mutex_t lock;  // Mutex for multi-threaded updating
    ref_entry *last;       // Last queried sequence
    int last_id;           // Used in cram_ref_decr_locked to delay free
} refs_t;

/*-----------------------------------------------------------------------------
 * CRAM index
 *
 * Detect format by number of entries per line.
 * 5 => 1.0 (refid, start, nseq, C offset, slice)
 * 6 => 1.1 (refid, start, span, C offset, S offset, S size)
 *
 * Indices are stored in a nested containment list, which is trivial to set
 * up as the indices are on sorted data so we're appending to the nclist
 * in sorted order. Basically if a slice entirely fits within a previous
 * slice then we append to that slices list. This is done recursively.
 *
 * Lists are sorted on two dimensions: ref id + slice coords.
 */
typedef struct cram_index {
    int nslice, nalloc;   // total number of slices
    struct cram_index *e; // array of size nslice

    int     refid;  // 1.0                 1.1
    int     start;  // 1.0                 1.1
    int     end;    //                     1.1
    int     nseq;   // 1.0 - undocumented
    int     slice;  // 1.0 landmark index, 1.1 landmark value
    int     len;    //                     1.1 - size of slice in bytes
    int64_t offset; // 1.0                 1.1
} cram_index;

typedef struct {
    int refid;
    int start;
    int end;
} cram_range;

/*-----------------------------------------------------------------------------
 */
/* CRAM File handle */

typedef struct spare_bams {
    bam_seq_t **bams;
    struct spare_bams *next;
} spare_bams;

#if defined(CRAM_IO_CUSTOM_BUFFERING)
typedef size_t (*cram_io_C_FILE_fread_t)(void *ptr, size_t size, size_t nmemb, void *stream);
typedef size_t (*cram_io_C_FILE_fwrite_t)(void *ptr, size_t size, size_t nmemb, void *stream);
typedef int    (*cram_io_C_FILE_fseek_t)(void * fd, off_t offset, int whence);
typedef off_t  (*cram_io_C_FILE_ftell_t)(void * fd);

typedef struct {
    void                   *user_data;
    cram_io_C_FILE_fread_t  fread_callback;
    cram_io_C_FILE_fseek_t  fseek_callback;
    cram_io_C_FILE_ftell_t  ftell_callback;
} cram_io_input_t;

typedef struct {
    void                   *user_data;
    cram_io_C_FILE_fwrite_t fwrite_callback;
    cram_io_C_FILE_ftell_t  ftell_callback;
} cram_io_output_t;

typedef cram_io_input_t * (*cram_io_allocate_read_input_t)(char const * filename, int const decompress);
typedef cram_io_input_t * (*cram_io_deallocate_read_input_t)(cram_io_input_t * obj);

typedef cram_io_output_t * (*cram_io_allocate_write_output_t)(char const * filename);
typedef cram_io_output_t * (*cram_io_deallocate_write_output_t)(cram_io_output_t * obj);



// FIXME: make cram_fd_input_buffer and cram_fd_input_buffer the same thing.
// Ie cram_fd_io_buffer and internals fp_io_*.

typedef struct {
    /* input buffer size */
    size_t         fp_in_buf_size;
    /* input buffer base pointer */
    char          *fp_in_buffer;
    /* position of buffer start in file */
    uint64_t       fp_in_buf_start;
    /* start of window pointer; same as fp_in_buffer */
    char          *fp_in_buf_pa;
    /* window current pointer */
    char          *fp_in_buf_pc;
    /* window end pointer;  same as fp_in_buffer + fp_in_buf_size (no seeks) */
    char          *fp_in_buf_pe;    
} cram_fd_input_buffer;

typedef struct {
    /* output buffer size */
    size_t         fp_out_buf_size;
    /* output buffer base pointer */
    char          *fp_out_buffer;
    /* position of buffer start in file */
    uint64_t       fp_out_buf_start;
    /* start of window pointer; same as fp_out_buffer */
    char          *fp_out_buf_pa;
    /* window current pointer */
    char          *fp_out_buf_pc;
    /* window end pointer */
    char          *fp_out_buf_pe;    
} cram_fd_output_buffer;
#endif

typedef struct {
    FILE                 *fp_in;
#if defined(CRAM_IO_CUSTOM_BUFFERING)
    cram_fd_input_buffer            *fp_in_buffer;
    cram_io_input_t                 *fp_in_callbacks;
    cram_io_allocate_read_input_t    fp_in_callback_allocate_function;
    cram_io_deallocate_read_input_t  fp_in_callback_deallocate_function;

    cram_fd_output_buffer            *fp_out_buffer;
    cram_io_output_t                 *fp_out_callbacks;
    cram_io_allocate_write_output_t   fp_out_callback_allocate_function;
    cram_io_deallocate_write_output_t fp_out_callback_deallocate_function;
#endif
    
    FILE          *fp_out;
    int            mode;     // 'r' or 'w'
    int            version;
    cram_file_def *file_def;
    SAM_hdr       *header;

    char          *prefix;
    int64_t        record_counter;
    int            err;

    // Most recent compression header decoded
    //cram_block_compression_hdr *comp_hdr;
    //cram_block_slice_hdr       *slice_hdr;

    // Current container being processed.
    cram_container *ctr;

    // positions for encoding or decoding
    int first_base, last_base;

    // cached reference portion
    refs_t *refs;              // ref meta-data structure
    char *ref, *ref_free;      // current portion held in memory
    int   ref_id;
    int   ref_start;
    int   ref_end;
    char *ref_fn;   // reference fasta filename

    // compression level and metrics
    int level;
    cram_metrics *m[DS_END];
    HashTable *tags_used; // cram_metrics[], per tag types in use.

    // options
    int decode_md; // Whether to export MD and NM tags
    int verbose;
    int seqs_per_slice;
    int bases_per_slice;
    int slices_per_container;
    int embed_ref;
    int no_ref;
    int ignore_md5;
    int use_bz2;
    int use_rans;
    int use_lzma;
    int shared_ref;
    enum quality_binning binning;
    unsigned int required_fields;
    cram_range range;

    // lookup tables, stored here so we can be trivially multi-threaded
    unsigned int bam_flag_swap[0x1000]; // cram -> bam flags
    unsigned int cram_flag_swap[0x1000];// bam -> cram flags
    unsigned char L1[256];              // ACGT{*} ->0123{4}
    unsigned char L2[256];              // ACGTN{*}->01234{5}
    char cram_sub_matrix[32][32];	// base substituion codes

    int         index_sz;
    cram_index *index;                  // array, sizeof index_sz
    off_t first_container;
    int eof;
    int last_slice;                     // number of recs encoded in last slice
    int multi_seq;
    int unsorted;
    int empty_container; 		// Marker for EOF block
    
    // thread pool
    int own_pool;
    t_pool *pool;
    t_results_queue *rqueue;
    pthread_mutex_t *metrics_lock;
    pthread_mutex_t *ref_lock;
    spare_bams *bl;
    pthread_mutex_t *bam_list_lock;
    void *job_pending;

    int ooc;                            // out of containers.
    int ignore_chksum;
    int lossy_read_names;
    int preserve_aux_order;             // if set implies emitting RG, MD and NM
    int preserve_aux_size;              // does not replace 'i' with 'c' etc in aux.
} cram_fd;

#if defined(CRAM_IO_CUSTOM_BUFFERING)
extern size_t cram_io_input_buffer_read(void *ptr, size_t size, size_t nmemb, cram_fd * fd);
extern int cram_io_input_buffer_seek(cram_fd * fd, off_t offset, int whence);
extern int cram_io_input_buffer_underflow(cram_fd * fd);
extern char * cram_io_input_buffer_fgets(char * s, int size, cram_fd * fd);
extern int cram_io_flush_output_buffer(cram_fd *fd);
#endif

#if defined(CRAM_IO_CUSTOM_BUFFERING)
#define CRAM_IO_GETC(fd) ((fd->fp_in_buffer->fp_in_buf_pc!=fd->fp_in_buffer->fp_in_buf_pe) ? ((int)((unsigned char)(*(fd->fp_in_buffer->fp_in_buf_pc++)))) : cram_io_input_buffer_underflow(fd))
#define CRAM_IO_READ(ptr, size, nmemb, fd) cram_io_input_buffer_read(ptr,size,nmemb,fd)
#define CRAM_IO_SEEK(fd, offset, whence) cram_io_input_buffer_seek(fd, offset, whence)
#define CRAM_IO_TELLO(fd) (fd->fp_in_buffer->fp_in_buf_start +(fd->fp_in_buffer->fp_in_buf_pc-fd->fp_in_buffer->fp_in_buf_pa))
#define CRAM_IO_FGETS(s,size,fd) cram_io_input_buffer_fgets(s,size,fd)
#define CRAM_IO_PUTC(c,fd) cram_io_output_buffer_putc(c,fd)
#define CRAM_IO_WRITE(ptr, size, nmemb, fd) cram_io_output_buffer_write(ptr,size,nmemb,fd)
#define CRAM_IO_FLUSH(fd) cram_io_flush_output_buffer((fd))

#else // ! CRAM_IO_CUSTOM_BUFFERING
#define CRAM_IO_GETC(fd) getc(fd->fp_in)
#define CRAM_IO_READ(ptr, size, nmemb, fd) fread(ptr,size,nmemb,fd->fp_in)
#define CRAM_IO_TELLO(fd) ftello(fd->fp_in)
#define CRAM_IO_SEEK(fd, offset, whence) fseeko(fd->fp_in,offset,whence)
#define CRAM_IO_FGETS(s,size,fd) fgets(s,size,fd->fp_in)
#define CRAM_IO_PUTC(c,fd) putc(c,fd->fp_out)
#define CRAM_IO_WRITE(ptr, size, nmemb, fd) fwrite(ptr,size,nmemb,fd->fp_out)
#define CRAM_IO_FLUSH(fd) (fd->fp_out ? fflush(fd->fp_out) : 0)

#endif // end CRAM_IO_CUSTOM_BUFFERING

// REQUIRED_FIELDS
enum sam_fields {
    SAM_QNAME = 0x00000001,
    SAM_FLAG  = 0x00000002,
    SAM_RNAME = 0x00000004,
    SAM_POS   = 0x00000008,
    SAM_MAPQ  = 0x00000010,
    SAM_CIGAR = 0x00000020,
    SAM_RNEXT = 0x00000040,
    SAM_PNEXT = 0x00000080,
    SAM_TLEN  = 0x00000100,
    SAM_SEQ   = 0x00000200,
    SAM_QUAL  = 0x00000400,
    SAM_AUX   = 0x00000800,
    SAM_RGAUX = 0x00001000,
};

// Translation of required fields to cram data series
enum cram_fields {
    CRAM_BF = 0x00000001,
    CRAM_AP = 0x00000002,
    CRAM_FP = 0x00000004,
    CRAM_RL = 0x00000008,
    CRAM_DL = 0x00000010,
    CRAM_NF = 0x00000020,
    CRAM_BA = 0x00000040,
    CRAM_QS = 0x00000080,
    CRAM_FC = 0x00000100,
    CRAM_FN = 0x00000200,
    CRAM_BS = 0x00000400,
    CRAM_IN = 0x00000800,
    CRAM_RG = 0x00001000,
    CRAM_MQ = 0x00002000,
    CRAM_TL = 0x00004000,
    CRAM_RN = 0x00008000,
    CRAM_NS = 0x00010000,
    CRAM_NP = 0x00020000,
    CRAM_TS = 0x00040000,
    CRAM_MF = 0x00080000,
    CRAM_CF = 0x00100000,
    CRAM_RI = 0x00200000,
    CRAM_RS = 0x00400000,
    CRAM_PD = 0x00800000,
    CRAM_HC = 0x01000000,
    CRAM_SC = 0x02000000,
    CRAM_BB = 0x04000000,
    CRAM_BB_len = 0x08000000,
    CRAM_QQ = 0x10000000,
    CRAM_QQ_len = 0x20000000,
    CRAM_aux= 0x40000000,
    CRAM_ALL= 0x7fffffff,
};

// A CIGAR opcode, but not necessarily the implications of it. Eg FC/FP may
// encode a base difference, but we don't need to know what it is for CIGAR.
// If we have a soft-clip or insertion, we do need SC/IN though to know how
// long that array is.
#define CRAM_CIGAR (CRAM_FN | CRAM_FP | CRAM_FC | CRAM_DL | CRAM_IN | \
		    CRAM_SC | CRAM_HC | CRAM_PD | CRAM_RS | CRAM_RL | CRAM_BF)

#define CRAM_SEQ (CRAM_CIGAR | CRAM_BA | CRAM_BS | \
		  CRAM_RL    | CRAM_AP | CRAM_BB)

#define CRAM_QUAL (CRAM_CIGAR | CRAM_RL | CRAM_AP | CRAM_QS | CRAM_QQ)

enum cram_option {
    CRAM_OPT_DECODE_MD,
    CRAM_OPT_PREFIX,
    CRAM_OPT_VERBOSITY,
    CRAM_OPT_SEQS_PER_SLICE,
    CRAM_OPT_SLICES_PER_CONTAINER,
    CRAM_OPT_RANGE,
    CRAM_OPT_VERSION,
    CRAM_OPT_EMBED_REF,
    CRAM_OPT_IGNORE_MD5,
    CRAM_OPT_REFERENCE,
    CRAM_OPT_MULTI_SEQ_PER_SLICE,
    CRAM_OPT_NO_REF,
    CRAM_OPT_USE_BZIP2,
    CRAM_OPT_SHARED_REF,
    CRAM_OPT_NTHREADS,
    CRAM_OPT_THREAD_POOL,
    CRAM_OPT_BINNING,
    CRAM_OPT_USE_ARITH,
    CRAM_OPT_USE_LZMA,
    CRAM_OPT_REQUIRED_FIELDS,
    CRAM_OPT_USE_RANS,
    CRAM_OPT_IGNORE_CHKSUM,
    CRAM_OPT_BASES_PER_SLICE,
    CRAM_OPT_LOSSY_READ_NAMES,
    CRAM_OPT_PRESERVE_AUX_ORDER,
    CRAM_OPT_PRESERVE_AUX_SIZE
};

/* BF bitfields */
/* Corrected in 1.1. Use bam_flag_swap[bf] and BAM_* macros for 1.0 & 1.1 */
#define CRAM_FPAIRED      256
#define CRAM_FPROPER_PAIR 128
#define CRAM_FUNMAP        64
#define CRAM_FREVERSE      32
#define CRAM_FREAD1        16
#define CRAM_FREAD2         8
#define CRAM_FSECONDARY     4
#define CRAM_FQCFAIL        2
#define CRAM_FDUP           1

#define DS_aux_S "\001"
#define DS_aux_OQ_S "\002"
#define DS_aux_BQ_S "\003"
#define DS_aux_BD_S "\004"
#define DS_aux_BI_S "\005"
#define DS_aux_FZ_S "\006"
#define DS_aux_oq_S "\007"
#define DS_aux_os_S "\010"
#define DS_aux_oz_S "\011"

#define CRAM_M_REVERSE  1
#define CRAM_M_UNMAP    2


/* CF bitfields */
#define CRAM_FLAG_PRESERVE_QUAL_SCORES (1<<0)
#define CRAM_FLAG_DETACHED             (1<<1)
#define CRAM_FLAG_MATE_DOWNSTREAM      (1<<2)
#define CRAM_FLAG_NO_SEQ               (1<<3)
#define CRAM_FLAG_MASK                 ((1<<4)-1)

/* Internal only */
#define CRAM_FLAG_STATS_ADDED  (1<<30)
#define CRAM_FLAG_DISCARD_NAME (1<<31)

#ifdef __cplusplus
}
#endif

#endif /* _CRAM_STRUCTS_H_ */
