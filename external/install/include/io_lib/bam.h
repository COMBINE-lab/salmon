/*
 * Copyright (c) 2013 Genome Research Ltd.
 * Author(s): James Bonfield, Rob Davies
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

/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2010-3
 */

/*! \file
 * The primary SAM/BAM API.
 *
 * Consider using scram.h if you wish to also have support for CRAM.
 */ 

#ifndef _BAM_H_
#define _BAM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
#include <zlib.h>

#include "io_lib/os.h"
#include "io_lib/hash_table.h"
#include "io_lib/sam_header.h"
#include "io_lib/thread_pool.h"
#include "io_lib/binning.h"

/* BAM header structs */
typedef struct tag_list {
    char *value; /* NULL => end of tags */
    int key;
    int length;  
} tag_list_t;

/* The main bam sequence struct */
typedef struct bam_seq_s {
    uint32_t alloc; /* total size of this struct + 'data' onwards */
    uint32_t blk_size;

    /* The raw bam block follows, in same order as on the disk */
    /* This is the analogue of a bam1_core_t in samtools */
    int32_t  ref;
    int32_t  pos;

    union {
	struct {
	    uint32_t name_len:8, map_qual:8, bin:16;
	};
	uint32_t bin_packed;
    };
    union {
	struct {
	    uint32_t cigar_len:16, flag:16;
	};
	uint32_t flag_packed;
    };

    int32_t  len;
    int32_t  mate_ref;
    int32_t  mate_pos;
    int32_t  ins_size;

    /* Followed by arbitrary bytes of packed data */
    unsigned char data; /* unknown size */
} bam_seq_t;


/* Auxillary field handling */
typedef union {
    char    *s;
    int      i;
    uint64_t i64;
    float    f;
    double   d;
    struct {
	int n, t;
	unsigned char *s;
    } B;
} bam_aux_t;

/* Struct for making arrays of aux tags */

typedef struct {
    char     tag[2];
    char     type;
    uint32_t array_len;
    union {
	char    *z;
	uint8_t *h;
	char     a;
	int32_t  i;
	uint32_t ui;
	float    f;
	double   d;
	void    *array;
    } value;
} bam_aux_tag_t;

/*
 * Our bam stream consists of a zlib gzFile stream and a buffer for it to
 * output to. This allows us to call small bam_read requests while
 * translating these to fewer, larger, gzread calls. The overhead is
 * therefore minimal.
 */
#define Z_BUFF_SIZE 65536    /* Max size of a zlib block */
#define BGZF_BUFF_SIZE 65273 // 65535 - MIN_LOOKAHEAD to avoid fill_window()
typedef struct {
    FILE *fp;
    int mode, binary, level;
    z_stream s;

    unsigned char comp[Z_BUFF_SIZE];
    unsigned char *comp_p;
    size_t comp_sz;

    unsigned char uncomp[Z_BUFF_SIZE];
    unsigned char *uncomp_p;
    size_t uncomp_sz;

    /* BAM specifics */
    int32_t next_len;

    SAM_hdr *header;      /* Parsed SAM header */

    /* Cached bam_seq_t, to avoid excessive mallocs */
    bam_seq_t *bs;
    int bs_size;

    /* Boolean to indicate if we've finished the most recent z stream */
    int z_finish;

    /* Indicates whether gzipped, and if so with bgzf extra fields */
    int gzip;
    int bgzf; 

    /* Whether BAM or SAM format */
    int bam;

    /* If true, skip auxillary field parsing while reading SAM */
    int no_aux;

    /* line number (when in SAM mode) */
    int line;

    /* EOF block present in BAM */
    int eof_block;

    /* Static avoidance: used in sam_next_seq() */
    unsigned char *sam_str;
    size_t alloc_l;

    /* Thread pool for encoding */
    t_pool *pool;
    t_results_queue *equeue;

    /* Decoding queue */
    t_results_queue *dqueue;
    void *job_pending;
    int eof;
    int nd_jobs, ne_jobs;

    /* Quality binning */
    enum quality_binning binning;
} bam_file_t;

/* BAM flags */
#define BAM_FPAIRED           1
#define BAM_FPROPER_PAIR      2
#define BAM_FUNMAP            4
#define BAM_FMUNMAP           8
#define BAM_FREVERSE         16
#define BAM_FMREVERSE        32
#define BAM_FREAD1           64
#define BAM_FREAD2          128
#define BAM_FSECONDARY      256
#define BAM_FQCFAIL         512
#define BAM_FDUP           1024
#define BAM_FSUPPLEMENTARY 2048

#define BAM_CIGAR32    32768

/* Decoding the bam_seq_t struct */
#define bam_blk_size(b)       ((b)->blk_size)
#define bam_set_blk_size(b,v) ((b)->blk_size = (v))

#define bam_ref(b)       ((b)->ref)
#define bam_pos(b)       ((b)->pos)
#define bam_mate_ref(b)  ((b)->mate_ref)
#define bam_mate_pos(b)  ((b)->mate_pos)
#define bam_ins_size(b)  ((b)->ins_size)
#define bam_cigar_len(b) (((b)->flag & BAM_CIGAR32 ? ((b)->bin<<16) : 0) + (b)->cigar_len)
#define bam_name_len(b)  ((b)->name_len)
#define bam_map_qual(b)  ((b)->map_qual)
#define bam_bin(b)       ((b)->flag & BAM_CIGAR32 ? 0 : (b)->bin)
#define bam_flag(b)      ((b)->flag)
#define bam_seq_len(b)   ((b)->len)
#define bam_strand(b)    ((bam_flag((b)) & BAM_FREVERSE) != 0)
#define bam_mstrand(b)   ((bam_flag((b)) & BAM_FMREVERSE) != 0)

#define bam_set_ref(b,v)        ((b)->ref = (v))
#define bam_set_pos(b,v)        ((b)->pos = (v))
#define bam_set_mate_ref(b,v)   ((b)->mate_ref = (v))
#define bam_set_mate_pos(b,v)   ((b)->mate_pos = (v))
#define bam_set_ins_size(b,v)   ((b)->ins_size = (v))
#define bam_set_cigar_len(b, v) (((v)>>16) ? (((b)->flag |= BAM_CIGAR32), (b)->bin = ((v)>>16), (b)->cigar_len = (v)&0xffff) : ((b)->cigar_len = (v)))
#define bam_set_name_len(b,v)   ((b)->name_len = (v))
#define bam_set_map_qual(b,v)   ((b)->map_qual = (v))
static inline void bam_set_bin(bam_seq_t *b, uint32_t v) {
    if (!(b->flag & BAM_CIGAR32)) b->bin = v;
}
#define bam_set_flag(b,v)       ((b)->flag = (v))
#define bam_set_seq_len(b,v)    ((b)->len = (v))

#define bam_name(b)      ((char *)(&(b)->data))
#ifdef ALLOW_UAC
#define bam_cigar(b)     ((uint32_t *)(bam_name((b)) + bam_name_len((b))))
#else
#define bam_cigar(b)     ((uint32_t *)(bam_name((b)) + round4(bam_name_len((b)))))
#endif
#define bam_seq(b)       (((char *)bam_cigar((b))) + 4*bam_cigar_len(b))
#define bam_qual(b)      (bam_seq(b) + (int)(((b)->len+1)/2))
#define bam_aux(b)       (bam_qual(b) + (b)->len)

/* Rounds up to the next multiple of 4 */
#define round4(v) (((v-1)&~3)+4)

/* CIGAR operations, taken from samtools bam.h */
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)

enum cigar_op {
    BAM_UNKNOWN=-1,
    BAM_CMATCH=0,
    BAM_CINS=1,
    BAM_CDEL=2,
    BAM_CREF_SKIP=3,
    BAM_CSOFT_CLIP=4,
    BAM_CHARD_CLIP=5,
    BAM_CPAD=6,
    BAM_CBASE_MATCH=7,
    BAM_CBASE_MISMATCH=8
};

/*
 * Whether this cigar op marches along ref, seq or both
 *
 * Op   Ref   Seq
 * M    1     1
 * I    0     1
 * D    1     0
 * N    1     0
 * S    0     1
 * H    0     0
 * P    0     0
 * =    1     1
 * X    1     1
 */
#define BAM_CONSUME_REF(op) ((0x18d>>(op))&1)
#define BAM_CONSUME_SEQ(op) ((0x193>>(op))&1)


/* ----------------------------------------------------------------------
 * Function prototypes
 * We only support reading, so basically we have open, read, close along
 * with some utility functions for querying aux records.
 */

/*! Opens a SAM or BAM file.
 *
 * The mode parameter indicates the file
 * type (if not auto-detecting) and whether it is for reading or
 * writing. Use "rb" or "wb" for reading or writing BAM and "r" or
 * "w" or reading or writing SAM. When writing BAM, the mode may end
 * with a digit from 0 to 9 to indicate the compression to use with 0
 * indicating uncompressed data.
 *
 * @param fn The filename to open or create.
 * @param mode The input/output mode, similar to fopen().
 *
 * @return
 * Returns a bam_file_t pointer on success;
 *         NULL on failure.
 */
bam_file_t *bam_open(const char *fn, const char *mode);

bam_file_t *bam_open_block(const char *blk, size_t blk_size, SAM_hdr *sh);

/*! Closes a SAM or BAM file.
 * 
 * @param b The file to close.
 *
 * @return
 * Retrurns 0 on success;
 *         -1 on failure.
 */
int bam_close(bam_file_t *b);

/*! Deprecated: please use bam_get_seq() instead.
 */
int bam_next_seq(bam_file_t *b, bam_seq_t **bsp);

/*! Reads the next sequence.
 *
 * Fills out the next bam_seq_t struct.
 * This function will alloc and/or grow the memory accordingly, allowing for
 * efficient reuse.
 *
 * @param bsp Must be non-null, but *bsp may be NULL or an existing
 * bam_seq_t pointer.
 *
 * @return
 * Returns 1 on success;
 *         0 on eof;
 *        -1 on error.
 */
int bam_get_seq(bam_file_t *b, bam_seq_t **bsp);

/*!Looks for aux field 'key' and returns the value.
 * The type is the first char and the value is the 2nd character onwards.
 *
 * @return
 * Returns the value for key; NULL if not found.
 */
char *bam_aux_find(bam_seq_t *b, const char *key);

//!Converts an encoded integer value return by bam_aux_find to an integer.
//
//Analogous to the bam_aux2i functions in samtools.
int32_t bam_aux_i(const uint8_t *dat);
float bam_aux_f(const uint8_t *s);
double bam_aux_d(const uint8_t *s);
char bam_aux_A(const uint8_t *s);
char *bam_aux_Z(const uint8_t *s);

//int bam_aux_del(bam_seq_t *b, uint8_t *s); // not implemented yet

/*! Add auxiliary tags to a bam_seq_t structure.
 *
 * Appends a tag onto the end of a bam_seq_t structure.  The tag name
 * is supplied in the 'tag' parameter, and the data type code in 'type'.
 * Valid type codes are [AcCsSiIfdHZ], as described in the SAM specification.
 *
 * The array_len parameter is used for B type (i.e. array) tags.  If
 * array_len is 0, an ordinary non-array tag is added.  If it is greater
 * than zero, a B type tag of the apropriate data type is made.  Arrays
 * of types H and Z are not allowed, so array_len must be zero if these types
 * are specified.
 *
 * data should point to the data to be added.  This should be of appropriate
 * size for the tag type, i.e. (u)int8_t for A, C and c; (u)int16_t for S and s;
 * (u)int32_t for I and i; a float for f; a double for d; a NUL-terminated
 * string for H and Z.  If array_len is greater than zero, then data should
 * point to an array of the given type.  All data should be in the native
 * format for the machine - it will be converted to little-endian if necessary
 * by the function.
 *
 * @param b         Points to the location of a bam_seq_t *.  If (*b)->alloc
 *                  is too small, the bam_seq_t struct will be reallocated.
 *                  Neither b nor *b should be NULL.
 * @param tag       The tag name (RG, NM, etc.)
 * @param type      The tag data type.
 * @param array_len Array length for array tags, zero for ordinary ones.
 * @param data      Pointer to the tag value.
 *
 * @return
 * Returns  0 on success;
 *         -1 on error
 */

int bam_aux_add(bam_seq_t **b, const char tag[2], char type,
		uint32_t array_len, const void *data);

/*! Add multiple auxiliary tags to a bam_seq_t structure
 *
 * bam_aux_add_vec adds one or more tags listed in the bam_aux_tag_t array
 * to a bam_seq_t structure.  The bam_aux_tag_t struct has four elements,
 * the tag name (char[2]), the type code (char), the array length for B tags and
 * the value which is a union.
 * 
 * The type code determines both the data type of the tag, and the member
 * of value used to access the data.  If array_len is zero, data types 'H' and
 * 'Z' are accessed via member value.h or value.z; 'f' via value.f; 'd' via
 * value.d; 'A' via value.a.  Signed integers should have type 'i', and are
 * accessed through value.i.  Unsigned integers should use type 'I' and
 * value.ui.  The actual type used to store the integer will be the smallest
 * that it will fit in, so signed and unsigned integers that fit in one byte
 * will be stored as type 'c' or 'C'.  Similarly, integers that fit in two
 * bytes will be stored as type 's' or 'S' for unsigned.
 *
 * If array_len is non-zero, a B-type (array) tag is stored.  In this case
 * value.array should point to the data to be stored.  The type of array
 * is interpreted according to the requested tag type, i.e. char for 'A';
 * int8_t for 'c'; uint8_t for 'C'; int16_t for 's'; uint16_t for 'S';
 * int32_t for 'i'; uint32_t for 'I'; float for 'f'; double for 'd'.  No
 * attempt is made to adjust the size of integers when storing arrays. All
 * data should be in the native format for the machine - it will be converted
 * to little-endian if necessary by the function.
 *
 * @param b         Points to the location of a bam_seq_t *.  If (*b)->alloc
 *                  is too small, the bam_seq_t struct will be reallocated.
 *                  Neither b nor *b should be NULL.
 * @param count     The number of elements in the tags array
 * @param tags      Array of bam_aux_tag_t structs, listing the tags to add.
 * 
 * @return 
 * Returns  0 on success;
 *         -1 on error
 */

int bam_aux_add_vec(bam_seq_t **b, uint32_t count, bam_aux_tag_t *tags);

/*! Calculate the amount of space needed to store auxiliary tags
 *
 * The tags array should be filled out as described in bam_aux_add_vec.
 * This function iterates through the list of items in the tags array
 * to work out how much space will be needed to store them.  This value
 * can be passed to bam_construct_seq in the extra_len parameter to
 * ensure that enough space is allocated to store the tags.
 *
 * @param count     The number of elements in the tags array
 * @param tags      Array of bam_aux_tag_t structs, listing the tags to add.
 *
 * @return
 * Returns the total space required for the tags on success;
 *         -1 on failure.
 */

ssize_t bam_aux_size_vec(uint32_t count, bam_aux_tag_t *tags);

/*! Helper for filling in the bam_aux_tag_t struct
 * @param tag   Pointer to bam_aux_tag_t struct
 * @param name  Tag name
 * @param val   Tag value
 */
static inline void bam_aux_tag_char(bam_aux_tag_t *tag, char *name, char val) {
    tag->tag[0] = name[0];
    tag->tag[1] = name[1];
    tag->type = 'A';
    tag->array_len = 0;
    tag->value.a = val;
}
/*! Helper for filling in the bam_aux_tag_t struct
 * @param tag   Pointer to bam_aux_tag_t struct
 * @param name  Tag name
 * @param val   Tag value
 */
static inline void bam_aux_tag_int(bam_aux_tag_t *tag,
				   char *name, int32_t val) {
    tag->tag[0] = name[0];
    tag->tag[1] = name[1];
    tag->type = 'i';
    tag->array_len = 0;
    tag->value.i = val;
}
/*! Helper for filling in the bam_aux_tag_t struct
 * @param tag   Pointer to bam_aux_tag_t struct
 * @param name  Tag name
 * @param val   Tag value
 */
static inline void bam_aux_tag_uint(bam_aux_tag_t *tag,
				    char *name, uint32_t val) {
    tag->tag[0] = name[0];
    tag->tag[1] = name[1];
    tag->type = 'I';
    tag->array_len = 0;
    tag->value.ui = val;
}
/*! Helper for filling in the bam_aux_tag_t struct
 * @param tag   Pointer to bam_aux_tag_t struct
 * @param name  Tag name
 * @param val   Tag value
 */
static inline void bam_aux_tag_float(bam_aux_tag_t *tag,
				     char *name, float val) {
    tag->tag[0] = name[0];
    tag->tag[1] = name[1];
    tag->type = 'f';
    tag->array_len = 0;
    tag->value.f = val;
}
/*! Helper for filling in the bam_aux_tag_t struct
 * @param tag   Pointer to bam_aux_tag_t struct
 * @param name  Tag name
 * @param val   Tag value
 */
static inline void bam_aux_tag_double(bam_aux_tag_t *tag,
				      char *name, double val) {
    tag->tag[0] = name[0];
    tag->tag[1] = name[1];
    tag->type = 'd';
    tag->array_len = 0;
    tag->value.d = val;
}
/*! Helper for filling in the bam_aux_tag_t struct
 * @param tag   Pointer to bam_aux_tag_t struct
 * @param name  Tag name
 * @param val   Tag value
 */
static inline void bam_aux_tag_string(bam_aux_tag_t *tag,
				      char *name, char *str) {
    tag->tag[0] = name[0];
    tag->tag[1] = name[1];
    tag->type = 'Z';
    tag->array_len = 0;
    tag->value.z = str;
}
/*! Helper for filling in the bam_aux_tag_t struct
 * @param tag   Pointer to bam_aux_tag_t struct
 * @param name  Tag name
 * @param val   Tag value
 */
static inline void bam_aux_tag_hexstring(bam_aux_tag_t *tag,
					 char *name, uint8_t *str) {
    tag->tag[0] = name[0];
    tag->tag[1] = name[1];
    tag->type = 'H';
    tag->array_len = 0;
    tag->value.h = str;
}
/*! Helper for filling in the bam_aux_tag_t struct
 * @param tag   Pointer to bam_aux_tag_t struct
 * @param name  Tag name
 * @param val   Tag value
 */
static inline void bam_aux_tag_array(bam_aux_tag_t *tag, char *name, char type,
				     uint32_t array_len, void *data) {
    tag->tag[0] = name[0];
    tag->tag[1] = name[1];
    tag->type = type;
    tag->array_len = array_len;
    tag->value.array = data;
}

/*! Add SAM formatted aux tags to a bam_seq_t.
 *
 * Appends one or more SAM-format tags onto the end of a bam_seq_t structure.
 * If multiple tags are present, they should be separated by tabs, as in
 * a SAM file.
 *
 * @param bsp       Points to the location of a bam_seq_t *.  If (*bsp)->alloc
 *                  is too small, the bam_seq_t struct will be reallocated.
 *                  Neither bsp nor *bsp should be NULL.
 * @param sam       SAM-foratted tag string.
 *
 * @return
 * Returns  0 on success;
 *         -1 on error
 */

int bam_aux_add_from_sam(bam_seq_t **bsp, char *sam);

/*! Add preformated raw aux data to the bam_seq.
 *
 * This interface is similar to the samtools bam_aux_append function.
 * It creates a tag with the given name and type, and then appends the
 * supplied data to it as the value.  It is up to the caller to ensure that
 * the data has been formatted correctly as given in the SAM specification.
 * This function makes no checks on the data, but simply copies it.
 *
 * Consider using bam_aux_add instead if you have information in a more
 * integer or string form.
 *
 * @param b         Points to the location of a bam_seq_t *.  If (*b)->alloc
 *                  is too small, the bam_seq_t struct will be reallocated.
 *                  Neither b nor *b should be NULL.
 * @param tag       The tag name (RG, NM, etc.)
 * @param type      The tag data type.
 * @param len       The number of bytes of data present.
 * @param data      Pre-formatted data.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int bam_aux_add_data(bam_seq_t **b, const char tag[2],
		     char type, size_t len, const uint8_t *data);

/*! Add raw data to a bam structure.
 *
 * This could be useful if you wish to manually construct your own bam
 * entries or if you need to append an entire block of preformatting
 * aux data.
 *
 * @param b         Points to the location of a bam_seq_t *.  If (*b)->alloc
 *                  is too small, the bam_seq_t struct will be reallocated.
 *                  Neither b nor *b should be NULL.
 * @param len       The number of bytes of data present.
 * @param data      Pre-formatted data.
 *
 * @retrun
 * Returns 0 on success;
 *        -1 on failure
 */
int bam_add_raw(bam_seq_t **b, size_t len, const uint8_t *data);

/*! An iterator on bam_aux_t fields.
 *
 * NB: This code is not reentrant or multi-thread capable. The values
 * returned are valid until the next call to this function.
 *
 * @param key  points to an array of 2 characters (eg "RG", "NM")
 * @param type points to an address of 1 character (eg 'Z', 'i')
 * @paran val  points to an address of a bam_aux_t union.
 * @param iter_handle NULL to initialise the search, and then the
 * returned (modified) iter_handle on each subsequent call to continue
 * the iteration.
 *
 * @return
 * Returns 0 if the next value is valid, setting key, type and val;
 *        -1 when no more found.
 */
int bam_aux_iter(bam_seq_t *b, char **iter_handle,
		 char *key, char *type, bam_aux_t *val);

/* Taken from samtools/bam.h */
#define bam_seqi(s, i) ((s)[(i)/2] >> 4*(1-(i)%2) & 0xf)
#define bam_nt16_rev_table "=ACMGRSVTWYHKDBN"

/* Output code */

/*! Writes a single bam sequence object.
 *
 * @param fp The SAM/BAM file handle.
 * @param b  The bam_seq_t pointer
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int bam_put_seq(bam_file_t *fp, bam_seq_t *b);

/*! Constructs a bam_seq_t from separate components.
 *
 * Note: ignores auxiliary tags for now. These need to be appended
 * manually by the calling function.
 *
 * @param b         Points to the location of a bam_seq_t *.  If *b is NULL
 *                  or (*b)->alloc is too small, the bam_seq_t struct will
 *                  be reallocated.
 * @param extra_len Extra space to allocate for auxiliary tags
 * @param qname     Query name
 * @param qname_len Query name length
 * @param flag      BAM flags
 * @param rname     Reference ID
 * @param pos       Mapped position (N.B. 1-based)
 * @param end       Last aligned base
 * @param mapq      Mapping quality
 * @param ncigar    Number of CIGAR elements
 * @param cigar     CIGAR alignment information
 * @param mrnm      Mate reference ID
 * @param mpos      Mate position (N.B. 1-based)
 * @param isize     Insert size
 * @param len       Sequence length
 * @param seq       Sequence (ASCII format)
 * @param qual      Quality values (phred scale, 8-bit binary, no offset)
 *                  Passing in NULL to qual will cause all quality values
 *                  to be treated as absent (i.e. set to 0xff).
 *
 * @return
 * Returns number of bytes written to bam_seq_t on success (ie tag offset);
 *         -1 on error.
 */
int bam_construct_seq(bam_seq_t **b, size_t extra_len,
		      const char *qname, size_t qname_len,
		      int flag,
		      int rname,      // Ref ID
		      int pos,
		      int end, // aligned start/end coords
		      int mapq,
		      uint32_t ncigar, const uint32_t *cigar,
		      int mrnm,       // Mate Ref ID
		      int mpos,
		      int isize,
		      int len,
		      const char *seq,
		      const char *qual);

/*! Duplicates a bam_seq_t structure.
 *
 * @return
 * Returns the new bam_seq_t pointer on success;
 *         NULL on failure.
 */
bam_seq_t *bam_dup(bam_seq_t *b);

/*! Writes a SAM header block.
 *
 * @return
 * Returns 0 for success;
 *        -1 for failure
 */
int bam_write_header(bam_file_t *out);

enum bam_option {
    BAM_OPT_THREAD_POOL,
    BAM_OPT_BINNING
};

/*! Sets options on the bam_file_t.
 *
 * Sets options on the bam_file_t. See BAM_OPT_* definitions in bam.h.
 * Use this immediately after opening.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int bam_set_option(bam_file_t *fd, enum bam_option opt, ...);

/*! Sets options on the bam_file_t.
 *
 * Sets options on the bam_file_t. See BAM_OPT_* definitions in bam.h.
 * Use this immediately after opening.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int bam_set_voption(bam_file_t *fd, enum bam_option opt, va_list args);



/*
 * Signed and unsigned fast functions to act as equiv to sprintf(cp, "%d", i)
 */
unsigned char *append_int(unsigned char *cp, int32_t i);
unsigned char *append_uint(unsigned char *cp, uint32_t i);

#ifdef __cplusplus
}
#endif

#endif /* _BAM_H_ */
