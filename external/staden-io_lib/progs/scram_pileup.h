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

#ifndef _PILEUP_H_
#define _PILEUP_H_

#include "io_lib/scram.h"

typedef struct pileup {
    struct pileup *next;  // A link list, for active seqs
    void *cd;		  // General purpose per-seq client-data

    bam_seq_t *b;	  // Bam entry associated with struct
    unsigned char *b_qual;// cached bam_qual
    unsigned char *b_seq; // cached bam_seq
    uint32_t *b_cigar;    // cached bam_cigar
    int  b_strand;        // 0 => fwd, 1 => rev

    int  pos;             // Current unpadded position in seq
    int  nth;		  // nth base at unpadded position 'pos'
    int  seq_offset;      // Current base position in s->seq[] array.

    int  cigar_ind;       // Current location in s->alignment cigar str
    int  cigar_op;        // Current cigar operation
    int  cigar_len;       // Remaining length of this cigar op

    int  first_del;       // Used when first base is a deletion

    int  eof;		  // True if this sequence has finished
    int  qual;            // Current qual (for active seq only)
    char base;		  // Current base (for active seq only)
    char start;		  // True if this is a new sequence
    char ref_skip;        // True if the cause of eof or start is cigar N
    char padding;         // True if the base was added due to another seq
} pileup_t;

int pileup_loop(scram_fd *fp,
		int (*seq_init)(void *client_data,
				scram_fd *fp,
				pileup_t *p),
		int (*seq_add)(void *client_data,
			       scram_fd *fp,
			       pileup_t *p,
			       int depth,
			       int pos,
			       int nth,
			       int is_insert),
		void *client_data);

#endif
