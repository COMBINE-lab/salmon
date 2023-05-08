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

/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2011-2013
 *
 * A basic pileup command.
 * This needs a lot of work still and to be split into library + program.
 * For now we have all in one place.
 *
 * Compatibility wise it is not meant to be a full blown replacement for
 * samtools mpileup or even the old samtools pileup. Primarily it is a test
 * harness for the pileup_loop API code.
 *
 * Speed wise it is approaching double the performance of samtools mpileup
 * (with no extra arguments) when run on BAM. CRAM performance is approx 40%
 * slower than the BAM version.
 */


#include <io_lib_config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include "scram_pileup.h"

/*
 * START_WITH_DEL is the mode that Gap5 uses when building this. It prepends
 * all cigar strings with 1D and decrements the position by one. (And then
 * has code to reverse this operation in the pileup handler.)
 *
 * The reason for this is that it means reads starting with an insertion work.
 * Otherwise the inserted bases are silently lost. (Try it with "samtools
 * mpileup" and you can see it has the same issue.)
 *
 * However it's probably not want most people expect.
 */
//#define START_WITH_DEL

/* --------------------------------------------------------------------------
 * The pileup code itself. 
 *
 * This consists of the external pileup_loop() function, which takes a
 * sam/bam samfile_t pointer and a callback function. The callback function
 * is called once per column of aligned data (so once per base in an
 * insertion).
 *
 * Current known issues.
 * 1) zero length matches, ie 2S2S cause failures.
 * 2) Insertions at starts of sequences get included in the soft clip, so
 *    2S2I2M is treated as if it's 4S2M
 * 3) From 1 and 2 above, 1S1I2S becomes 2S2S which fails.
 */

/*
 * Fast conversion from encoded SAM base nibble to a printable character
 */
static char tab[256][2];
static void init_tab(void) {
    int i, j;
    unsigned char b2;
    static int done = 0;

    if (done)
	return;

    for (i = 0; i < 16; i++) {
	for (j = 0; j < 16; j++) {
	    b2 = (i<<4) | j;
	    tab[b2][0] = "NACMGRSVTWYHKDBN"[i];
	    tab[b2][1] = "NACMGRSVTWYHKDBN"[j];
	}
    }

    done = 1;
}


/*
 * Fetches the next base => the nth base at unpadded position pos. (Nth can
 * be greater than 0 if we have an insertion in this column). Do not call this
 * with pos/nth lower than the previous query, although higher is better.
 * (This allows it to be initialised at base 0.)
 *
 * Stores the result in base and also updates is_insert to indicate that
 * this sequence still has more bases in this position beyond the current
 * nth parameter.
 *
 * Returns 1 if a base was fetched
 *         0 if not (eg ran off the end of sequence)
 */
static int get_next_base(pileup_t *p, int pos, int nth, int *is_insert) {
    bam_seq_t *b = p->b;
    enum cigar_op op = p->cigar_op;

    if (p->first_del && op != BAM_CPAD)
	p->first_del = 0;

    *is_insert = 0;

    /* Find pos first */
    while (p->pos < pos) {
	p->nth = 0;

	if (p->cigar_len == 0) {
	    if (p->cigar_ind >= bam_cigar_len(b)) {
		p->eof = 1;
		return 0;
	    }

	    op=p->cigar_op  = p->b_cigar[p->cigar_ind] & BAM_CIGAR_MASK;
	    p->cigar_len = p->b_cigar[p->cigar_ind] >> BAM_CIGAR_SHIFT;
	    p->cigar_ind++;
	}
	
	if ((op == BAM_CMATCH ||
	     op == BAM_CBASE_MATCH ||
	     op == BAM_CBASE_MISMATCH) && p->cigar_len <= pos - p->pos) {
	    p->seq_offset += p->cigar_len;
	    p->pos += p->cigar_len;
	    p->cigar_len = 0;
	} else {
	    switch (op) {
	    case BAM_CMATCH:
	    case BAM_CBASE_MATCH:
	    case BAM_CBASE_MISMATCH:
		p->seq_offset++;
		/* Fall through */
	    case BAM_CDEL:
	    case BAM_CREF_SKIP:
		p->pos++;
		p->cigar_len--;
		break;

	    case BAM_CINS:
	    case BAM_CSOFT_CLIP:
		p->seq_offset += p->cigar_len;
		/* Fall through */
	    case BAM_CPAD:
	    case BAM_CHARD_CLIP:
		p->cigar_len = 0;
		break;

	    default:
		fprintf(stderr, "Unhandled cigar_op %d\n", op);
		return -1;
	    }
	}
    }

    /* Now at pos, find nth base */
    while (p->nth < nth) {
	if (p->cigar_len == 0) {
	    if (p->cigar_ind >= bam_cigar_len(b)) {
		p->eof = 1;
		return 0; /* off end of seq */
	    }

	    op=p->cigar_op  = p->b_cigar[p->cigar_ind] & BAM_CIGAR_MASK;
	    p->cigar_len = p->b_cigar[p->cigar_ind] >> BAM_CIGAR_SHIFT;
	    p->cigar_ind++;
	}

	switch (op) {
	case BAM_CMATCH:
	case BAM_CBASE_MATCH:
	case BAM_CBASE_MISMATCH:
	case BAM_CSOFT_CLIP:
	case BAM_CDEL:
	case BAM_CREF_SKIP:
	    goto at_nth; /* sorry, but it's fast! */

	case BAM_CINS:
	    p->seq_offset++;
	    /* Fall through */
	case BAM_CPAD:
	    p->cigar_len--;
	    p->nth++;
	    break;

	case BAM_CHARD_CLIP:
	    p->cigar_len = 0;
	    break;

	default:
	    fprintf(stderr, "Unhandled cigar_op %d\n", op);
	    return -1;
	}
    }
 at_nth:

    /* Fill out base & qual fields */
    p->ref_skip = 0;
    if (p->nth < nth && op != BAM_CINS) {
	//p->base = '-';
	p->base = '*';
	p->padding = 1;
	if (p->seq_offset < b->len)
	    p->qual = (p->qual + p->b_qual[p->seq_offset+1])/2;
	else
	    p->qual = 0;
    } else {
	p->padding = 0;
	switch(op) {
	case BAM_CDEL:
	    p->base = '*';
	    if (p->seq_offset+1 < b->len)
		p->qual = (p->qual + p->b_qual[p->seq_offset+1])/2;
	    else
		p->qual = (p->qual + p->b_qual[p->seq_offset])/2;
	    break;

	case BAM_CPAD:
	    //p->base = '+';
	    p->base = '*';
	    if (p->seq_offset+1 < b->len)
		p->qual = (p->qual + p->b_qual[p->seq_offset+1])/2;
	    else
		p->qual = (p->qual + p->b_qual[p->seq_offset])/2;
	    break;

	case BAM_CREF_SKIP:
	    p->base = '.';
	    p->qual = 0;
	    /* end of fragment, but not sequence */
	    p->eof = p->eof ? 2 : 3;
	    p->ref_skip = 1;
	    break;

	default:
	    if (p->seq_offset < b->len) {
		p->qual = p->b_qual[p->seq_offset];
	    /*
	     * If you need to label inserted bases as different from
	     * (mis)matching bases then this is where we'd make that change.
	     * The reason could be to allow the consensus algorithm to easily
	     * distinguish between reference bases and non-reference bases.
	     *
	     * Eg:
	     * if (nth)
	     *     p->base = tolower(tab[p->b_seq[p->seq_offset/2]][p->seq_offset&1]);
	     * else
	     */
		p->base = tab[p->b_seq[p->seq_offset/2]][p->seq_offset&1];
	    } else {
		p->base = 'N';
		p->qual = 0xff;
	    }
		
	    break;
	}
    }

    /* Handle moving out of N (skip) into sequence again */
    if (p->eof && p->base != '.') {
	p->start = 1;
	p->ref_skip = 1;
	p->eof = 0;
    }

    /* Starting with an indel needs a minor fudge */
    if (p->start && p->cigar_op == BAM_CDEL) {
	p->first_del = 1;
    }

    /* Check if next op is an insertion of some sort */
    if (p->cigar_len == 0) {
	if (p->cigar_ind < bam_cigar_len(b)) {
	    op=p->cigar_op  = p->b_cigar[p->cigar_ind] & BAM_CIGAR_MASK;
	    p->cigar_len = p->b_cigar[p->cigar_ind] >> BAM_CIGAR_SHIFT;
	    p->cigar_ind++;
	    if (op == BAM_CREF_SKIP) {
		p->eof = 3;
		p->ref_skip = 1;
	    }
	} else {
	    p->eof = 1;
	}
    }

    switch (op) {
    case BAM_CPAD:
    case BAM_CINS:
	*is_insert = p->cigar_len;
	break;

    case BAM_CSOFT_CLIP:
	/* Last op 'S' => eof */
	p->eof = p->cigar_ind == bam_cigar_len(b) ? 1 : 0;
	break;

    case BAM_CHARD_CLIP:
	p->eof = 1;
	break;

    default:
	break;
    }
    
    return 1;
}

/*
 * Loops through a set of supplied ranges producing columns of data.
 * When found, it calls func with clientdata as a callback. Func should
 * return 0 for success and non-zero for failure. seq_init() is called
 * on each new entry before we start processing it. It should return 0 or 1
 * to indicate reject or accept status (eg to filter unmapped data).
 * If seq_init() returns -1 we abort the pileup_loop with an error.
 * seq_init may be NULL.
 *
 * Returns 0 on success
 *        -1 on failure
 */
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
		void *client_data) {
    int ret = -1;
    pileup_t *phead = NULL, *p, *pfree = NULL, *last, *next, *ptail = NULL;
    pileup_t *pnew = NULL;
    int is_insert, nth = 0;
    int col = 0, r;
    int last_ref = -1;

    /* FIXME: allow for start/stop boundaries rather than consuming all data */

    init_tab();
    if (NULL == (pnew = calloc(1, sizeof(*p))))
	return -1;
    
    do {
	bam_seq_t *b;
	int pos, last_in_contig;

	r = scram_next_seq(fp, &pnew->b);
	if (r == -1) {
	    //fprintf(stderr, "bam_next_seq() failure on line %d\n", fp->line);
	    if (!scram_eof(fp)) {
		fprintf(stderr, "bam_next_seq() failure.\n");
		return -1;
	    }
	}

	b = pnew->b;

	/* Force realloc */
	//fp->bs = NULL;
	//fp->bs_size = 0;

	//r = samread(fp, pnew->b);
	if (r >= 0) {
	    if (bam_flag(b) & BAM_FUNMAP)
		continue;
	    
	    if (b->ref == -1) {
		/* Another indicator for unmapped */
		continue;
	    } else if (b->ref == last_ref) {
		pos = b->pos+1;
		//printf("New seq at pos %d @ %d %s\n", pos, b->ref,
		//       bam_name(b));
		last_in_contig = 0;
	    } else {
		//printf("New ctg at pos %d @ %d\n",b->pos+1,b->ref);
		pos = (b->pos > col ? b->pos : col)+1;
		last_in_contig = 1;
	    }
	} else {
	    last_in_contig = 1;
	    pos = col+1;
	}

	if (col > pos) {
	    fprintf(stderr, "BAM/SAM file is not sorted by position. "
		    "Aborting\n");
	    return -1;
	}

	/* Process data between the last column and our latest addition */
	while (col < pos && phead) {
	    int v, ins, depth = 0;

	    //printf("Col=%d pos=%d nth=%d\n", col, pos, nth);

	    /* Pileup */
	    is_insert = 0;
	    for (p = phead; p; p = p->next) {
		if (!get_next_base(p, col, nth, &ins))
		    p->eof = 1;

		if (is_insert < ins)
		    is_insert = ins;
		
		depth++;
	    }

	    /* Call our function on phead linked list */
#ifdef START_WITH_DEL
	    v = seq_add(client_data, fp, phead, depth, col-1, nth, is_insert);
#else
	    v = seq_add(client_data, fp, phead, depth, col, nth, is_insert);
#endif

	    /* Remove dead seqs */
	    for (p = phead, last = NULL; p; p = next) {
		next = p->next;
		
		p->start = 0;
		if (p->eof == 1) {
		    if (last)
			last->next = p->next;
		    else
			phead = p->next;

		    p->next = pfree;
		    pfree = p;
		    
		    //printf("Del seq %s at pos %d\n", bam_name(p->b), col);
		} else {
		    last = p;
		}
	    }
	    if ((ptail = last) == NULL)
		ptail = phead;

	    if (v == 1)
		break; /* early abort */
	    
	    if (v != 0)
		goto error;

	    /* Next column */
	    if (is_insert) {
		nth++;
	    } else {
		nth = 0;
		col++;
	    }

	    /* Special case for the last sequence in a contig */
	    if (last_in_contig && phead)
		pos++;
	}

	/* May happen if we have a hole in the contig */
	col = pos;

	/* New contig */
	if (b && b->ref != last_ref) {
	    last_ref = b->ref;
	    pos = b->pos+1;
	    nth = 0;
	    col = pos;
	}

	/*
	 * Add this seq.
	 * Note: cigars starting with I or P ops (eg 2P3I10M) mean we have
	 * alignment instructions that take place before the designated
	 * starting location listed in the SAM file. They won't get included
	 * in the callback function until they officially start, which is
	 * already too late.
	 *
	 * So to workaround this, we prefix all CIGAR with 1D, move the
	 * position by 1bp, and then force the callback code to remove
	 * leaving pads (either P or D generated).
	 *
	 * Ie it's a level 10 hack!
	 */
	if (r >= 0) {
	    p = pnew;
	    p->next       = NULL;
	    p->cd         = NULL;
	    p->start      = 1;
	    p->eof        = 0;
#ifdef START_WITH_DEL
	    p->pos        = pos-1;
	    p->cigar_ind  = 0;
	    p->b_cigar    = bam_cigar(p->b);
	    if ((p->b_cigar[0] & BAM_CIGAR_MASK) == BAM_CHARD_CLIP) {
		p->cigar_len  = p->b_cigar[0] >> BAM_CIGAR_SHIFT;
		p->cigar_op   = BAM_CHARD_CLIP;
		if ((p->b_cigar[1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
		    /* xHxS... => xHxS1D... */
		    p->b_cigar[0] = p->b_cigar[1];
		    p->b_cigar[1] = (1 << BAM_CIGAR_SHIFT) | BAM_CDEL;
		} else {
		    /* xH... => xH1D... */
		    p->b_cigar[0] = (1 << BAM_CIGAR_SHIFT) | BAM_CDEL;
		}
	    } else {
		if ((p->b_cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
		    /* xS... => xS1D... */
		    p->cigar_len  = p->b_cigar[0] >> BAM_CIGAR_SHIFT;
		    p->cigar_op   = BAM_CSOFT_CLIP;
		    p->b_cigar[0] = (1 << BAM_CIGAR_SHIFT) | BAM_CDEL;
		} else {
		    /* ... => 1D... */
		    p->cigar_len  = 1;        /* was  0  */
		    p->cigar_op   = BAM_CDEL; /* was 'X' */
		}
	    }
	    p->seq_offset = -1;
	    p->first_del  = 1;
#else
	    p->pos        = pos-1;
	    p->cigar_ind  = 0;
	    p->b_cigar    = bam_cigar(p->b);
	    p->cigar_len  = 0;
	    p->cigar_op   = -1;
	    p->seq_offset = -1;
	    p->first_del  = 0;
#endif
	    p->b_strand   = bam_strand(p->b) ? 1 : 0;
	    p->b_qual     = (uc *)bam_qual(p->b);
	    p->b_seq      = (uc *)bam_seq(p->b);

	    if (seq_init) {
		int v;
		v = seq_init(client_data, fp, p);
		if (v == -1)
		    return -1;
		
		if (v == 1) {
		    /* Keep this seq */
		    if (phead) {
			ptail->next = p;
		    } else {
			phead = p;
		    }
		    ptail = p;
		} else {
		    /* Push back on free list */
		    p->next = pfree;
		    pfree = p;
		}
	    } else {
		if (phead)
		    ptail->next = p;
		else
		    phead = p;
		ptail = p;
	    }

	    /* Allocate the next pileup rec */
	    if (pfree) {
		pnew = pfree;
		pfree = pfree->next;
	    } else {
		if (NULL == (pnew = calloc(1, sizeof(*pnew))))
		    goto error;
	    }
	}
    } while (r >= 0);

    ret = 0;
 error:

    if (pnew) {
	free(pnew->b);
	free(pnew);
    }

    /* Tidy up */
    for (p = pfree; p; p = next) {
	next = p->next;
	free(p->b);
	free(p);
    }

    return ret;
}

/* --------------------------------------------------------------------------
 * Example usage of the above pileup code
 */

#include <ctype.h>
#include <io_lib/bam.h>

char strand_char[2][256];
void strand_init(void) {
    int i;
    for (i = 0; i < 256; i++) {
	strand_char[0][i] = toupper((unsigned char)i);
	strand_char[1][i] = tolower((unsigned char)i);
    }
}

/*
 * Used to delay emitting pileup lines around insertions, so we can
 * stack multiple bases together.
 */
typedef struct {
    int alloc;
    char *base;       // first base call
    int  *seq_offset; // first seq_offset
    int  *seq_len;    // length of insertion
} sam_pileup_t;

static int sam_pileup(void *cd_v, scram_fd *fp, pileup_t *p,
		      int depth, int pos, int nth, int is_insert) {
    static unsigned char *seq = NULL, *qual = NULL, *buf = NULL;
    static size_t seq_alloc = 0, buf_alloc = 0;
    static int max_depth = 0;
    unsigned char *sp, *qp, *cp;
    int ref;
    sam_pileup_t *cd = (sam_pileup_t *)cd_v;
    size_t buf_len;

    if (max_depth < depth) {
	max_depth = depth;
	seq  = realloc(seq,  seq_alloc = max_depth*2);
	qual = realloc(qual, max_depth);
	//buf  = realloc(buf,  max_depth*2+1000);

	if (!seq || !qual)
	    return -1;
    }

    sp = seq; qp = qual; cp = buf;

    if (!p)
	return 0;

    if (is_insert) {
	pileup_t *p_orig = p;
	int i;

	if (nth == 0) {
	    /* Reference position pos with is_insert inserted bases to come */
	    if (cd->alloc < depth) {
		cd->alloc = depth;
		cd->base       = realloc(cd->base,
					 cd->alloc * sizeof(*cd->base));
		cd->seq_offset = realloc(cd->seq_offset,
					 cd->alloc * sizeof(*cd->seq_offset));
		cd->seq_len    = realloc(cd->seq_len,
					 cd->alloc * sizeof(*cd->seq_len));
	    }

	    /*
	     * FIXME: This assumes that the number of entries in the p list
	     * here matches the number of entries in the p list when is_inser
	     * becomes 0 (ie the last base of the insert).
	     *
	     * Given that sequences can end mid-way through an insert or
	     * even start mid-way, this assumption is false.
	     */
	    for (i = 0; p; p = p->next, i++) {
		cd->base[i]       = strand_char[p->b_strand][(uc)p->base];
		cd->seq_offset[i] = p->seq_offset+1;
		cd->seq_len[i]    = 0;
	    }
	    p = p_orig;
	} else {
	    for (i = 0; p; p = p->next, i++)
		if (p->base != '*')
		    cd->seq_len[i]++;
	}
	p = p_orig;

	return 0;
    }

    ref = p->b->ref;

    if (nth) {
	int n;
	for (n = 0; p; p = p->next, n++) {
	    int i, j;
	    uint8_t *b_seq = (uint8_t *)bam_seq(p->b);

	    while ((sp - seq + 5 + cd->seq_len[n]) > seq_alloc) {
		ptrdiff_t d = sp - seq;
		seq = realloc(seq, seq_alloc*=2);
		sp = seq + d;
	    }

	    if(p->base != '*')
		cd->seq_len[n]++;

	    if (p->start) {
		*sp++ = '^';
		*sp++ = MIN(p->b->map_qual,93) + '!';
	    }
	    *sp++ = cd->base[n];
	    if (cd->seq_len[n]) {
		*sp++ = '+';
		sp = append_int(sp, cd->seq_len[n]);

		for (i = cd->seq_offset[n], j = 0; j < cd->seq_len[n]; i++, j++) {
		    uc call = bam_nt16_rev_table[bam_seqi(b_seq, i)];
		    *sp++ = strand_char[bam_strand(p->b)][call];
		}
	    }
	    if (p->eof)
		*sp++ = '$';
	    *qp++ = MIN(p->qual,93) + '!';
	}
    } else {
	for (; p; p = p->next) {
	    while ((sp - seq + 4) > seq_alloc) {
		ptrdiff_t d = sp - seq;
		seq = realloc(seq, seq_alloc*=2);
		sp = seq + d;
	    }
	    if (p->start) {
		*sp++ = '^';
		*sp++ = MIN(p->b->map_qual,93) + '!';
	    }
	    *sp++ = strand_char[p->b_strand][(uc)p->base];
	    //*sp++ = strand_char[bam_strand(p->b)][p->base];
	    if (p->eof)
		*sp++ = '$';
	    *qp++ = MIN(p->qual,93) + '!';
	}
    }

    /* Equivalent to the printf below, but faster */
    buf_len = strlen(scram_get_header(fp)->ref[ref].name) + 1 // name
	+ 10 + 1                                              // pos
	+ 1  + 1                                              // base
	+ 10 + 1                                              // depth
	+ sp - seq + 1                                        // seq
	+ qp - qual + 1;                                      // qual
    if (buf_len > buf_alloc)
	buf = realloc(buf, buf_alloc = buf_len);

    cp = buf;
    strcpy((char *) cp, scram_get_header(fp)->ref[ref].name);
    cp += strlen((char *) cp);
    *cp++ = '\t';
    cp = append_int(cp, pos);   *cp++ = '\t';
    *cp++ = 'N';
    *cp++ = '\t';
    cp = append_int(cp, depth); *cp++ = '\t';
    memcpy(cp, seq,  sp-seq);  cp += sp-seq;  *cp++ = '\t';
    memcpy(cp, qual, qp-qual); cp += qp-qual; *cp++ = '\0';
    puts((char *) buf);

    //*sp++ = 0;
    //*qp++ = 0;
    //printf("ref\t%d+%d\tN\t%d\t%s\t%s\n", pos, nth, depth, seq, qual);

    return 0;
}

static int basic_pileup(void *cd, scram_fd *fp, pileup_t *p,
			int depth, int pos, int nth, int is_insert) {
    static unsigned char *seq = NULL, *qual = NULL, *buf = NULL;
    static int max_depth = 0;
    unsigned char *qp, *cp, *rp;
    int ref;

    if (max_depth < depth) {
	max_depth = depth;
	seq  = realloc(seq,  max_depth*3);
	qual = realloc(qual, max_depth);
	buf  = realloc(buf,  max_depth*2+1000);

	if (!seq || !qual || !buf)
	    return -1;
    }

    qp = qual; cp = buf;

    if (!p)
	return 0;

    /* Ref, pos, depth */
    ref = p->b->ref;
    rp = (unsigned char *) scram_get_header(fp)->ref[ref].name;
    while ((*cp++ = *rp++))
	;
    cp--;
    *cp++ = '\t';
    cp = append_int(cp, pos);   *cp++ = '+';
    cp = append_int(cp, nth);   *cp++ = '\t';
    *cp++ = 'N';
    *cp++ = '\t';
    cp = append_int(cp, depth); *cp++ = '\t';

    /* Seq + qual at predetermined offsets */
    qp = cp + depth + 1;
    for (; p; p = p->next) {
	*cp++ = p->base;
	*qp++ = MIN(p->qual,93) + '!';
    }
    *cp++ = '\t';
    *qp++ = '\0';
    puts((char *) buf);

    return 0;
}

static int depth_pileup(void *cd, scram_fd *fp, pileup_t *p,
			int depth, int pos, int nth, int is_insert) {
    unsigned char buf[1024], *cp = buf, *rp;

    if (nth)
	return 0;

    rp = (unsigned char *) scram_get_header(fp)->ref[p->b->ref].name;
    while ((*cp++ = *rp++))
	;
    cp[-1] = '\t';
    cp = append_int(cp, pos);
    *cp++=  '\t';
    cp = append_int(cp, depth);
    *cp++ = '\0';
    puts((char *) buf);
    
    return 0;
}

int main(int argc, char **argv) {
    scram_fd *fp;
    sam_pileup_t *p;
    int mode = 0;

    if (argc < 2) {
	fprintf(stderr, "Usage: scram_pileup [options] filename.{sam,bam,cram}\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, " -5          Gap5 pileup format.\n");
	fprintf(stderr, " -d          Depth format.\n");
	fprintf(stderr, " (otherwise) Samtools pileup format.\n");
	fprintf(stderr, "\n\nNOTE: This program is still under development "
		"and should be considered a proof\nof concept only.\n");
	return 1;
    }

    if (argc >= 2 && strcmp(argv[1], "-5") == 0) {
	mode = '5';
	argc--;
	argv++;
    }

    if (argc >= 2 && strcmp(argv[1], "-d") == 0) {
	mode = 'd';
	argc--;
	argv++;
    }

    if (argc != 2) {
	fprintf(stderr, "sam_pileup filename\n");
	return 1;
    }

    strand_init();

    fp = scram_open(argv[1], "r");
    if (!fp) {
	perror(argv[1]);
	return 1;
    }
    
    if (!(p = calloc(1, sizeof(*p))))
	return 1;

    switch(mode) {
    case '5':
	pileup_loop(fp, NULL, basic_pileup, NULL);
	break;

    case 'd':
	pileup_loop(fp, NULL, depth_pileup, NULL);
	break;

    default:
	pileup_loop(fp, NULL, sam_pileup, p);
	break;
    }

    if (p)
	free(p);

    if (0 != scram_close(fp))
	return 1;

    return 0;
}

