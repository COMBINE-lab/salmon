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
 * Multi-threading.

Trying to multi-thread BAM input is problematic as the BGZF boundaries
share no common ground with the sequence boundaries. Header and
sequences may be larger than a BGZF block, and the header may also be
in the same block as the sequence data.

Therefore the granularity of multi-threading needs to be the
compression and uncompression code.

Threadable sections:

1) zlib portion bgzf_more_output?
2) decode (may merge with 1)
3) bgzf_write


Similary for cram writing:
1) Building crecs
2) Encode + basic huff
3) Block compression (may merge with 2)


Ie input is hard to multi-thread
Output can be distributed to multiple threads and aggregated together
before sending on to FILE *fp.
 */


/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013
 *
 * A wrapper around SAM, BAM and CRAM I/O to give a unified interface.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <string.h>
#include <assert.h>

#include "io_lib/scram.h"

#define SCRAM_BUF_SIZE (1024*1024)

/*
 * Expands the input buffer.
 * Returns 0 on sucess
 *        -1 on failure
 */
static int scram_more_input(scram_fd *fd) {
    size_t avail = fd->alloc - fd->used;
    size_t l;

    l = fread(&fd->buf[fd->used], 1, avail, fd->fp);
    if (l <= 0)
	return -1;
    
    fd->used += l;
    return 0;
}

/*
 * Consumes a block of data from the input stream and returns a malloced
 * copy of it. The input buffer is then copied down (FIXME: inefficient).
 */
static unsigned char *scram_input_next_block(scram_fd *fd, size_t max_size,
					     size_t *out_size) {
    ssize_t l = MIN(max_size, fd->used);
    ssize_t i;
    unsigned char *r = NULL;

    if (max_size > fd->used) {
	scram_more_input(fd);
	if (fd->used == 0)
	    return NULL;
    }

    if (fd->b->binary) {
	uint32_t bsize;

	if (l < 19)
	    return NULL;
	bsize = fd->buf[16] + 256*fd->buf[17] + 1;
	fprintf(stderr, "block_size=%d\n", bsize);
	
	l = MIN(bsize, l);
    } else {
	for (i = l-1; i >= 0; i--) {
	    while (fd->buf[i] != '\n')
		i--;
	}
	assert(i >= 0);

	l = i;
    }

    if (!(r = malloc(l)))
	return NULL;
    memcpy(r, fd->buf, l);
    memcpy(fd->buf, &fd->buf[l], fd->used - l);
    fd->used -= l;

    if (out_size)
	*out_size = l;

    return r;
}

int scram_input_bam_block(scram_fd *fd) {
    size_t sz;
    unsigned char *r;

    if (!fd->is_bam)
	return -1;
    r = scram_input_next_block(fd, Z_BUFF_SIZE, &sz);
    if (!r)
	return -1;

//    if (fd->b->comp_p && fd->b->comp_p != fd->b->comp)
//	free(fd->b->comp_p);
    fd->b->comp_p = r;
    fd->b->comp_sz = sz;
    
    return 0;
}

/*
 * Opens filename.
 * If reading we initially try cram first and then bam/sam if that fails.
 * The exception is when reading from stdin, where bam/sam is first.
 *
 * If writing we look at the mode parameter:
 *     w  => SAM
 *     ws => SAM
 *     wb => BAM
 *     wc => CRAM
 *
 * Returns scram pointer on success
 *         NULL on failure
 */
scram_fd *scram_open(const char *filename, const char *mode) {
    char mode2[10];
    scram_fd *fd = calloc(1, sizeof(*fd));
    if (!fd)
	return NULL;

    fd->eof = 0;

    /* I/O buffer */
    fd->fp = NULL;
    fd->buf = NULL;
    fd->alloc = fd->used = 0;
    fd->pool = NULL;

    if (strcmp(filename, "-") == 0 && mode[0] == 'r'
	&& mode[1] != 'b' && mode[1] != 'c' && mode[1] != 's') { 
	int c;
	/*
	 * Crude auto-detection.
	 * First char @ = sam, 0x1f = bam (gzip), C = cram
	 * Headerless SAM will need explicit mode setting.
	 */
	c = fgetc(stdin);
	ungetc(c, stdin);

	if (c == '@')
	    sprintf(mode2, "rs%.7s", mode+1), mode = mode2;
	else if (c == 0x1f)
	    sprintf(mode2, "rb%.7s", mode+1), mode = mode2;
	else if (c == 'C')
	    sprintf(mode2, "rc%.7s", mode+1), mode = mode2;
    }

    if (*mode == 'r') {
	if (mode[1] != 'b' && mode[1] != 's') {
	    if ((fd->c = cram_open(filename, mode))) {
		cram_load_reference(fd->c, NULL);
		fd->is_bam = 0;
		return fd;
	    }
	}

	if ((fd->b = bam_open(filename, mode))) {
	    fd->is_bam = 1;
	    return fd;
	}
	
	free(fd);
	return NULL;
    }

    /* For writing we cannot auto detect, so create the file type based
     * on the format in the mode string.
     */
    if (strncmp(mode, "wc", 2) == 0) {
	if (!(fd->c = cram_open(filename, mode))) {
	    free(fd);
	    return NULL;
	}
	fd->is_bam = 0;
	return fd;
    }

    /* Otherwise assume bam/sam */
    if (!(fd->b = bam_open(filename, mode))) {
	free(fd);
	return NULL;
    }
    fd->is_bam = 1;
    return fd;
}

#if defined(CRAM_IO_CUSTOM_BUFFERING)
/*
 * Open CRAM file for reading via callbacks
 *
 * Returns scram pointer on success
 *         NULL on failure
 */
scram_fd *scram_open_cram_via_callbacks(
    char const * filename,
    cram_io_allocate_read_input_t   callback_allocate_function,
    cram_io_deallocate_read_input_t callback_deallocate_function,
    size_t const bufsize            
)
{
    scram_fd *fd = calloc(1, sizeof(*fd));
    if (!fd)
	return NULL;

    fd->eof = 0;

    /* I/O buffer */
    fd->fp = NULL;
    fd->buf = NULL;
    fd->alloc = fd->used = 0;
    fd->pool = NULL;

    if ((fd->c = cram_open_by_callbacks(filename,
					callback_allocate_function,
					callback_deallocate_function,
					bufsize))) 
    {
	cram_load_reference(fd->c, NULL);
	fd->is_bam = 0;
	return fd;
    }

    return NULL;
}
#endif

int scram_close(scram_fd *fd) {
    int r;

    if (fd->is_bam) {
	r = bam_close(fd->b);
    } else {
	r = cram_close(fd->c);
    }

    if (fd->pool)
	t_pool_destroy(fd->pool, 0);


    free(fd);
    return r;
}

SAM_hdr *scram_get_header(scram_fd *fd) {
#ifdef __INTEL_COMPILER
    // avoids cmovne generation from icc 2015 (bug)
    return fd->is_bam && fd->b ? fd->b->header : fd->c->header;
#else
    return fd->is_bam ? fd->b->header : fd->c->header;
#endif
}

refs_t *scram_get_refs(scram_fd *fd) {
    return fd->is_bam ? NULL : fd->c->refs;
}

void scram_set_refs(scram_fd *fd, refs_t *refs) {
    if (fd->is_bam)
	return;
    if (fd->c->refs)
	refs_free(fd->c->refs);
    fd->c->refs = refs;
    if (refs)
	refs->count++;
}

void scram_set_header(scram_fd *fd, SAM_hdr *sh) {
    if (fd->is_bam) {
	fd->b->header = sh;
    } else {
	fd->c->header = sh;
    }

    sam_hdr_incr_ref(sh);
}

int scram_write_header(scram_fd *fd) {
    return fd->is_bam
	? bam_write_header(fd->b)
	: cram_write_SAM_hdr(fd->c, fd->c->header);
}

int scram_get_seq(scram_fd *fd, bam_seq_t **bsp) {
    if (fd->is_bam) {
	switch (bam_get_seq(fd->b, bsp)) {
	case 1:
	    return 0;

	case 0:
	    // FIXME: if we ever implement range queries for BAM this will
	    // need amendments to not claim a sub-range is invalid EOF.
	    fd->eof = fd->b->eof_block ? 1 : 2;
	    return -1;

	default:
	    fd->eof = -1; // err
	    return -1;
	}
    }

    if (-1 == cram_get_bam_seq(fd->c, bsp)) {
	fd->eof = cram_eof(fd->c);
	return -1;
    }
    return 0;
}

int scram_next_seq(scram_fd *fd, bam_seq_t **bsp) {
    return scram_get_seq(fd, bsp);
}

int scram_put_seq(scram_fd *fd, bam_seq_t *s) {
    return fd->is_bam
	? bam_put_seq(fd->b, s)
	: cram_put_bam_seq(fd->c, s);
}

int scram_set_option(scram_fd *fd, enum cram_option opt, ...) {
    int r = 0;
    va_list args;

    va_start(args, opt);

    if (opt == CRAM_OPT_THREAD_POOL) {
	t_pool *p = va_arg(args, t_pool *);
	if (fd->is_bam)
	    return bam_set_option(fd->b, BAM_OPT_THREAD_POOL, p);
	else
	    return cram_set_option(fd->c, CRAM_OPT_THREAD_POOL, p);
    } else if (opt == CRAM_OPT_NTHREADS) {
	int nthreads = va_arg(args, int);
	if (nthreads > 1) {
	    if (!(fd->pool = t_pool_init(nthreads*2, nthreads)))
		return -1;

	    if (fd->is_bam)
		return bam_set_option(fd->b, BAM_OPT_THREAD_POOL, fd->pool);
	    else
		return cram_set_option(fd->c, CRAM_OPT_THREAD_POOL, fd->pool);
	} else {
	    fd->pool = NULL;
	    return 0;
	}
    } else if (opt == CRAM_OPT_BINNING) {
	int bin = va_arg(args, int);

	return fd->is_bam
	    ? bam_set_option (fd->b,  BAM_OPT_BINNING, bin)
	    : cram_set_option(fd->c, CRAM_OPT_BINNING, bin);
    }

    if (!fd->is_bam) {
	r = cram_set_voption(fd->c, opt, args);
    }

    va_end(args);

    return r;
}

/*! Returns the line number when processing a SAM file
 *
 * @return
 * Returns line number if input is SAM;
 *         0 for CRAM / BAM input.
 */
int scram_line(scram_fd *fd) {
    if (fd->is_bam)
	return fd->b->line;
    else
	return 0;
}
