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

/*! \file
 * Generic SAM/BAM/CRAM interface.
 *
 * This file implements a higher level scram_*() API for programs that
 * wish to be file format agnostic.
 */

#ifndef _SCRAM_H_
#define _SCRAM_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include "io_lib/bam.h"
#include "io_lib/cram.h"

/*! The primary file handle for reading and writing. */
typedef struct {
    int is_bam;
    int eof;
    union {
	bam_file_t *b;
	cram_fd    *c;
    };

    /* Primary Input/Output buffer */
    unsigned char *buf;
    size_t alloc;
    size_t used;
    FILE *fp;   // copy of file handle.

    t_pool *pool;
} scram_fd;

/*
 * An input stream in SCRAM is a large block of memory which we periodically
 * fread into.
 *
 * This input stream is then broken down into chunks of appropriate size
 * as used by the underlying format. The only tricky bit here is the first
 * portion (opening the underlying format) can use an unknown amount of 
 * buffer due to the BAM header being variable length.
 *
 * Once we have this, scram_next_input() will return the next natural
 * chunk from the input buffer. This permits a single input buffer being
 * divided into multiple scram_buffers to pass to separate threads for
 * decoding.
 */
typedef struct {
    unsigned char *buf;
    size_t alloc; // allocated size of buf
    size_t size;  // size loaded
    size_t usize; // size usable by the underlying format
} scram_buffer_t;

/*!@return
 * Returns 0 if not at end of file
 *         1 if we hit an expected EOF (end of range or EOF block)
 *         2 for other EOF (end of stream without EOF block)
 */
#define scram_eof(fd) ((fd)->eof)


/*! Opens a file.
 *
 * If reading we look for the following mode parameters:
 * -    r  => Try SAM/BAM first, if fail try CRAM
 * -    rb => BAM
 * -    rc => CRAM
 *
 * If writing we look at the mode parameter:
 * -    w  => SAM
 * -    wb => BAM
 * -    wc => CRAM
 *
 * Additionally we can specify the compression level when writing
 * after the file type character, as 0 to 9. Eg "wb9" for maximum
 * compression of BAM or "wc0" for uncompressed CRAM.
 *
 * @return
 * Returns scram pointer on success
 *         NULL on failure
 */
scram_fd *scram_open(const char *filename, const char *mode);

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
);
#endif

/*! Closes a scram_fd handle
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int scram_close(scram_fd *fd);


/*! Returns the SAM_hdr struct.
 *
 * @return
 * The SAM_hdr struct on success; NULL on failure.
 */
SAM_hdr *scram_get_header(scram_fd *fd);


/*! Sets the SAM_hdr struct.
 *
 * Note that this sets the raw pointer and does not take an internal
 * copy of it. If you need to do this call sam_hdr_dup() first.
 */
void scram_set_header(scram_fd *fd, SAM_hdr *sh);


/*! Writes the SAM hdr.
 *
 * This calls the appropriate SAM, BAM or CRAM I/O function to write
 * out the SAM_hdr currently associated with this fd.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int scram_write_header(scram_fd *fd);


/*! Returns the reference sequence array.
 *
 * Note: this only works for CRAM files.
 *
 * @return
 * Returns the refs structure on success;
 *         NULL on failure.
 * 
 * After failure, check with scram_eof(fd) to see whether an genuine
 * error occurred or whether we hit the end of file.
 */
refs_t *scram_get_refs(scram_fd *fd);


/*! Sets the reference sequence array.
 *
 * Note: this only works for CRAM files.
 */
void scram_set_refs(scram_fd *fd, refs_t *refs);


/*!
 * Replaces the FILE* input interface with an explicit buffer to decode
 * from.
 *
 * @Returns 0 on success;
 *         -1 on failure
 */
int scram_input_buffer(scram_fd *fd, unsigned char *buf, size_t size);


/*! Fetches the next sequence and returns it in BAM format.
 *
 * This reads a new sequence line from fd and returns it in the BAM
 * in-memory format, regardless of whether the input file was SAM, BAM
 * or CRAM.
 *
 * @param bsp bsp is a pointer to a bam_seq_t*, as our usual bam_seq_t
 * structure pointer may be reallocated internally by this
 * function. It is permitted to pass in the address of a bam_seq_t*
 * that points to NULL. This behaviour differs to the Samtools API due
 * to the bam_seq_t structure being a single contiguous block of
 * memory instead of in two halves; the static and variable "data"
 * component.
 *
 * Note: For maximum speed of CRAM I/O you may wish to use the cram
 * specific layer and return cram_record objects instead.
 *
 * @return
 * Returns 0 on success and fills out bsp;
 *        -1 on failure
 */
int scram_get_seq(scram_fd *fd, bam_seq_t **bsp);

/*! Deprecated: please use scram_get_seq() instead */
int scram_next_seq(scram_fd *fd, bam_seq_t **bsp);


/*! Writes a BAM encoded bam_seq_t to fd.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int scram_put_seq(scram_fd *fd, bam_seq_t *s);


/*! Sets a CRAM option on fd.
 *
 * This is only supported for CRAM files currently.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int scram_set_option(scram_fd *fd, enum cram_option opt, ...);

/*! Returns the line number when processing a SAM file
 *
 * @return
 * Returns line number if input is SAM;
 *         0 for CRAM / BAM input.
 */
int scram_line(scram_fd *fd);
#ifdef __cplusplus
}
#endif

#endif /* _SCRAM_H_ */
