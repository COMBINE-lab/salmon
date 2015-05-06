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

#ifndef _CRAM_INDEX_H_
#define _CRAM_INDEX_H_

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Loads a CRAM .crai index into memory.
 * Returns 0 for success
 *        -1 for failure
 */
int cram_index_load(cram_fd *fd, char const * fn);

#if defined(CRAM_IO_CUSTOM_BUFFERING)
/*
 * Loads a CRAM .crai index into memory using callbacks. fn denotes the name of the cram file.
 * The suffix .crai will be appended to obtain the name of the index.
 *
 * Returns 0 for success
 *        -1 for failure
 */
extern int cram_index_load_via_callbacks(
    cram_fd *fd, char const * fn,
    cram_io_allocate_read_input_t   callback_allocate_function,
    cram_io_deallocate_read_input_t callback_deallocate_function        
);
#endif

void cram_index_free(cram_fd *fd);

/*
 * Searches the index for the first slice overlapping a reference ID
 * and position.
 *
 * Returns the cram_index pointer on sucess
 *         NULL on failure
 */
cram_index *cram_index_query(cram_fd *fd, int refid, int pos, cram_index *frm);

/*
 * Skips to a container overlapping the start coordinate listed in
 * cram_range.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_seek_to_refpos(cram_fd *fd, cram_range *r);

/*
 * Seek within a cram file.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_seek(cram_fd *fd, off_t offset, int whence);

/*
 * Builds an index file.
 *
 * fd is a newly opened cram file that we wish to index.
 * fn_base is the filename of the associated CRAM file. Internally we
 * add ".crai" to this to get the index filename.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_index_build(cram_fd *fd, const char *fn_base);

#ifdef __cplusplus
}
#endif

#endif
