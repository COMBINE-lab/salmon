/*
 * Copyright (c) 2015 Genome Research Ltd.
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
 * This code contains an interface for Biobambam and Libmaus to use
 * when encoding CRAM files.
 *
 * Libmaus has its own thread-pool system.  We construct work packages
 * to be dispatched to libmaus.
 *
 * The basic model is that we have a single encoder context, created /
 * destroyed by cram_allocate_encoder() and cram_deallocate_encoder().
 * The context can just be a cram_fd pointer.
 *
 * This context is then passed along with a block of uncompressed BAM
 * records to cram_enque_compression_block(), running in the top thread.
 * This procedure creates a work package and puts it on the work queue
 * passed in. 
 *
 * The cram_process_work_package() function is then called per work
 * package, within 1 or more threads.  Within this function it will
 * call the libmaus supplied write function and when finished will
 * call the libmaus supplied work-finished function.  The write
 * function takes the place of the CRAM_IO_PUTC and CRAM_IO_WRITE
 * functions already defined in cram_io.c.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "io_lib/cram.h"


//-----------------------------------------------------------------------------
// FIXME: move these to a public header.
// Public interface

typedef enum cram_data_write_block_type {
    cram_data_write_block_type_internal,
    cram_data_write_block_type_block_final,
    cram_data_write_block_type_file_final
} cram_data_write_block_type;

// Enqueue a package of work to compress a CRAM slice.
typedef void
(*cram_enque_compression_work_package_function_t)(void *userdata,
						  void *workpackage);

// Callback to indicate the block has been compressed
typedef void
(*cram_compression_work_package_finished_t)(void *userdata,
					    size_t const inblockid,
					    int const final);

// Write function for compressed blocks, provided by libmaus.
// Inblockid is the same as supplied by the dispatcher.
// Outblockid increments from 0 per unique inblockid.
typedef void
(*cram_data_write_function_t)(void *userdata,
			      ssize_t const inblockid,
			      size_t const outblockid,
			      char const *data,
			      size_t const n,
			      cram_data_write_block_type const blocktype);


// Temporary copy from biobambam (BSD licence verbally accepted) to
// help validate the interface via the compiler.
extern void *cram_allocate_encoder(void *userdata,
				   char const *sam_header,
				   size_t const sam_headerlength,
				   cram_data_write_function_t write_func);
extern void cram_deallocate_encoder(void *context);
extern int cram_enque_compression_block(
	void *userdata,
	void *context,
	size_t const inblockid,
	char const **block,
	size_t const *blocksize,
	size_t const *blockelements,
	size_t const numblocks,
	int const final,
	cram_enque_compression_work_package_function_t workenqueuefunction,
	cram_data_write_function_t writefunction,
	cram_compression_work_package_finished_t workfinishedfunction);
extern int cram_process_work_package(void *workpackage);


//-----------------------------------------------------------------------------
// Internally used structures


// The cram_enc_context is the primary encoder context used by
// libmaus.  It is allocated once and passed into all subsequent calls
// by libmaus.  We use it to track the total number of records so far
// (for the CRAM container header) and our CRAM in-memory file
// descriptor.
typedef struct {
    cram_fd *fd;
    void *userdata; // supplied by caller, pass back into write_func
    cram_data_write_function_t write_func;
    size_t num_records;
    pthread_mutex_t context_lock; // a lock for manipulating this struct
    pthread_mutex_t header_lock;  // a lock on fd->header
} cram_enc_context;

// A work package is a series of BAM blocks for conversion to CRAM
typedef struct {
    // encoder context
    cram_enc_context *context;

    // BAM block
    char const **block;
    size_t const *blocksize;
    size_t num_blocks;
    ssize_t inblockid;
    size_t outblockid;
    int final;

    size_t num_records;

    // Libmaus callback functions
    void *userdata;
    cram_data_write_function_t write_func;
    cram_compression_work_package_finished_t finished_func;
} cram_enc_work_package;


//-----------------------------------------------------------------------------
// The callbacks for making CRAM write to an in-memory data block.

/*
 * Cram buffered I/O writer.
 * Has an fwrite()-like interface.  It returns the number of items
 * (not bytes) written.  As dstring_nappend either writes everything
 * or nothing, this will return nmemb on success and 0 on failure.
 */
size_t cram_mem_write_callback(void *ptr,
			       size_t size,
			       size_t nmemb,
			       void *cram_io_userdata) {
    dstring_t *ds = (dstring_t *)cram_io_userdata;

    if (0 == dstring_nappend(ds, ptr, size*nmemb))
	return nmemb;

    return 0;
}

static cram_io_output_t *
cram_callback_allocate_func(char const *filename) {
    cram_io_output_t *io = malloc(sizeof(*io));
    dstring_t *ds = dstring_create(NULL);

    if (!io)
	return NULL;

    io->user_data = ds;
    io->fwrite_callback = cram_mem_write_callback;
    //io->ftell_callback = cram_mem_tell_callback;
    io->ftell_callback = NULL;

    return io;
}

static cram_io_output_t *
cram_callback_deallocate_func(cram_io_output_t *io) {
    if (io) {
	if (io->user_data)
	    dstring_destroy(io->user_data);
	free(io);
    }

    return NULL;
}

//-----------------------------------------------------------------------------
// The libmaus threading interface itself.

/**
 * Allocate cram encoder and return compression context
 *
 * @param samheader SAM header text
 * @param samheaderlength length of SAM header in bytes
 * @param workenqueuefunction function which will be called to enque
 *        compression work packages
 *
 * @return NULL on failure.
 **/
void *cram_allocate_encoder(void *userdata,
			    char const *sam_header,
			    size_t const sam_headerlength,
			    cram_data_write_function_t write_func) {
    cram_fd *fd = NULL;
    SAM_hdr *hdr = NULL;
    cram_enc_context *c = malloc(sizeof(*c));

    if (!c)
	goto err;

    if (!(hdr = sam_hdr_parse(sam_header, sam_headerlength)))
	goto err;

    fd = cram_openw_by_callbacks(NULL,
				 cram_callback_allocate_func,
				 cram_callback_deallocate_func,
				 1024*1024);
    if (!fd)
	goto err;

    //fd->inblockid = 0;
    //fd->outblockid = 0;

    fd->header = hdr;
    sam_hdr_incr_ref(hdr);
    if (cram_write_SAM_hdr(fd, hdr) != 0)
	goto err;

    cram_io_flush_output_buffer(fd);

    c->fd = fd;
    c->userdata = userdata;
    c->write_func = write_func;
    c->num_records = 0;

    // While the cram_fd itself does not have its own internal
    // multithreading, we manually create the mutexes it would use
    // to ensure that this cram_fd running in our own separate
    // threads can handle locking correctly.
    fd->metrics_lock = malloc(sizeof(pthread_mutex_t));
    fd->ref_lock = malloc(sizeof(pthread_mutex_t));
    fd->bam_list_lock = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(fd->metrics_lock, NULL);
    pthread_mutex_init(fd->ref_lock, NULL);
    pthread_mutex_init(fd->bam_list_lock, NULL);

    dstring_t *ds = (dstring_t *)fd->fp_out_callbacks->user_data;
    write_func(userdata, -1, 0,
	       DSTRING_STR(ds), DSTRING_LEN(ds),
	       cram_data_write_block_type_block_final);

    pthread_mutex_init(&c->context_lock, NULL);
    pthread_mutex_init(&c->header_lock, NULL);

    return c;

 err:
    if (c)
	free(c);

    if (fd)
	cram_close(fd);

    if (hdr)
	sam_hdr_free(hdr);

    return NULL;
}

void cram_deallocate_encoder(void *context) {
    cram_enc_context *c = (cram_enc_context *)context;
    cram_fd *fd;

    if (!c)
	return;

    fd = c->fd;

    pthread_mutex_destroy(&c->context_lock);
    pthread_mutex_destroy(&c->header_lock);

    pthread_mutex_destroy(fd->metrics_lock);
    pthread_mutex_destroy(fd->ref_lock);
    pthread_mutex_destroy(fd->bam_list_lock);
    free(fd->metrics_lock);
    free(fd->ref_lock);
    free(fd->bam_list_lock);

    if (fd->header)
	sam_hdr_free(fd->header);

    if (fd)
	cram_close(fd);

    free(c);
}


/**
 * Notify cram encoder there is more data to be compressed
 *
 * @param userdata pointer supplied back to callback functions
 * @param context compression context returned by
 *        cram_allocate_encoder
 * @param inblockid running id of input block
 * @param block of alignments (uncompressed bam format)
 * @param blocksize length of block in sizeof(char)
 * @param final 1 if this is the last block passed, 0 otherwise
 * @param workenqueuefunction callback for queueing work in the thread
 *        pool
 * @param workfinishedfunction callback notifying caller that
 *        compression of this block is done
 * @param writefunction function called from writing compressed data
 *
 * @return 0 on success;
 *        -1 on failure
 **/
int cram_enque_compression_block(
	void *userdata,
	void *context,
	size_t const inblockid,
	char const **block,
	size_t const *blocksize,
	size_t const *blockelements,
	size_t const numblocks,
	int const final,
	cram_enque_compression_work_package_function_t workenqueuefunction,
	cram_data_write_function_t writefunction,
	cram_compression_work_package_finished_t workfinishedfunction)
{
    cram_enc_context *c = (cram_enc_context *)context;
    cram_enc_work_package *pkg = malloc(sizeof(*pkg));
    size_t n, numrecs;

    if (!pkg)
	return -1;

    numrecs = 0;
    for (n = 0; n < numblocks; n++)
	numrecs += blockelements[n];

    pthread_mutex_lock(&c->context_lock);
    pkg->num_records = c->num_records;
    c->num_records += numrecs;
    pthread_mutex_unlock(&c->context_lock);

#if defined(IO_LIB_CRAM_BAMBAM_DEBUG)
    {
	fprintf(stderr, "Enqueue block %d, rec_start %d+%d, final %d\n",
		(int)inblockid, (int)pkg->num_records, (int)numrecs, final);
	fprintf(stderr, "blocksize[]={");
	int tot;
	for (tot = n = 0; n < numblocks; n++) {
	    fprintf(stderr, "%d%c", 
		    (int)blocksize[n],
		    "},"[n < numblocks-1]);
	    tot += blocksize[n];
	}
	fprintf(stderr, "; // sum %d\n", tot);
    }
#endif

    pkg->context       = c;
    pkg->block         = block;
    pkg->blocksize     = blocksize;
    pkg->num_blocks    = numblocks;
    pkg->inblockid     = inblockid;
    pkg->outblockid    = 0;
    pkg->final         = final;
    pkg->userdata      = userdata;
    pkg->write_func    = writefunction;
    pkg->finished_func = workfinishedfunction;

    workenqueuefunction(userdata, pkg); // FIXME: should probably have
					// err code.

    return 0;
}

static cram_fd *cram_dup_fd(cram_fd *orig) {
    int bufsize = 65536; // FIXME
    cram_fd *fd = malloc(sizeof(*fd));

    if (!fd)
	return NULL;

    memcpy(fd, orig, sizeof(*fd));
    fd->ctr = NULL;

    fd->fp_out_buffer = cram_io_allocate_output_buffer(bufsize);
    fd->fp_out_callbacks = cram_callback_allocate_func(NULL);

    fd->fp_out = NULL;

    return fd;
}

static void cram_dup_close(cram_fd *fd) {
    spare_bams *bl, *next;

    if (!fd)
	return;

    if (fd->fp_out_buffer)
	cram_io_deallocate_output_buffer(fd->fp_out_buffer);

    if (fd->fp_out_callbacks)
	cram_callback_deallocate_func(fd->fp_out_callbacks);

    for (bl = fd->bl; bl; bl = next) {
	int i, max_rec = fd->seqs_per_slice * fd->slices_per_container;

	next = bl->next;
	for (i = 0; i < max_rec; i++) {
	    if (bl->bams[i])
		free(bl->bams[i]);
	}
	free(bl->bams);
	free(bl);
    }

    if (fd->ctr)
	cram_free_container(fd->ctr);

    free(fd);
}


/**
 * Work package dispatch function for cram encoding
 *
 * @param Workpackage containing task to perform and all function
 *        pointers necessary to communicate back to dispatcher.
 *
 * @return 0 on success;
 *        -1 on failure
 **/
int cram_process_work_package(void *workpackage) {
    cram_enc_work_package *pkg = (cram_enc_work_package *)workpackage;
    cram_enc_context *c;
    cram_fd *fd;
    size_t bnum;

    if (!pkg)
	return -1;

    if (!(c = pkg->context))
	return -1;

    // Each work package can be running in a separate thread, so we
    // need to make sure writing to CRAM isn't clobbering over shared
    // memory.
    //
    // The reference sequences work fine with reference counting, but
    // the output buffer is one per cram_fd.  Therefore we create a
    // temporary local copy of cram_fd with pointers to share as much
    // as we can.
    //
    // FIXME: consider having a free-list of previously used cram_fd.
    pthread_mutex_lock(&c->context_lock);
    fd = cram_dup_fd(c->fd);
    pthread_mutex_unlock(&c->context_lock);

    fd->record_counter = pkg->num_records;


    // We create a fake bam_file_t containing the entire BAM block and
    // then use the standard bam_get_seq() API to iterate over
    // sequences within the BAM block.
    for (bnum = 0; bnum < pkg->num_blocks; bnum++) {
	bam_file_t *bf;
	bam_seq_t *bsp = NULL;

	pthread_mutex_lock(&c->header_lock);
	bf = bam_open_block(pkg->block[bnum],
			    pkg->blocksize[bnum],
			    fd->header);
	pthread_mutex_unlock(&c->header_lock);
	if (!bf)
	    return -1;

	while (bam_get_seq(bf, &bsp)) {
	    if (cram_put_bam_seq(fd, bsp) != 0) {
		fprintf(stderr, "Failed to write CRAM record\n");
		pthread_mutex_lock(&c->header_lock);
		bam_close(bf);
		pthread_mutex_unlock(&c->header_lock);

		cram_dup_close(fd);
		return -1;
	    }
	}

	pthread_mutex_lock(&c->header_lock);
	bam_close(bf);
	pthread_mutex_unlock(&c->header_lock);

	if (bsp)
	    free(bsp);
    }

    cram_flush(fd);

    if (pkg->final) {
	// The final package needs the EOF block adding too.
	cram_write_eof_block(fd);
    }

    // Write the block
    dstring_t *ds = (dstring_t *)fd->fp_out_callbacks->user_data;
#if defined(IO_LIB_CRAM_BAMBAM_DEBUG)
    fprintf(stderr, "Writing work package %d,%d "
	    "from rec %d, length %d, final %d\n",
	    (int)pkg->inblockid, (int)pkg->outblockid,
	    (int)pkg->num_records,
	    (int)DSTRING_LEN(ds),
	    pkg->final);
#endif

    pkg->write_func(pkg->userdata, 
		    pkg->inblockid,
		    pkg->outblockid++,
		    DSTRING_STR(ds),
		    DSTRING_LEN(ds),
		    pkg->final
		    ? cram_data_write_block_type_file_final
		    : cram_data_write_block_type_block_final);

    pkg->finished_func(pkg->userdata, pkg->inblockid, pkg->final);

    // Free the work package
    free(pkg);

    // FIXME: do we also need to do something to decr reference seqs?
    cram_dup_close(fd);

    return 0;
}

