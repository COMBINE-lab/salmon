/*
 * Copyright (c) 2007-2008, 2013 Genome Research Ltd.
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
 * This adds a hash table index (".hsh" v1.01 format) to an SRF archive.
 * It does this either inline on the file itself (provided it doesn't already
 * have an index) or by producing an external index file.
 */

/* ---------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <io_lib/hash_table.h>
#include <io_lib/os.h>
#include <io_lib/array.h>
#include <io_lib/srf.h>

/* ------------------------------------------------------------------------ */
void usage(int code) {
    printf("Usage: srf_index_hash [-c] srf_file\n");
    printf(" Options:\n");
    printf("    -c       check an existing index, don't re-index\n");
    exit(code);
}

int main(int argc, char **argv) {
    srf_t *srf;
    uint64_t pos;
    char name[512];
    int i, type;
    char *archive;
    int dbh_pos_stored_sep = 0;
    int check = 0;
    off_t old_index = 0;
    srf_index_t *idx;
    
    /* Parse args */
    for (i = 1; i < argc && argv[i][0] == '-'; i++) {
	if (!strcmp(argv[i], "-")) {
	    break;
	} else if (!strcmp(argv[i], "-c")) {
	    check = 1;
	} else if (!strcmp(argv[i], "-h")) {
	    usage(0);
	} else {
	    usage(1);
	}
    }

    if (argc != (i+1)) {
      usage(1);
    }

    archive = argv[i];

    if( check ){
        srf = srf_open(archive, "rb");
    } else {
        srf = srf_open(archive, "r+b");
    }
    if (NULL == srf ){
 	perror(argv[i]);
	return 1;
    }
    
    
    idx = srf_index_create(NULL, NULL, dbh_pos_stored_sep);
    if (NULL == idx)
	return 1;

    /* Scan through file gathering the details to index in memory */
    while ((type = srf_next_block_details(srf, &pos, name)) >= 0) {
	/* Only want this set if the last block in the file is an index */
	old_index = 0;

	switch (type) {
	case SRFB_CONTAINER:
	    if (srf_index_add_cont_hdr(idx, pos))
		return 1;
	    break;

	case SRFB_TRACE_HEADER:
	    if (srf_index_add_trace_hdr(idx, pos))
		return 1;
	    break;

	case SRFB_TRACE_BODY:
	    if (srf_index_add_trace_body(idx, name, pos))
		return 1;
	    break;

	case SRFB_INDEX:
	    /* An old index */
	    old_index = pos;
	    break;

        case SRFB_NULL_INDEX:
	    old_index = pos;
            break;

	default:
	    abort();
	}
    }

    /* the type should be -1 (EOF) */
    if( type != -1 )
        abort();

    /* are we really at the end of the srf file */
    pos = ftell(srf->fp);
    fseek(srf->fp, 0, SEEK_END);
    if( pos != ftell(srf->fp) ){
        fprintf(stderr, "srf file is corrupt\n");
	return 1;
    }
    
    if (check) {
	srf_index_destroy(idx);
	srf_destroy(srf, 1);
	return 0;
    }

    /* Write out the index */
    if (old_index)
	fseeko(srf->fp, old_index, SEEK_SET);

    srf_index_stats(idx, NULL);
    srf_index_write(srf, idx);

    /* Truncate incase we've somehow overwritten an old longer index */
    if (ftello(srf->fp) != -1 &&
	ftruncate(fileno(srf->fp), ftello(srf->fp)) == -1)
	return -1;

    srf_index_destroy(idx);
    srf_destroy(srf, 1);
    
    return 0;
}
