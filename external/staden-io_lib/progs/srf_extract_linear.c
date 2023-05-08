/*
 * Copyright (c) 2007-2008 Genome Research Ltd.
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
 * This performs a linear (non-indexed) search for a trace in an SRF archive.
 *
 * It's not intended as a suitable production program or as a library of code
 * to use, but as a test and benchmark statistic.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <io_lib/Read.h>
#include <io_lib/misc.h>
#include <io_lib/ztr.h>
#include <io_lib/srf.h>

/* ------------------------------------------------------------------------ */
/*
 * Looks for a trace name in an SRF archive and returns the binary contents
 * if found, or NULL if not.
 */
mFILE *find_reading(srf_t *srf, char *tr_name) {
    do {
	int type;

	switch(type = srf_next_block_type(srf)) {
	case -1:
	    /* EOF */
	    return NULL;

	case SRFB_CONTAINER:
	    if (0 != srf_read_cont_hdr(srf, &srf->ch))
		return NULL;
	    break;

	case SRFB_XML:
	    if (0 != srf_read_xml(srf, &srf->xml))
		return NULL;
	    break;
	    

	case SRFB_TRACE_HEADER: {
	    /* off_t pos = ftell(srf->fp); */

	    if (0 != srf_read_trace_hdr(srf, &srf->th))
		return NULL;

#if 0
	    /*
	     * If the name prefix doesn't match tr_name then skip this entire
	     * block.
	     */
	    if (0 != strncmp(tr_name, srf->th.id_prefix,
			     strlen(srf->th.id_prefix)) &&
		0 != srf->th.next_block_offset) {
		fseek(srf->fp, pos + srf->th.next_block_offset, SEEK_SET);
	    }
#endif
	    break;
	}

	case SRFB_TRACE_BODY: {
	    mFILE *mf = mfcreate(NULL, 0);
	    srf_trace_body_t tb;
	    char name[512];

	    if (!mf || 0 != srf_read_trace_body(srf, &tb, 0))
		return NULL;

	    sprintf(name, "%s%s", srf->th.id_prefix, tb.read_id);
	    if (strcmp(name, tr_name)) {
		mfdestroy(mf);
		if (tb.trace)
		    free(tb.trace);
		continue;
	    }

	    if (srf->th.trace_hdr_size)
		mfwrite(srf->th.trace_hdr, 1,
			srf->th.trace_hdr_size, mf);
	    if (tb.trace_size)
		mfwrite(tb.trace, 1, tb.trace_size, mf);
	    if (tb.trace)
		free(tb.trace);
	    mrewind(mf);
	    return mf;
	}

	case SRFB_INDEX: {
	    off_t pos = ftello(srf->fp);
	    srf_read_index_hdr(srf, &srf->hdr, 1);

	    /* Skip the index body */
	    fseeko(srf->fp, pos + srf->hdr.size, SEEK_SET);
	    break;
	}

 	case SRFB_NULL_INDEX:
 	    break;

	default:
	    fprintf(stderr, "Block of unknown type '%c'. Aborting\n", type);
	    return NULL;
	}
    } while (1);

    return NULL;
}

/* ------------------------------------------------------------------------ */
int main(int argc, char **argv) {
    char *ar_name, *tr_name;
    mFILE *mf;
    srf_t *srf;

    if (argc != 3) {
	fprintf(stderr, "Usage: extract_linear_srf archive_name trace_name\n");
	return 1;
    }
    ar_name = argv[1];
    tr_name = argv[2];

    if (NULL == (srf = srf_open(ar_name, "r"))) {
	perror(ar_name);
	return 4;
    }

    if (NULL == (mf = find_reading(srf, tr_name))) {
	fprintf(stderr, "%s not found in archive\n", tr_name);
	return 3;
    }

#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif

    fwrite(mf->data, 1, mf->size, stdout);
    return 0;
}
