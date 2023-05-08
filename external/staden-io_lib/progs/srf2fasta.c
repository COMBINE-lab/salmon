/*
 * Copyright (c) 2008-2010 Genome Research Ltd.
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
#include <math.h>
#include <string.h>
#include <fcntl.h>

#include <io_lib/Read.h>
#include <io_lib/misc.h>
#include <io_lib/ztr.h>
#include <io_lib/srf.h>

/* ------------------------------------------------------------------------ */

#define MAX_READ_LEN 10000
void ztr2fasta(ztr_t *z, char *name) {
    int i, nc;
    char buf[MAX_READ_LEN*2 + 512 + 6];
    char *seq = buf;
    ztr_chunk_t **chunks;

    /* Extract the sequence only */
    chunks = ztr_find_chunks(z, ZTR_TYPE_BASE, &nc);
    if (nc != 1) {
	fprintf(stderr, "Zero or greater than one BASE chunks found.\n");
	if (chunks)
	    free(chunks);
	return;
    }

    uncompress_chunk(z, chunks[0]);

    /* Construct fasta entry */
    *seq++ = '>';
    while (*name)
	*seq++ = *name++;
    *seq++ = '\n';

    for (i = 1; i < chunks[0]->dlength; i++) {
	char base = chunks[0]->data[i];
	if (base == '.')
	    *seq++ = 'N';
	else
	    *seq++ = base;
    }
    *seq++ = '\n';

    fwrite(buf, 1, seq - buf, stdout);
    free(chunks);

    return;
}

/* ------------------------------------------------------------------------ */
void usage(void) {
    fprintf(stderr, "Usage: srf2fasta [-C] archive_name\n");
    exit(1);
}

int main(int argc, char **argv) {
    char *ar_name;
    srf_t *srf;
    char name[512];
    ztr_t *ztr;
    int mask = 0, i;

    /* Parse args */
    for (i = 1; i < argc && argv[i][0] == '-'; i++) {
	if (!strcmp(argv[i], "-")) {
	    break;
	} else if (!strcmp(argv[i], "-C")) {
	    mask = SRF_READ_FLAG_BAD_MASK;
	} else {
	    usage();
	}
    }    

    if (i == argc) {
	usage();
    }
    ar_name = argv[i];

    if (NULL == (srf = srf_open(ar_name, "r"))) {
	perror(ar_name);
	return 4;
    }

    read_sections(READ_BASES);

#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif

    while (NULL != (ztr = srf_next_ztr(srf, name, mask))) {
	ztr2fasta(ztr, name);
	delete_ztr(ztr);
    }

    srf_destroy(srf, 1);

    return 0;
}
