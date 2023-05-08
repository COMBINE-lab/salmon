/*
 * Copyright (c) 2005, 2007, 2010, 2013 Genome Research Ltd.
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
 * Concatenates multiple SFF files together. It also strips out any indexing
 * so this will need to be added again afterwards.
 *
 * The first argument is the archive to append to. All subsequent arguments
 * are the archives to append to the first argument.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <io_lib/sff.h>
#include <io_lib/os.h>

#define BSIZE 1024*1024

int main(int argc, char **argv) {
    int i;
    FILE *fpin, *fpout;
    unsigned char *block;
    unsigned char chdr[31];
    sff_common_header *ch;

    block = (unsigned char *)malloc(BSIZE);

    if (argc < 3) {
	fprintf(stderr, "Usage: append_sff sff_file add_me.sff ...\n");
	return 1;
    }

    /* Open and decode the common header of the archive we'll append to */
    if (NULL == (fpout = fopen(argv[1], "rb+"))) {
	perror(argv[1]);
	return 1;
    }

    if (31 != fread(chdr, 1, 31, fpout)) {
	fprintf(stderr, "Couldn't read common header\n");
	return 1;
    }
    ch = decode_sff_common_header(chdr);

    /*
     * Jump to the end of the archive or the start of the index.
     * NB: If the index is not at the end we cannot use a simple append
     * method.
     */
    fseek(fpout, 0, SEEK_END);
    if (ch->index_len) {
	if (ch->index_offset + ch->index_len != ftell(fpout)) {
	    fprintf(stderr, "Index is not at the end of file => cannot append\n");
	    return 1;
	}
	fseek(fpout, ch->index_offset, SEEK_SET);
    }

    /* Iterate around other archives, appending to the first one */
    for (i = 2; i < argc; i++) {
	sff_common_header *h;
	uint64_t pos;
	size_t len;
	char *sff = argv[i];
	int skipped;

	printf("Copying %s\n", sff);

	if (NULL == (fpin = fopen(sff, "rb"))) {
	    perror(sff);
	    return 1;
	}

	if (31 != fread(chdr, 1, 31, fpin)) {
	    fprintf(stderr, "Couldn't read common header\n");
	    return 1;
	}

	h = decode_sff_common_header(chdr);

	/* Check if headers are compatible */
	if (ch->flow_len        != h->flow_len ||
	    ch->key_len         != h->key_len ||
	    ch->flowgram_format != h->flowgram_format)
	    fprintf(stderr, "*** Error: incompatible SFF headers ***\n");
	fseek(fpin, h->header_len - 31, SEEK_CUR);
	ch->nreads += h->nreads;

	/* Copy all data from fpin to fpout, skipping index if present */
	skipped = h->index_len ? 0 : 1;
	for (pos = ftell(fpin); ; pos += len) {
	    len = BSIZE;
	    if (!skipped) {
		if (pos == h->index_offset) {
		    fseek(fpin, h->index_len, SEEK_CUR);
		    pos += h->index_len;
		    skipped = 1;
		} else if (pos + BSIZE > h->index_offset) {
		    len = h->index_offset - pos;
		}
	    }
	    if (0 == (len = fread(block, 1, len, fpin)))
		break;
	    fwrite(block, 1, len, fpout);
	}

	free_sff_common_header(h);
	fclose(fpin);
    }

    /* Seek back and update the header with the new nreads */
    fseek(fpout, 0, SEEK_SET);
    if (31 != fread(chdr, 1, 31, fpout))
	exit(1);
    *(uint64_t *)(chdr+8)  = be_int8(0);
    *(uint32_t *)(chdr+16) = be_int4(0);
    *(uint32_t *)(chdr+20) = be_int4(ch->nreads);
    fseek(fpout, 0, SEEK_SET);
    fwrite(chdr, 1, 31, fpout);
    fclose(fpout);

    free_sff_common_header(ch);
    free(block);

    return 0;
}
