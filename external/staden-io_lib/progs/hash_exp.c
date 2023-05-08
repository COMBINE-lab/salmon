/*
 * Copyright (c) 2006-2007, 2010-2011 Genome Research Ltd.
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

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <io_lib/hash_table.h>
#include <io_lib/os.h>
#include <io_lib/mFILE.h>

HashFile *build_index(FILE *fp) {
    char line[8192];
    char rname[8192];
    size_t pos = 0, last = 0;
    char *c = "magic";
    HashFile *hf;

    /* Create the hash table */
    hf = HashFileCreate(0, HASH_DYNAMIC_SIZE);
    hf->headers = (HashFileSection *)malloc(sizeof(*hf->headers));
    hf->nheaders = 0;

    for (*rname = 0; c;) {
	c = fgets(line, 8192, fp);

	if (c == NULL || strncmp("ID   ", line, 5) == 0) {
	    /* Add this entry; it extends from 'last' to 'pos' */
	    if (*rname) {
		HashData hd;
		HashFileItem *hfi = (HashFileItem *)calloc(1, sizeof(*hfi));
		
		hfi->header = 0;
		hfi->footer = 0;
		hfi->pos = last;
		hfi->size = pos - last;
		hd.p = hfi;
		
		HashTableAdd(hf->h, rname, strlen(rname), hd, NULL);
	    }

	    /* Remember this ID line for when we meet the next */
	    if (c) {
		char *nl;
		if ((nl = strchr(c, '\n')))
		    *nl = 0;
		if ((nl = strchr(c, '\r')))
		    *nl = 0;

		strcpy(rname, c+5);
	    }

	    last = pos;
	}
	pos = ftell(fp);
    }

    HashTableStats(hf->h, stderr);

    return hf;
}

int main(int argc, char **argv) {
    HashFile *hf;
    FILE *fp;

    if (argc != 2) {
	fprintf(stderr, "Usage: hash_exp exp_file_ball > exp.hash\n");
	return 1;
    }
    if (NULL == (fp = fopen(argv[1], "rb+"))) {
	perror(argv[1]);
	return 1;
    }

    hf = build_index(fp);
    //hf->archive = NULL;

    HashFileSave(hf, fp, 0);

    return 0;
}
