/*
 * Copyright (c) 2005-2008, 2010, 2013 Genome Research Ltd.
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
#include <stdlib.h>

#include <io_lib/hash_table.h>
#include <io_lib/os.h>

/*
 * Dumps a textual represenation of the hash table to stdout.
 */
void HashTableLongDump(HashFile *hf, FILE *fp, int long_format) {
    HashTable *h = hf->h;
    int i;

    for (i = 0; i < h->nbuckets; i++) {
	HashItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    HashFileItem *hfi;
	    hfi = (HashFileItem *)hi->data.p;
	    if (long_format) {
		char *aname;
		
		if (hf->archives &&
		    hfi->archive < hf->narchives) {
		    aname = hf->archives[hfi->archive];
		} else {
		    aname = "?";
		}
		fprintf(fp, "%10"PRId64" %6"PRId32" %.*s %s\n",
			hfi->pos, hfi->size, hi->key_len, hi->key,
			aname);
		/*
		fprintf(fp, "%10ld %6d %.*s\n",
			hfi->pos, hfi->size, hi->key_len, hi->key);
		*/
	    } else {
		fprintf(fp, "%.*s\n", hi->key_len, hi->key);
	    }
	}
    }
}


/*
 * Lists the contents of a .hash file
 */
int main(int argc, char **argv) {
    FILE *fp;
    HashFile *hf;
    int long_format = 0;

    /* process command line arguments of the form -arg */
    if (argc >= 2 && strcmp(argv[1], "-l") == 0) {
	long_format = 1;
	argc--;
	argv++;
    }
    if (argc >= 2) {
	fp = fopen(argv[1], "rb");
	if (NULL == fp) {
	    perror(argv[1]);
	    return 1;
	}
    } else {
	fp = stdin;
    }

    hf = HashFileLoad(fp);
    if (hf) {
	HashTableLongDump(hf, stdout, long_format);
	HashFileDestroy(hf);
    }

    return 0;
}
