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
#include <fcntl.h>
#include <io_lib/hash_table.h>

/*
 * Copies a single named file to stdout.
 * Returns 0 on success
 *         1 on failure
 */
int extract(HashFile *hf, char *file) {
    size_t len;
    char *data;

    if ((data = HashFileExtract(hf, file, &len))) {
	fwrite(data, len, 1, stdout);
	free(data);
	return 0;
    }
    return 1;
}

int main(int argc, char **argv) {
    char *fofn = NULL;
    char *hash;
    HashFile *hf;
    int ret = 0;

    /* process command line arguments of the form -arg */
    for (argc--, argv++; argc > 0; argc--, argv++) {
	if (**argv != '-' || strcmp(*argv, "--") == 0)
	    break;

	if (strcmp(*argv, "-I") == 0) {
	    argv++;
	    fofn = *argv;
	    argc--;
	}
    }

    if (argc < 2 && !fofn) {
	fprintf(stderr, "Usage: hash_extract [-I fofn] hashfile [name ...]\n");
	return 1;
    }
    hash = argv[0];
    argc--;
    argv++;

    if (NULL == (hf = HashFileOpen(hash))) {
	perror(hash);
	return 1;
    }

    if (fofn) {
	FILE *fofnfp;
	char file[256];

	if (strcmp(fofn, "-") == 0) {
	    fofnfp = stdin;
	} else {
	    if (NULL == (fofnfp = fopen(fofn, "r"))) {
		perror(fofn);
		return 1;
	    }
	}

	while (fgets(file, 255, fofnfp)) {
	    char *c;
	    if ((c = strchr(file, '\n')))
		*c = 0;

	    ret |= extract(hf, file);
	}

	fclose(fofnfp);
    }

#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    for (; argc; argc--, argv++) {
	ret |= extract(hf, *argv);
    }

    HashFileDestroy(hf);

    return ret;
}
