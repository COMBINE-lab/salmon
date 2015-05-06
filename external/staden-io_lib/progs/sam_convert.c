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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <fcntl.h>
#include <zlib.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>

#include <io_lib/bam.h>
#include <io_lib/os.h>

int main(int argc, char **argv) {
    bam_file_t *in, *out;
    bam_seq_t *s;
    int nr = 0;
    int nm[2];
    char *mode = "rb";


    if (argc != 3) {
	fprintf(stderr, "Usage: sam_convert input_file output_file\n");
	fprintf(stderr, "Input and output files should end in .sam or .bam.\n");
	return 1;
    }
    
    {
	size_t l = strlen(argv[1]);
	if (l >= 4 && argv[1][l-3] == 's')
	    mode = "r";
    }

    in = bam_open(argv[1], mode);
    if (!in) {
	fprintf(stderr, "Failed to open bam file %s\n", argv[1]);
	return 1;
    }

    {
	size_t l = strlen(argv[2]);
	if (l >= 4 && argv[2][l-3] == 's')
	    mode = "w";
	else
	    mode = "wb";
    }

    out = bam_open(argv[2], mode);
    if (!out) {
	fprintf(stderr, "Failed to open bam file %s\n", argv[2]);
	return 1;
    }

    /* Copy header and refs from in to out, for writing purposes */
    out->header     = in->header;

    if (in->header) {
	if (bam_write_header(out))
	    return 1;
    }

    nm[0] = nm[1] = 0;
    s = NULL;
    while (bam_get_seq(in, &s) > 0) {
	if (-1 == bam_put_seq(out, s))
	    return 1;
	nr++;
	nm[(bam_flag(s) & BAM_FUNMAP) > 0]++;
    }

    printf("Mapped      = %d\n",nm[0]);
    printf("Unmpped     = %d\n",nm[1]);
    printf("Total reads = %d\n",nr);

    out->header = NULL;

    bam_close(in);
    bam_close(out);

    if (s)
	free(s);

    return 0;
}
