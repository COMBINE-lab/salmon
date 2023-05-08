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

/*
 * Author: James Bonfield, Sanger Institute, 2013.
 *
 * Converts a SAM or BAM file into a CRAM file.
 *
 * Usage:
 *     sam_to_cram [-level] input.sam reference.fasta [output.cram]
 */

#include "io_lib_config.h"

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>

#include <io_lib/cram.h>

void usage(FILE *fp) {
    fprintf(fp, "Usage: sam_to_cram [-r ref.fa] [-0..9] [-u] [-v] [-s int] "
	    "[-S int] in.sam/bam [output.cram]\n\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "    -r ref.fa      Specifies the reference file.\n");
    fprintf(fp, "    -1 to -9       Set zlib compression level for CRAM\n");
    fprintf(fp, "    -0 or -u       No zlib compression.\n");
    fprintf(fp, "    -v             Verbose output.\n");
    fprintf(fp, "    -s integer     Sequences per slice, default %d.\n",
	    SEQS_PER_SLICE);
    fprintf(fp, "    -S integer     Slices per container, default %d.\n",
	    SLICE_PER_CNT);
    fprintf(fp, "    -V version     Specify the CRAM format version to write (eg 1.1, 2.0)\n");
    fprintf(fp, "    -X             Embed reference sequence.\n");
}

int main(int argc, char **argv) {
    cram_fd *out;
    bam_file_t *in;
    bam_seq_t *s = NULL;
    char *out_fn;
    int level = '\0'; // nul terminate string => auto level
    char out_mode[4];
    int c, verbose = 0;
    int s_opt = 0, S_opt = 0, embed_ref = 0;
    char *arg_list, *ref_fn = NULL;

    while ((c = getopt(argc, argv, "u0123456789hvs:S:V:r:X")) != -1) {
	switch (c) {
	case '0': case '1': case '2': case '3': case '4':
	case '5': case '6': case '7': case '8': case '9':
	    level = c;
	    break;
	    
	case 'u':
	    level = '0';
	    break;

	case 'h':
	    usage(stdout);
	    return 0;

	case 'v':
	    verbose++;
	    break;

	case 's':
	    s_opt = atoi(optarg);
	    break;

	case 'S':
	    S_opt = atoi(optarg);
	    break;

	case 'V':
	    cram_set_option(NULL, CRAM_OPT_VERSION, optarg);
	    break;

	case 'r':
	    ref_fn = optarg;
	    break;

	case 'X':
	    embed_ref = 1;
	    break;

	case '?':
	    fprintf(stderr, "Unrecognised option: -%c\n", optopt);
	    usage(stderr);
	    return 1;
	}
    }

    if (argc - optind != 1 && argc - optind != 2) {
	usage(stderr);
	return 1;
    }

    /* opening */
    if (NULL == (in = bam_open(argv[optind], "rb"))) {
	perror(argv[optind]);
	return 1;
    }

    out_fn = argc - optind == 2 ? argv[optind+1] : "-";
    sprintf(out_mode, "wb%c", level);
    if (NULL == (out = cram_open(out_fn, out_mode))) {
	fprintf(stderr, "Error opening CRAM file '%s'.\n", out_fn);
	return 1;
    }

    /* SAM Header */
    if (!(arg_list = stringify_argv(argc, argv)))
	return 1;
    sam_hdr_add_PG(in->header, "sam_to_cram",
		   "VN", PACKAGE_VERSION,
		   "CL", arg_list, NULL);
    free(arg_list);

    /* Find and load reference */
    if (!ref_fn) {
	SAM_hdr_type *ty = sam_hdr_find(in->header, "SQ", NULL, NULL);
	if (ty) {
	    SAM_hdr_tag *tag;

	    if ((tag = sam_hdr_find_key(in->header, ty, "UR", NULL))) {
		ref_fn  = tag->str + 3;
		if (strncmp(ref_fn, "file:", 5) == 0)
		    ref_fn += 5;
	    }
	}
    }

    out->header = in->header;
    if (ref_fn)
	cram_load_reference(out, ref_fn);

    if (!out->refs) {
	fprintf(stderr, "Unable to open reference.\n"
		"Please specify a valid reference with -r ref.fa option.\n");
	return 1;
    }
    refs2id(out->refs, out->header);

    if (-1 == cram_write_SAM_hdr(out, in->header))
	return 1;

    cram_set_option(out, CRAM_OPT_VERBOSITY, verbose);
    if (s_opt)
	cram_set_option(out, CRAM_OPT_SEQS_PER_SLICE, s_opt);

    if (S_opt)
	cram_set_option(out, CRAM_OPT_SLICES_PER_CONTAINER, S_opt);

    if (embed_ref)
	cram_set_option(out, CRAM_OPT_EMBED_REF, embed_ref);

    /* Sequence iterators */
    while (bam_get_seq(in, &s) > 0) {
	if (-1 == cram_put_bam_seq(out, s)) {
	    fprintf(stderr, "Failed in cram_put_bam_seq()\n");
	    return 1;
	}
    }

    bam_close(in);
    out->header = NULL; // freed by bam_close()
    if (-1 == cram_close(out)) {
	fprintf(stderr, "Failed in cram_close()\n");
	return 1;
    }

    if (s)
	free(s);

    return 0;
}
