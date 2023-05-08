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
 * Converts a CRAM file to a SAM or BAM file.
 */

#include "io_lib_config.h"

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>

#include <io_lib/cram.h>
#include <io_lib/bam.h>

void usage(FILE *fp) {
    fprintf(fp, "Usage: cram_to_sam [-r ref.fa] [-m] [-b] [-0..9] [-u] "
	    "filename.cram [output_filename]\n\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "    -r ref.fa      Specifies the reference file.\n");
    fprintf(fp, "    -m             Generate MD and NM tags:\n");
    fprintf(fp, "    -b             Output in BAM (defaults to SAM)\n");
    fprintf(fp, "    -1 to -9       Set zlib compression level for BAM\n");
    fprintf(fp, "    -0 or -u       Output uncompressed, if BAM.\n");
    fprintf(fp, "    -p str         Set the prefix for auto-generated seq. names\n");
    fprintf(fp, "    -R region	    Extract region 'ref:start-end', eg -R chr1:1000-2000\n");
    fprintf(fp, "    -X             Extract using the embedded reference (if present).\n");
}

int main(int argc, char **argv) {
    cram_fd *fd;
    bam_file_t *bfd;
    bam_seq_t *bam = NULL;
    char mode[4] = {'w', '\0', '\0', '\0'};
    char *prefix = NULL;
    int decode_md = 0;
    int C;
    int start, end;
    char ref_name[1024] = {0}, *arg_list, *ref_fn = NULL;
    int embed_ref = 0;

    while ((C = getopt(argc, argv, "bu0123456789mp:hr:R:X")) != -1) {
	switch (C) {
	case 'b':
	    mode[1] = 'b';
	    break;

	case 'u':
	    mode[2] = '0';
	    break;

	case '0': case '1': case '2': case '3': case '4':
	case '5': case '6': case '7': case '8': case '9':
	    mode[2] = C;
	    break;

	case 'm':
	    decode_md = 1;
	    break;

	case 'p':
	    prefix = optarg;
	    break;

	case 'h':
	    usage(stdout);
	    return 0;

	case 'r':
	    ref_fn = optarg;
	    break;

	case 'X':
	    embed_ref = 1;
	    break;

	case 'R': {
	    char *cp = strchr(optarg, ':');
	    if (cp) {
		*cp = 0;
		switch (sscanf(cp+1, "%d-%d", &start, &end)) {
		case 1:
		    end = start;
		    break;
		case 2:
		    break;
		default:
		    fprintf(stderr, "Malformed range format\n");
		    return 1;
		}
	    } else {
		start = INT_MIN;
		end   = INT_MAX;
	    }
	    strncpy(ref_name, optarg, 1023);
	    break;
	}

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

    if (argc - optind == 1) {
	if (NULL == (bfd = bam_open("-", mode))) {
	    fprintf(stderr, "Failed to open SAM/BAM output\n.");
	    return 1;
	}
    } else {
	if (NULL == (bfd = bam_open(argv[optind+1], mode))) {
	    fprintf(stderr, "Failed to open SAM/BAM output\n.");
	    perror(argv[optind+1]);
	    return 1;
	}
    }

    if (NULL == (fd = cram_open(argv[optind], "rb"))) {
	fprintf(stderr, "Error opening CRAM file '%s'.\n", argv[optind]);
	return 1;
    }

    if (*ref_name != 0)
	cram_index_load(fd, argv[optind]);

    if (prefix)
	cram_set_option(fd, CRAM_OPT_PREFIX, prefix);

    if (decode_md)
	cram_set_option(fd, CRAM_OPT_DECODE_MD, decode_md);

    if (embed_ref)
	cram_set_option(fd, CRAM_OPT_EMBED_REF, embed_ref);

    /* Find and load reference */
    cram_load_reference(fd, ref_fn);
    if (!fd->refs && !embed_ref) {
	fprintf(stderr, "Unable to find an appropriate reference.\n"
		"Please specify a valid reference with -r ref.fa option.\n");
	return 1;
    }

    bfd->header = fd->header;

    if (*ref_name != 0) {
	cram_range r;
	int refid = sam_hdr_name2ref(fd->header, ref_name);

	if (refid == -1 && *ref_name != '*') {
	    fprintf(stderr, "Unknown reference name '%s'\n", ref_name);
	    return 1;
	}
	r.refid = refid;
	r.start = start;
	r.end = end;
	cram_set_option(fd, CRAM_OPT_RANGE, &r);
    }

    /* SAM Header */
    if (!(arg_list = stringify_argv(argc, argv)))
	return 1;
    sam_hdr_add_PG(bfd->header, "cram_to_sam",
		   "VN", PACKAGE_VERSION,
		   "CL", arg_list, NULL);
    free(arg_list);

    bam_write_header(bfd);

    while (cram_get_bam_seq(fd, &bam) == 0) {
	bam_put_seq(bfd, bam);
    }

    if (!cram_eof(fd)) {
	fprintf(stderr, "Error while reading file\n");
	return 1;
    }

    cram_close(fd);

    bfd->header = NULL;
    bam_close(bfd);

    free(bam);

    return 0;
}
