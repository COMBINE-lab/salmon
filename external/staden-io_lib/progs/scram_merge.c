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
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013
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

#if defined(__MINGW32__) || defined(__FreeBSD__) || defined(__APPLE__)
#   include <getopt.h>
#endif

#include <io_lib/scram.h>
#include <io_lib/os.h>

/*
 * Return 1 for compatible
 *        0 for incompatible
 */
static int hdr_compare(SAM_hdr *h1, SAM_hdr *h2) {
    int i;
    if (h1->nref != h2->nref)
	return 0;

    for (i = 0; i < h1->nref; i++) {
	if (strcmp(h1->ref[i].name, h2->ref[i].name) != 0)
	    return 0;
	if (h1->ref[i].len != h2->ref[i].len)
	    return 0;
    }

    return 1;
}

static char *parse_format(char *str) {
    if (strcmp(str, "sam") == 0 || strcmp(str, "SAM") == 0)
	return "";

    if (strcmp(str, "bam") == 0 || strcmp(str, "BAM") == 0)
	return "b";

    if (strcmp(str, "cram") == 0 || strcmp(str, "CRAM") == 0)
	return "c";

    fprintf(stderr, "Unrecognised file format '%s'\n", str);
    exit(1);
}

static char *detect_format(char *fn) {
    char *cp = strrchr(fn, '.');
    
    if (!cp)
	return "";

    if (strcmp(cp, ".sam") == 0 || strcmp(cp, ".SAM") == 0)
	return "";
    if (strcmp(cp, ".bam") == 0 || strcmp(cp, ".BAM") == 0)
	return "b";
    if (strcmp(cp, ".cram") == 0 || strcmp(cp, ".CRAM") == 0)
	return "c";

    return "";
}

static void usage(FILE *fp) {
    fprintf(fp, "  -=- scram_merge -=-     version %s\n", PACKAGE_VERSION);
    fprintf(fp, "Author: James Bonfield, Wellcome Trust Sanger Institute. 2013\n\n");

    fprintf(fp, "Usage:    scram_merge [options] input_file ...\n");

    fprintf(fp, "Options:\n");
    fprintf(fp, "    -I format      Set input format:  \"bam\", \"sam\" or \"cram\".\n");
    fprintf(fp, "    -O format      Set output format: \"bam\", \"sam\" or \"cram\".\n");
    fprintf(fp, "    -1 to -9       Set zlib compression level.\n");
    fprintf(fp, "    -0 or -u       No zlib compression.\n");
    //fprintf(fp, "    -v             Verbose output.\n");
    fprintf(fp, "    -R range       [Cram] Specifies the refseq:start-end range\n");
    fprintf(fp, "    -r ref.fa      [Cram] Specifies the reference file.\n");
    fprintf(fp, "    -s integer     [Cram] Sequences per slice, default %d.\n",
	    SEQS_PER_SLICE);
    fprintf(fp, "    -S integer     [Cram] Slices per container, default %d.\n",
	    SLICE_PER_CNT);
    fprintf(fp, "    -V version     [Cram] Specify the file format version to write (eg 1.1, 2.0)\n");
    fprintf(fp, "    -X             [Cram] Embed reference sequence.\n");
}

int main(int argc, char **argv) {
    scram_fd **in, *out;
    int n_input, i;
    bam_seq_t **s;
    char imode[10], *in_f = "", omode[10], *out_f = "";
    int level = '\0'; // nul terminate string => auto level
    int c, verbose = 0;
    int s_opt = 0, S_opt = 0, embed_ref = 0;
    char *ref_fn = NULL;
    int start, end;
    char ref_name[1024] = {0};
    refs_t *refs = NULL;
    int max_reads = -1;

    /* Parse command line arguments */
    while ((c = getopt(argc, argv, "u0123456789hvs:S:V:r:XI:O:R:N:")) != -1) {
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

	case 'I':
	    in_f = parse_format(optarg);
	    break;

	case 'O':
	    out_f = parse_format(optarg);
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

	case 'N': // For debugging
	    max_reads = atoi(optarg);
	    break;

	case '?':
	    fprintf(stderr, "Unrecognised option: -%c\n", optopt);
	    usage(stderr);
	    return 1;
	}
    }    

    /* Open output file */
    sprintf(omode, "w%s%c", out_f, level);
    if (!(out = scram_open("-", omode))) {
	fprintf(stderr, "Failed to open bam file %s\n", argv[optind+1]);
	return 1;
    }

    /* Open multiple input files */
    sprintf(imode, "r%s%c", in_f, level);
    n_input = argc - optind;
    if (!n_input) {
	fprintf(stderr, "No input files specified.\n");
	return 1;
    }
    if (!(in = malloc(n_input * sizeof(*in))))
	return 1;
    if (!(s = malloc(n_input * sizeof(*s))))
	return 1;
    for (i = 0; i < n_input; i++, optind++) {
	s[i] = NULL;
	if (*in_f == 0)
	    sprintf(imode, "r%s%c", detect_format(argv[optind]), level);
	if (!(in[i] = scram_open(argv[optind], imode))) {
	    fprintf(stderr, "Failed to open bam file %s\n", argv[optind]);
	    return 1;
	}
	if (i && !hdr_compare(scram_get_header(in[0]),
			      scram_get_header(in[i]))) {
	    fprintf(stderr, "Incompatible reference sequence list.\n");
	    fprintf(stderr, "Currently the @SQ lines need to be identical"
		    " in all files.\n");
	    return 1;
	}

	if (!refs && scram_get_refs(in[i]))
	    refs = scram_get_refs(in[i]);

	if (refs && scram_set_option(in[i], CRAM_OPT_SHARED_REF, refs))
	    return 1;

	/* Support for sub-range queries, currently implemented for CRAM only */
	if (*ref_name != 0) {
	    cram_range r;
	    int refid;

	    if (in[i]->is_bam) {
		fprintf(stderr, "Currently the -R option is only implemented for CRAM indices\n");
		return 1;
	    }
	    
	    cram_index_load(in[i]->c, argv[optind]);

	    refid = sam_hdr_name2ref(in[i]->c->header, ref_name);


	    if (refid == -1 && *ref_name != '*') {
		fprintf(stderr, "Unknown reference name '%s'\n", ref_name);
		return 1;
	    }
	    r.refid = refid;
	    r.start = start;
	    r.end = end;

	    if (scram_set_option(in[i], CRAM_OPT_RANGE, &r))
	    	return 1;
	}
    }

    /* Set any format specific options */
    if (refs)
	scram_set_option(out, CRAM_OPT_SHARED_REF, refs);

    if (scram_set_option(out, CRAM_OPT_VERBOSITY, verbose))
	return 1;
    if (s_opt)
	if (scram_set_option(out, CRAM_OPT_SEQS_PER_SLICE, s_opt))
	    return 1;

    if (S_opt)
	if (scram_set_option(out, CRAM_OPT_SLICES_PER_CONTAINER, S_opt))
	    return 1;

    if (embed_ref)
	if (scram_set_option(out, CRAM_OPT_EMBED_REF, embed_ref))
	    return 1;
    
    /* Copy header and refs from in to out, for writing purposes */
    // FIXME: do proper merging of @PG lines
    // FIXME: track mapping of old PG aux name to new PG aux name per seq
    scram_set_header(out, sam_hdr_dup(scram_get_header(in[0])));

    // Needs doing after loading the header.
    if (ref_fn)
	if (scram_set_option(out, CRAM_OPT_REFERENCE, ref_fn))
	    return 1;

    if (scram_get_header(in[0])) {
	if (scram_write_header(out))
	    return 1;
    }


    /* Do the actual file format conversion */
    fprintf(stderr, "Opening and loading initial seqs\n");
    for (i = 0; i < n_input; i++) {
	if (scram_get_seq(in[i], &s[i]) < 0) {
	    if (scram_close(in[i]))
		return 1;
	    in[i] = NULL;
	    free(s[i]);
	    continue;
	}
    }

    fprintf(stderr, "Merging...\n");
    for (;;) {
	int64_t best_val = INT64_MAX;
	int best_j = 0, j;

	for (j = 0; j < n_input; j++) {
	    bam_seq_t *b = s[j];
	    uint64_t x;
	    if (!in[j])
		continue;

	    x = (((uint64_t)bam_ref(b))<<33)
		| (bam_pos(b)<<2)
		| (bam_strand(b)<<1)
		| !(bam_flag(b) & BAM_FREAD1);
	    if (best_val > x) {
		best_val = x;
		best_j = j;
	    }
	}
	
	if (best_val == INT64_MAX) { // all closed
	    break;
	}

	if (-1 == scram_put_seq(out, s[best_j]))
	    return 1;
	
	if (scram_get_seq(in[best_j], &s[best_j]) < 0) {
	    if (scram_close(in[best_j]))
		return 1;
	    in[best_j] = NULL;
	    free(s[best_j]);
	}

	if (max_reads >= 0)
	    if (--max_reads == 0)
		break;
    }

    for (i = 0; i < n_input; i++) {
	if (!in[i])
	    continue;
	scram_close(in[i]);
	if (s[i])
	    free(s[i]);
    }

    /* Finally tidy up and close files */
    if (scram_close(out))
	return 1;
    free(in);
    free(s);

    return 0;
}
