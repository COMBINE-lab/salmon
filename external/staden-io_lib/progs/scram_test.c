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
    fprintf(fp, "  -=- scram_flagstat -=-     version %s\n", PACKAGE_VERSION);
    fprintf(fp, "Author: James Bonfield, Wellcome Trust Sanger Institute. 2013\n\n");

    fprintf(fp, "Usage:    scram_flagstat [options] [input_file]\n");

    fprintf(fp, "Options:\n");
    fprintf(fp, "    -I format      Set input format:  \"bam\", \"sam\" or \"cram\".\n");
    fprintf(fp, "    -R range       [Cram] Specifies the refseq:start-end range\n");
    fprintf(fp, "    -r ref.fa      [Cram] Specifies the reference file.\n");
    fprintf(fp, "    -t N           Use N threads (availability varies by format)\n");
}

typedef struct {
    long long n_reads[2], n_mapped[2], n_pair_all[2], n_pair_map[2], n_pair_good[2];
    long long n_sgltn[2], n_read1[2], n_read2[2];
    long long n_dup[2];
    long long n_diffchr[2], n_diffhigh[2];
} bam_flagstat_t;

int main(int argc, char **argv) {
    scram_fd *in;
    bam_seq_t *s;
    char imode[10], *in_f = "";
    int level = '\0'; // nul terminate string => auto level
    int c;
    char *ref_fn = NULL;
    int start, end, ignore_md5 = 0;
    char ref_name[1024] = {0};
    bam_flagstat_t st;
    int nthreads = 1;
    int col = SAM_FLAG, req = 0;

    memset(&st, 0, sizeof(st));

    /* Parse command line arguments */
    while ((c = getopt(argc, argv, "hI:R:r:!t:c:C:")) != -1) {
	switch (c) {
	case 'h':
	    usage(stdout);
	    return 0;

	case 'r':
	    ref_fn = optarg;
	    break;

	case 'I':
	    in_f = parse_format(optarg);
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

	case '!':
	    ignore_md5 = 1;
	    break;

	case 't':
	    nthreads = atoi(optarg);
	    if (nthreads < 1) {
		fprintf(stderr, "Number of threads needs to be >= 1\n");
		return 1;
	    }
	    break;

	case 'c':
	    req = col = strtol(optarg, NULL, 0);
	    break;

	case 'C':
	    req = strtol(optarg, NULL, 0);
	    break;

	case '?':
	    fprintf(stderr, "Unrecognised option: -%c\n", optopt);
	    usage(stderr);
	    return 1;
	}
    }    

    if (argc - optind > 2) {
	fprintf(stderr, "Usage: scramble [input_file [output_file]]\n");
	return 1;
    }
    

    /* Open up input and output files */
    sprintf(imode, "r%s%c", in_f, level);
    if (argc - optind > 0) {
	if (*in_f == 0)
	    sprintf(imode, "r%s%c", detect_format(argv[optind]), level);
	if (!(in = scram_open(argv[optind], imode))) {
	    fprintf(stderr, "Failed to open file %s\n", argv[optind]);
	    return 1;
	}
    } else {
	if (!(in = scram_open("-", imode))) {
	    fprintf(stderr, "Failed to open file %s\n", argv[optind]);
	    return 1;
	}
    }
    if (!in->is_bam && ref_fn)
	cram_load_reference(in->c, ref_fn);

    if (nthreads > 1) 
	if (scram_set_option(in,  CRAM_OPT_NTHREADS, nthreads))
	    return 1;


    if (ignore_md5)
	if (scram_set_option(in, CRAM_OPT_IGNORE_MD5, ignore_md5))
	    return 1;

    scram_set_option(in, CRAM_OPT_REQUIRED_FIELDS, req);

    /* Support for sub-range queries, currently implemented for CRAM only */
    if (*ref_name != 0) {
	cram_range r;
	int refid;

	if (in->is_bam) {
	    fprintf(stderr, "Currently the -R option is only implemented for CRAM indices\n");
	    return 1;
	}
	    
	cram_index_load(in->c, argv[optind]);

	refid = sam_hdr_name2ref(in->c->header, ref_name);

	if (refid == -1 && *ref_name != '*') {
	    fprintf(stderr, "Unknown reference name '%s'\n", ref_name);
	    return 1;
	}
	r.refid = refid;
	r.start = start;
	r.end = end;
	if (scram_set_option(in, CRAM_OPT_RANGE, &r))
	    return 1;
    }

    /* Do the actual file format conversion */
    s = NULL;
    while (scram_get_seq(in, &s) >= 0) {
	if (col & SAM_QNAME)
	    printf("Name:\t%.*s\n", bam_name_len(s), bam_name(s));

	if (col & SAM_FLAG)
	    printf("Flag:\t%d\n", bam_flag(s));

	if (col & SAM_RNAME)
	    printf("Ref:\t%d\n", bam_ref(s));

	if (col & SAM_POS)
	    printf("Pos:\t%d\n", bam_pos(s));

	if (col & SAM_MAPQ)
	    printf("MapQ:\t%d\n", bam_map_qual(s));

	if (col & SAM_CIGAR) {
	    int i;
	    uint32_t *cig = bam_cigar(s);
	    printf("Cigar:\t");
	    for (i = 0; i < bam_cigar_len(s); i++)
		printf("%d%c", cig[i]>>4, "MIDNSHP=X"[cig[i]&0xf]);
	    putchar('\n');
	}

	if (col & SAM_RNEXT)
	    printf("Rnext:\t%d\n", bam_mate_ref(s));

	if (col & SAM_PNEXT)
	    printf("Pnext:\t%d\n", bam_mate_pos(s));

	if (col & SAM_TLEN)
	    printf("Tlen:\t%d\n", bam_ins_size(s));

	if (col & SAM_SEQ) {
	    unsigned char *seq = (unsigned char *) bam_seq(s);
	    int i;
	    printf("Seq:\t");
	    for (i = 0; i < bam_seq_len(s); i++)
		putchar("=ACMGRSVTWYHKDBN"[i&1 ? seq[i/2]&15 : seq[i/2]>>4]);
	    putchar('\n');
	}

	if (col & SAM_QUAL) {
	    char *qual = bam_qual(s);
	    int i;
	    printf("Qual:\t");
	    for (i = 0; i < bam_seq_len(s); i++)
		putchar(qual[i]+'!');
	    putchar('\n');
	}

	if (col & SAM_AUX) {
	    char *aux = bam_aux(s);
	    int i, aux_size = ((char *)s + bam_blk_size(s)) - aux;
	    printf("aux:\t");
	    for (i = 0; i < aux_size; i++)
		printf("%02x", aux[i]);
	    putchar('\n');
	}
    }

    if (s)
	free(s);

    if (!scram_eof(in))
	return 1;

    if (scram_close(in))
	return 1;

    return 0;
}
