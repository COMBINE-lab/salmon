#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <string>
#include "ksw2.h"
#include "kseq.h"
#include "KSW2Aligner.hpp"
KSEQ_INIT(gzFile, gzread)

#ifdef HAVE_GABA
#include "gaba.h"
#endif

#ifdef HAVE_PARASAIL
#include "parasail.h"
#endif

using ksw2pp::KSW2Aligner;
using ksw2pp::KSW2Config;
using ksw2pp::EnumToType;
using ksw2pp::KSW2AlignmentType;

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
}

static void global_aln(const char *algo, void *km, const char *qseq_, const char *tseq_, int8_t m, const int8_t *mat, int8_t q, int8_t e, int8_t q2, int8_t e2,
                       int w, int zdrop, int flag, ksw_extz_t *ez, KSW2Aligner& aligner)
{
	int i, qlen, tlen;
	uint8_t *qseq, *tseq;
	ez->max = 0, ez->mqe = ez->mte = KSW_NEG_INF;
	ez->n_cigar = 0;
	qlen = strlen(qseq_);
	tlen = strlen(tseq_);
	qseq = (uint8_t*)calloc(qlen + 33, 1); // 32 for gaba
	tseq = (uint8_t*)calloc(tlen + 33, 1);
	for (i = 0; i < qlen; ++i)
		qseq[i] = seq_nt4_table[(uint8_t)qseq_[i]];
	for (i = 0; i < tlen; ++i)
		tseq[i] = seq_nt4_table[(uint8_t)tseq_[i]];
	if (strcmp(algo, "gg") == 0) {
		if (flag & KSW_EZ_SCORE_ONLY) ez->score = ksw_gg(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, 0, 0, 0);
		else ez->score = ksw_gg(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
	} else if (strcmp(algo, "gg2") == 0) {
		if (flag & KSW_EZ_SCORE_ONLY) ez->score = ksw_gg2(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, 0, 0, 0);
		else ez->score = ksw_gg2(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
	} else if (strcmp(algo, "gg2_sse_aligner") == 0) {
    ez->score = aligner(qseq, qlen, tseq, tlen, ez, EnumToType<KSW2AlignmentType::GLOBAL>());
  } else if (strcmp(algo, "gg2_sse") == 0)     {
    ez->score = ksw_gg2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
  }
	else if (strcmp(algo, "extz") == 0)        ksw_extz(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, zdrop, flag, ez);
	else if (strcmp(algo, "extz2_sse") == 0)   ksw_extz2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, zdrop, flag, ez);
  else if (strcmp(algo, "extz2_sse_aligner") == 0) aligner(qseq, qlen, tseq, tlen, ez, EnumToType<KSW2AlignmentType::EXTENSION>());
	else if (strcmp(algo, "extd") == 0)        ksw_extd(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, q2, e2, w, zdrop, flag, ez);
	else if (strcmp(algo, "extd2_sse") == 0)   ksw_extd2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, q2, e2, w, zdrop, flag, ez);
	else if (strcmp(algo, "extf2_sse") == 0)   ksw_extf2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, mat[0], mat[1], e, w, zdrop, ez);
	else if (strcmp(algo, "exts2_sse") == 0) {
		int8_t mat[25];
		ksw_gen_simple_mat(5, mat, 1, 2);
		ksw_exts2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, 2, 1, 32, 4, zdrop, flag|KSW_EZ_SPLICE_FOR, ez);
	}
	else if (strcmp(algo, "test") == 0) ksw_extd2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, 4, 2, 24, 1, 751, 400, 8, ez);
#ifdef HAVE_GABA
	else if (strcmp(algo, "gaba") == 0) { // libgaba. Note that gaba may not align to the end
		int buf_len = 0x10000;
		char *cigar = (char*)calloc(buf_len, 1);
		for (i = 0; i < qlen; ++i)
			qseq[i] = qseq[i] < 4? 1<<qseq[i] : 15;
		for (i = 0; i < tlen; ++i)
			tseq[i] = tseq[i] < 4? 1<<tseq[i] : 15;
		struct gaba_section_s qs = gaba_build_section(0, qseq, qlen);
		struct gaba_section_s ts = gaba_build_section(2, tseq, tlen);

		gaba_t *ctx = gaba_init(GABA_PARAMS(.xdrop = (zdrop<120?zdrop:120), GABA_SCORE_SIMPLE(mat[0], abs(mat[1]), q, e)));
		void const *lim = (void const *)0x800000000000;
		gaba_dp_t *dp = gaba_dp_init(ctx, lim, lim);
		struct gaba_fill_s *f = gaba_dp_fill_root(dp, &qs, 0, &ts, 0);
		f = gaba_dp_fill(dp, f, &qs, &ts);
		gaba_alignment_t *r = gaba_dp_trace(dp, 0, f, 0);
		ez->score = r->score;
		gaba_dp_dump_cigar_forward(cigar, buf_len, r->path->array, 0, r->path->len);
		printf("gaba cigar: %s\n", cigar);
		gaba_dp_clean(dp);
		gaba_clean(ctx);
		free(cigar);
	}
#endif
#ifdef HAVE_PARASAIL
	else if (strncmp(algo, "ps_", 3) == 0) {
		parasail_matrix_t *ps_mat;
		parasail_function_t *func;
		parasail_result_t *res;
		ps_mat = parasail_matrix_create("ACGTN", mat[0], mat[1]); // TODO: call parasail_matrix_set_value() to change N<->A/C/G/T scores
		func = parasail_lookup_function(algo + 3);
		if (func == 0) {
			fprintf(stderr, "ERROR: can't find parasail function '%s'\n", algo + 3);
			exit(1);
		}
		res = func(qseq_, qlen, tseq_, tlen, q+e, e, ps_mat);
		ez->score = res->score;
		parasail_matrix_free(ps_mat);
		parasail_result_free(res);
	}
#endif
	else {
		fprintf(stderr, "ERROR: can't find algorithm '%s'\n", algo);
		exit(1);
	}
	free(qseq); free(tseq);
}

static void print_aln(const char *tname, const char *qname, ksw_extz_t *ez)
{
	printf("%s\t%s\t%d", tname, qname, ez->score);
	printf("\t%d\t%d\t%d", ez->max, ez->max_t, ez->max_q);
	if (ez->n_cigar > 0) {
		int i;
		putchar('\t');
		for (i = 0; i < ez->n_cigar; ++i)
			printf("%d%c", ez->cigar[i]>>4, "MID"[ez->cigar[i]&0xf]);
	}
	putchar('\n');
}

int main(int argc, char *argv[])
{
	void *km = 0;
	int8_t a = 2, b = 4, q = 4, e = 2, q2 = 13, e2 = 1;
	int c, i, pair = 1, w = -1, flag = 0, rep = 1, zdrop = -1, no_kalloc = 0;
	char *algo = nullptr, *s;
	int8_t mat[25];
	ksw_extz_t ez;
	gzFile fp[2];
  KSW2Config config;
  config.gapo = q;
  config.gape = e;


	while ((c = getopt(argc, argv, "t:w:R:rsgz:A:B:O:E:K")) >= 0) {
		if (c == 't') algo = optarg;
		else if (c == 'w') {w = atoi(optarg); config.bandwidth = w;}
		else if (c == 'R') rep = atoi(optarg);
		else if (c == 'z') {zdrop = atoi(optarg); config.dropoff = zdrop;}
		else if (c == 'r') {flag |= KSW_EZ_RIGHT; config.flag |= KSW_EZ_RIGHT; }
    else if (c == 's') {flag |= KSW_EZ_SCORE_ONLY; config.flag |= KSW_EZ_SCORE_ONLY; }
    else if (c == 'g') {flag |= KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP; config.flag |= KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP; }
		else if (c == 'K') no_kalloc = 1;
		else if (c == 'A') a = atoi(optarg);
		else if (c == 'B') b = atoi(optarg);
		else if (c == 'O') {
			q = q2 = strtol(optarg, &s, 10);
			if (*s == ',') q2 = strtol(s+1, &s, 10);
		} else if (c == 'E') {
			e = e2 = strtol(optarg, &s, 10);
			if (*s == ',') e2 = strtol(s+1, &s, 10);
		}
	}
  if (algo == nullptr) {
    std::string algos = "extd";
    algo = new char[algos.length()+1];
    strcpy(algo, algos.c_str());
  }
	if (argc - optind < 2) {
		fprintf(stderr, "Usage: ksw2-test [options] <DNA-target> <DNA-query>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -t STR        algorithm: gg, gg2, gg2_sse, extz, extz2_sse, extd, extd2_sse [%s]\n", algo);
		fprintf(stderr, "  -R INT        repeat INT times (for benchmarking) [1]\n");
		fprintf(stderr, "  -w INT        band width [inf]\n");
		fprintf(stderr, "  -z INT        Z-drop [%d]\n", zdrop);
		fprintf(stderr, "  -r            gap right alignment\n");
		fprintf(stderr, "  -s            score only\n");
		fprintf(stderr, "  -A INT        match score [%d]\n", a);
		fprintf(stderr, "  -B INT        mismatch penalty [%d]\n", b);
		fprintf(stderr, "  -O INT[,INT]  gap open penalty [%d,%d]\n", q, q2);
		fprintf(stderr, "  -E INT[,INT]  gap extension penalty [%d,%d]\n", e, e2);
		return 1;
	}


  KSW2Aligner aligner(static_cast<int8_t>(a), static_cast<int8_t>(b));
  aligner.config() = config;

#ifdef HAVE_KALLOC
	km = no_kalloc? 0 : km_init();
#endif
	memset(&ez, 0, sizeof(ksw_extz_t));
	ksw_gen_simple_mat(5, mat, a, -b);
	fp[0] = gzopen(argv[optind], "r");
	fp[1] = gzopen(argv[optind+1], "r");


	if (fp[0] == 0 && fp[1] == 0) {
		global_aln(algo, km, argv[optind+1], argv[optind], 5, mat, q, e, q2, e2, w, zdrop, flag, &ez, aligner);
		print_aln("first", "second", &ez);
	} else if (fp[0] && fp[1]) {
		kseq_t *ks[2];
		ks[0] = kseq_init(fp[0]);
		ks[1] = kseq_init(fp[1]);
		if (pair) {
			while (kseq_read(ks[0]) > 0) {
				if (kseq_read(ks[1]) <= 0) break;
				for (i = 0; i < rep; ++i) 
					global_aln(algo, km, ks[1]->seq.s, ks[0]->seq.s, 5, mat, q, e, q2, e2, w, zdrop, flag, &ez, aligner);
				print_aln(ks[0]->name.s, ks[1]->name.s, &ez);
			}
		}
		kseq_destroy(ks[0]);
		kseq_destroy(ks[1]);
		gzclose(fp[0]);
		gzclose(fp[1]);
	}
	if (strcmp(algo, "gg2_sse_aligner") != 0 and strcmp(algo, "extz2_sse_aligner") != 0) {
    kfree(km, ez.cigar);
  } else {
    aligner.freeCIGAR(&ez);
  }
#ifdef HAVE_KALLOC
	km_destroy(km);
#endif
  delete [] algo;
	return 0;
}
