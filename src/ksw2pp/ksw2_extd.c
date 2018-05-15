#include <stdio.h> // for debugging only
#include "ksw2pp/ksw2.h"

typedef struct { int32_t h, e, e2; } eh_t;

void ksw_extd(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
			  int8_t gapo, int8_t gape, int8_t gapo2, int8_t gape2, int w, int zdrop, int flag, ksw_extz_t *ez)
{
	eh_t *eh;
	int8_t *qp; // query profile
	int32_t i, j, k, max_j = 0, gapoe = gapo + gape, gapoe2 = gapo2 + gape2, n_col, *off = 0, with_cigar = !(flag&KSW_EZ_SCORE_ONLY);
	uint8_t *z = 0; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be more complex

	ksw_reset_extz(ez);

	// allocate memory
	if (w < 0) w = tlen > qlen? tlen : qlen;
	n_col = qlen < 2*w+1? qlen : 2*w+1; // maximum #columns of the backtrack matrix
	qp = (int8_t*)kmalloc(km, qlen * m);
	eh = (eh_t*)kcalloc(km, qlen + 1, sizeof(eh_t));
	if (with_cigar) {
		z = (uint8_t*)kmalloc(km, (size_t)n_col * tlen);
		off = (int32_t*)kcalloc(km, tlen, 4);
	}

	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}

	// fill the first row
	eh[0].h = 0, eh[0].e = -gapoe - gapoe, eh[0].e2 = -gapoe2 - gapoe2;
	for (j = 1; j <= qlen && j <= w; ++j) {
		int tmp;
		eh[j].h = -(gapo + gape * j) > -(gapo2 + gape2 * j)? -(gapo + gape * j) : -(gapo2 + gape2 * j);
		tmp = -(gapoe + gape * j) > -(gapoe2 + gape2 * j)? -(gapoe + gape * j) : -(gapoe2 + gape2 * j);
		eh[j].e = tmp - gapoe;
		eh[j].e2 = tmp - gapoe2;
	}
	for (; j <= qlen; ++j) eh[j].h = eh[j].e = eh[j].e2 = KSW_NEG_INF; // everything is -inf outside the band

	// DP loop
	for (i = 0; i < tlen; ++i) { // target sequence is in the outer loop
		int32_t f, f2, h1, st, en, max = KSW_NEG_INF, tmp;
		int8_t *q = &qp[target[i] * qlen];
		st = i > w? i - w : 0;
		en = i + w < qlen - 1? i + w : qlen - 1;
		tmp = -(gapoe + gape * i) > -(gapoe2 + gape2 * i)? -(gapoe + gape * i) : -(gapoe2 + gape2 * i);
		h1 = st > 0? KSW_NEG_INF : tmp;
		f  = st > 0? KSW_NEG_INF : tmp - gapoe;
		f2 = st > 0? KSW_NEG_INF : tmp - gapoe2;
		if (!with_cigar) {
			for (j = st; j <= en; ++j) {
				eh_t *p = &eh[j];
				int32_t h = p->h, h2, e = p->e, e2 = p->e2;
				p->h = h1;
				h += q[j];
				h = h >= e?  h : e;
				h = h >= f?  h : f;
				h = h >= e2? h : e2;
				h = h >= f2? h : f2;
				h1 = h;
				max_j = max > h? max_j : j;
				max   = max > h? max   : h;
				h -= gapoe;
				e -= gape;
				e  = e > h? e : h;
				p->e = e;
				f -= gape;
				f  = f > h? f : h;
				h2 = h1 - gapoe2;
				e2-= gape2;
				e2 = e2 > h2? e2 : h2;
				p->e2 = e2;
				f2-= gape2;
				f2 = f2 > h2? f2 : h2;
			}
		} else if (!(flag&KSW_EZ_RIGHT)) {
			uint8_t *zi = &z[(long)i * n_col];
			off[i] = st;
			for (j = st; j <= en; ++j) {
				eh_t *p = &eh[j];
				int32_t h = p->h, h2, e = p->e, e2 = p->e2;
				uint8_t d; // direction
				p->h = h1;
				h += q[j];
				d = h >= e?  0 : 1;
				h = h >= e?  h : e;
				d = h >= f?  d : 2;
				h = h >= f?  h : f;
				d = h >= e2? d : 3;
				h = h >= e2? h : e2;
				d = h >= f2? d : 4;
				h = h >= f2? h : f2;
				h1 = h;
				max_j = max > h? max_j : j;
				max   = max > h? max   : h;
				h -= gapoe;
				e -= gape;
				d |= e > h? 1<<3 : 0;
				e  = e > h? e    : h;
				p->e = e;
				f -= gape;
				d |= f > h? 1<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
				f  = f > h? f    : h;
				h2 = h1 - gapoe2;
				e2-= gape2;
				d |= e2 > h2? 1<<5 : 0;
				e2 = e2 > h2? e2 : h2;
				p->e2 = e2;
				f2-= gape2;
				d |= f2 > h2? 1<<6 : 0;
				f2 = f2 > h2? f2 : h2;
				zi[j - st] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
			}
		} else {
			uint8_t *zi = &z[(long)i * n_col];
			off[i] = st;
			for (j = st; j <= en; ++j) {
				eh_t *p = &eh[j];
				int32_t h = p->h, h2, e = p->e, e2 = p->e2;
				uint8_t d; // direction
				p->h = h1;
				h += q[j];
				d = h > e?  0 : 1;
				h = h > e?  h : e;
				d = h > f?  d : 2;
				h = h > f?  h : f;
				d = h > e2? d : 3;
				h = h > e2? h : e2;
				d = h > f2? d : 4;
				h = h > f2? h : f2;
				h1 = h;
				max_j = max > h? max_j : j;
				max   = max > h? max   : h;
				h -= gapoe;
				e -= gape;
				d |= e >= h? 1<<3 : 0;
				e  = e >= h? e    : h;
				p->e = e;
				f -= gape;
				d |= f >= h? 1<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
				f  = f >= h? f    : h;
				h2 = h1 - gapoe2;
				e2-= gape2;
				d |= e2 >= h2? 1<<5 : 0;
				e2 = e2 >= h2? e2 : h2;
				p->e2 = e2;
				f2-= gape2;
				d |= f2 >= h2? 1<<6 : 0;
				f2 = f2 >= h2? f2 : h2;
				zi[j - st] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
			}
		}
		eh[j].h = h1, eh[j].e = KSW_NEG_INF;
		// update ez
		if (en == qlen - 1 && eh[qlen].h > ez->mqe)
			ez->mqe = eh[qlen].h, ez->mqe_t = i;
		if (i == tlen - 1)
			ez->mte = max, ez->mte_q = max_j;
		if (ksw_apply_zdrop(ez, 0, max, i, max_j, zdrop, gape2)) break;
		if (i == tlen - 1 && en == qlen - 1)
			ez->score = eh[qlen].h;
	}
	kfree(km, qp); kfree(km, eh);
	if (with_cigar) {
		int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
		if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY))
			ksw_backtrack(km, 0, rev_cigar, 0, z, off, 0, n_col, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		else if (ez->max_t >= 0 && ez->max_q >= 0)
			ksw_backtrack(km, 0, rev_cigar, 0, z, off, 0, n_col, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		kfree(km, z); kfree(km, off);
	}
}
