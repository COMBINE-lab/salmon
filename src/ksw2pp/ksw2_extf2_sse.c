#include <string.h>
#include "ksw2pp/ksw2.h"

#ifdef __SSE2__
#include <emmintrin.h>

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif

#ifdef KSW_CPU_DISPATCH
#ifdef __SSE4_1__
void ksw_extf2_sse41(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t mch, int8_t mis, int8_t e, int w, int xdrop, ksw_extz_t *ez)
#else
  void ksw_extf2_sse2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t mch, int8_t mis, int8_t e, int w, int xdrop, ksw_extz_t *ez)
#endif
#else
void ksw_extf2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t mch, int8_t mis, int8_t e, int w, int xdrop, ksw_extz_t *ez)
#endif // ~KSW_CPU_DISPATCH
{
	int32_t r, t, tlen_, qlen_, last_st, last_en, H0 = 0, last_H0_t = 0;
	uint8_t *qr, *sf, *mem;
	__m128i e2_, sc_mch_, sc_mis_, *u, *v, *s;

	ksw_reset_extz(ez);
	e2_ = _mm_set1_epi8(e * 2);
	sc_mch_ = _mm_set1_epi8(mch);
	sc_mis_ = _mm_set1_epi8(mis < 0? mis : -mis);
	tlen_ = (tlen + 15) / 16;
	qlen_ = (qlen + 15) / 16;
	if (w < 0) w = tlen > qlen? tlen : qlen;

	mem = (uint8_t*)kcalloc(km, tlen_ * 4 + qlen_ + 1, 16);
	u = (__m128i*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned
	v = u + tlen_, s = v + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 16;

	for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
	memcpy(sf, target, tlen);

	for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1, st0, en0, st_, en_;
		uint8_t v1, *qrr = qr + (qlen - 1 - r), *u8 = (uint8_t*)u, *v8 = (uint8_t*)v;
		__m128i v1_;
		// find the boundaries
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r-w+1)>>1) st = (r-w+1)>>1; // take the ceil
		if (en > (r+w)>>1) en = (r+w)>>1; // take the floor
		if (st > en) break;
		st0 = st, en0 = en;
		st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
		// set boundary conditions
		v1 = (st > 0 && st - 1 >= last_st && st - 1 <= last_en)? v8[st - 1] : 0;
		if (en >= r) u8[r] = 0;
		// core loop
		v1_ = _mm_cvtsi32_si128(v1);
		st_ = st / 16, en_ = en / 16;
		for (t = st0; t <= en0; t += 16) {
			__m128i sq, st, tmp;
			sq = _mm_loadu_si128((__m128i*)&sf[t]);
			st = _mm_loadu_si128((__m128i*)&qrr[t]);
			tmp = _mm_cmpeq_epi8(sq, st);
#ifdef __SSE4_1__
			tmp = _mm_blendv_epi8(sc_mis_, sc_mch_, tmp);
#else
			tmp = _mm_or_si128(_mm_andnot_si128(tmp, sc_mis_), _mm_and_si128(tmp, sc_mch_));
#endif
			_mm_storeu_si128((__m128i*)((uint8_t*)s + t), tmp);
		}
		for (t = st_; t <= en_; ++t) {
			__m128i z, vt1, ut, tmp;
			z = _mm_add_epi8(_mm_load_si128(&s[t]), e2_);
			vt1 = _mm_load_si128(&v[t]);                     // vt1 <- v[r-1][t..t+15]
			tmp = _mm_srli_si128(vt1, 15);                   // tmp <- v[r-1][t+15]
			vt1 = _mm_or_si128(_mm_slli_si128(vt1, 1), v1_); // vt1 <- v[r-1][t-1..t+14]
			v1_ = tmp;
			ut = _mm_load_si128(&u[t]);                      // ut <- u[t..t+15]
#ifdef __SSE4_1__
			z = _mm_max_epi8(z, vt1);                        // z = z > a? z : a (signed)
#else
			z = _mm_and_si128(z, _mm_cmpgt_epi8(z, _mm_setzero_si128()));  // z = z > 0? z : 0;
			z = _mm_max_epu8(z, vt1);                        // z = max(z, a); this works because both are non-negative
#endif
			z = _mm_max_epu8(z, ut);                         // z = max(z, b); this works because both are non-negative
			_mm_store_si128(&u[t], _mm_sub_epi8(z, vt1));    // u[r][t..t+15] <- z - v[r-1][t-1..t+14]
			_mm_store_si128(&v[t], _mm_sub_epi8(z, ut));     // v[r][t..t+15] <- z - u[r-1][t..t+15]
		}
		if (r > 0) {
			if (last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0) {
				int32_t d0 = v8[last_H0_t] - e, d1 = u8[last_H0_t + 1] - e;
				if (d0 > d1) H0 += d0;
				else H0 += d1, ++last_H0_t;
			} else if (last_H0_t >= st0 && last_H0_t <= en0) {
				H0 += v8[last_H0_t] - e;
			} else {
				++last_H0_t, H0 += u8[last_H0_t] - e;
			}
			if (H0 > ez->max) ez->max = H0, ez->max_t = last_H0_t, ez->max_q = r - last_H0_t;
			else if (xdrop >= 0 && ez->max - H0 > xdrop) break;
		} else H0 = v8[0] - e - e, last_H0_t = 0;
		last_st = st, last_en = en;
	}
	if (r == qlen + tlen - 1) ez->score = H0;
	else ez->zdropped = 1;
	kfree(km, mem);
}
#endif // __SSE2__
