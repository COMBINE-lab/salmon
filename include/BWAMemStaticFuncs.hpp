#ifndef BWAMEM_STATIC_FUNCS_HPP
#define BWAMEM_STATIC_FUNCS_HPP

extern unsigned char nst_nt4_table[256];
char const* bwa_pg = "cha";

/******* STUFF THAT IS STATIC IN BWAMEM THAT WE NEED HERE --- Just re-define it
 * *************/
#define intv_lt(a, b) ((a).info < (b).info)
KSORT_INIT(mem_intv, bwtintv_t, intv_lt)

typedef struct {
  bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;

static smem_aux_t* smem_aux_init() {
  smem_aux_t* a;
  a = static_cast<smem_aux_t*>(calloc(1, sizeof(smem_aux_t)));
  a->tmpv[0] = static_cast<bwtintv_v*>(calloc(1, sizeof(bwtintv_v)));
  a->tmpv[1] = static_cast<bwtintv_v*>(calloc(1, sizeof(bwtintv_v)));
  return a;
}

static void smem_aux_destroy(smem_aux_t* a) {
  free(a->tmpv[0]->a);
  free(a->tmpv[0]);
  free(a->tmpv[1]->a);
  free(a->tmpv[1]);
  free(a->mem.a);
  free(a->mem1.a);
  free(a);
}

static void mem_collect_intv(const SalmonOpts& sopt, const mem_opt_t* opt,
                             SalmonIndex* sidx, int len, const uint8_t* seq,
                             smem_aux_t* a) {
  const bwt_t* bwt = sidx->bwaIndex()->bwt;
  size_t i, k, x = 0, old_n;
  int start_width = (opt->flag & MEM_F_SELF_OVLP) ? 2 : 1;
  int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);
  a->mem.n = 0;

  // first pass: find all SMEMs
  if (sidx->hasAuxKmerIndex()) {
    KmerIntervalMap& auxIdx = sidx->auxIndex();
    int32_t klen = static_cast<int32_t>(auxIdx.k());
    while (x < static_cast<size_t>(len)) {
      if (seq[x] < 4) {
        // Make sure there are at least k bases left
        if (len - static_cast<int>(x) < klen) {
          x = len;
          continue;
        }
        // search for this key in the auxiliary index
        KmerKey kmer(const_cast<uint8_t*>(&(seq[x])), klen);
        auto it = auxIdx.find(kmer);
        // if we can't find it, move to the next key
        if (it == auxIdx.end()) {
          ++x;
          continue;
        }
        // otherwise, start the search using the initial interval @it->second
        // from the hash
        int xb = x;
        x = bwautils::bwt_smem1_with_kmer(bwt, len, seq, x, start_width,
                                          it->second, &a->mem1, a->tmpv);
        for (i = 0; i < a->mem1.n; ++i) {
          bwtintv_t* p = &a->mem1.a[i];
          int slen = (uint32_t)p->info - (p->info >> 32); // seed length
          if (slen >= opt->min_seed_len)
            kv_push(bwtintv_t, a->mem, *p);
        }
      } else
        ++x;
    }
  } else {
    while (x < static_cast<size_t>(len)) {
      if (seq[x] < 4) {
        x = bwt_smem1(bwt, len, seq, x, start_width, &a->mem1, a->tmpv);
        for (i = 0; i < a->mem1.n; ++i) {
          bwtintv_t* p = &a->mem1.a[i];
          int slen = (uint32_t)p->info - (p->info >> 32); // seed length
          if (slen >= opt->min_seed_len)
            kv_push(bwtintv_t, a->mem, *p);
        }
      } else
        ++x;
    }
  }

  // For sensitive / extra-sensitive mode only
  if (sopt.sensitive or sopt.extraSeedPass) {
    // second pass: find MEMs inside a long SMEM
    old_n = a->mem.n;
    for (k = 0; k < old_n; ++k) {
      bwtintv_t* p = &a->mem.a[k];
      int start = p->info >> 32, end = (int32_t)p->info;
      if (end - start < split_len || p->x[2] > static_cast<bwtint_t>(opt->split_width))
        continue;

      // int idx = (start + end) >> 1;
      bwt_smem1(bwt, len, seq, (start + end) >> 1, p->x[2] + 1, &a->mem1,
                a->tmpv);
      for (i = 0; i < a->mem1.n; ++i)
        if ((uint32_t)a->mem1.a[i].info - (a->mem1.a[i].info >> 32) >=
            static_cast<uint32_t>(opt->min_seed_len))
          kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
    }
  }

  // For extra-sensitive mode only
  // third pass: LAST-like
  if (sopt.extraSeedPass and opt->max_mem_intv > 0) {
    x = 0;
    while (x < static_cast<size_t>(len)) {
      if (seq[x] < 4) {
        if (1) {
          bwtintv_t m;
          x = bwt_seed_strategy1(bwt, len, seq, x, opt->min_seed_len,
                                 opt->max_mem_intv, &m);
          if (m.x[2] > 0)
            kv_push(bwtintv_t, a->mem, m);
        } else { // for now, we never come to this block which is slower
          x = bwt_smem1a(bwt, len, seq, x, start_width, opt->max_mem_intv,
                         &a->mem1, a->tmpv);
          for (i = 0; i < a->mem1.n; ++i)
            kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
        }
      } else
        ++x;
    }
  }
  // sort
  // ks_introsort(mem_intv, a->mem.n, a->mem.a);
}

/******* END OF STUFF THAT IS STATIC IN BWAMEM THAT WE NEED HERE --- Just
 * re-define it *************/

#endif // BWAMEM_STATIC_FUNCS_HPP
