#include "BWAUtils.hpp"

namespace bwautils {
static void bwt_reverse_intvs(bwtintv_v* p) {
  if (p->n > 1) {
    int j;
    for (j = 0; static_cast<bwtint_t>(j) < p->n >> 1; ++j) {
      bwtintv_t tmp = p->a[p->n - 1 - j];
      p->a[p->n - 1 - j] = p->a[j];
      p->a[j] = tmp;
    }
  }
}

// Function modified from bwt_smem1a:
// https://github.com/lh3/bwa/blob/eb428d7d31ced059ad39af2701a22ebe6d175657/bwt.c#L289
/**
 * Search for the k-mer of length @len starting at @q.
 * Return true if an interval is found for the k-mer and false
 * otherwise. The appropriate bwt interval will be placed
 * in @resInterval upon success.
 *
 */
bool getIntervalForKmer(const bwt_t* bwt, // the bwt index
                        int len,          // k-mer length
                        const uint8_t* q, // query
                        bwtintv_t& resInterval) {
  int i, j, c, ret;
  int x = 0;
  bwtintv_t ik, ok[4];
  //bwtintv_v a[2], *prev, *curr, *swap;
  bwtintv_v a[2], *curr, *swap;

  if (q[x] > 3)
    return false;
  // if (min_intv < 1) min_intv = 1; // the interval size should be at least 1
  kv_init(a[0]);
  kv_init(a[1]);
  //prev = &a[0]; // use the temporary vector if provided
  curr = &a[1];
  bwt_set_intv(bwt, q[x], ik); // the initial interval of a single base
  ik.info = x + 1;

  for (i = x + 1, curr->n = 0; i < len; ++i) { // forward search
    if (q[i] < 4) {                            // an A/C/G/T base
      c = 3 - q[i];                            // complement of q[i]
      bwt_extend(bwt, &ik, ok, 0);
      ik = ok[c];
      ik.info = i + 1;
    } else { // an ambiguous base
      break; // always terminate extension at an ambiguous base; in this case,
             // i<len always stands
    }
  }
  if (i == len) { // kv_push(bwtintv_t, *curr, ik); // push the last interval if
                  // we reach the end
    /// Copy over the final interval to our output
    resInterval.x[0] = ik.x[0];
    resInterval.x[1] = ik.x[1];
    resInterval.x[2] = ik.x[2];
    resInterval.info = ik.info;
    return true;
  } else { // we didn't find a k-mer of the requested length
    return false;
  }
}

// NOTE: $max_intv is not currently used in BWA-MEM
// NOTE: Modified from the original functions to take an initial interval for
// the search query
int bwt_smem1a_with_kmer(const bwt_t* bwt, int len, const uint8_t* q, int x,
                         int min_intv, uint64_t max_intv,
                         bwtintv_t initial_interval, bwtintv_v* mem,
                         bwtintv_v* tmpvec[2]) {
  int i, j, c, ret;
  bwtintv_t ik, ok[4];
  bwtintv_v a[2], *prev, *curr, *swap;

  mem->n = 0;
  if (q[x] > 3)
    return x + 1;
  if (min_intv < 1)
    min_intv = 1; // the interval size should be at least 1
  kv_init(a[0]);
  kv_init(a[1]);
  prev = tmpvec && tmpvec[0] ? tmpvec[0]
                             : &a[0]; // use the temporary vector if provided
  curr = tmpvec && tmpvec[1] ? tmpvec[1] : &a[1];
  // bwt_set_intv(bwt, q[x], ik); // the initial interval of a single base
  // ik.info = x + 1;

  // ROB: set ik to our initial interval
  ik.x[0] = initial_interval.x[0];
  ik.x[1] = initial_interval.x[1];
  ik.x[2] = initial_interval.x[2];
  // Is this last one right?
  int k = initial_interval.info;
  ik.info = x + k;

  for (i = x + k, curr->n = 0; i < len; ++i) { // forward search
    if (ik.x[2] < max_intv) {                  // an interval small enough
      kv_push(bwtintv_t, *curr, ik);
      break;
    } else if (q[i] < 4) { // an A/C/G/T base
      c = 3 - q[i];        // complement of q[i]
      bwt_extend(bwt, &ik, ok, 0);
      if (ok[c].x[2] != ik.x[2]) { // change of the interval size
        kv_push(bwtintv_t, *curr, ik);
        if (ok[c].x[2] < static_cast<bwtint_t>(min_intv))
          break; // the interval size is too small to be extended further
      }
      ik = ok[c];
      ik.info = i + 1;
    } else { // an ambiguous base
      kv_push(bwtintv_t, *curr, ik);
      break; // always terminate extension at an ambiguous base; in this case,
             // i<len always stands
    }
  }
  if (i == len)
    kv_push(bwtintv_t, *curr, ik); // push the last interval if we reach the end
  bwt_reverse_intvs(
      curr); // s.t. smaller intervals (i.e. longer matches) visited first
  ret = curr->a[0].info; // this will be the returned value
  swap = curr;
  curr = prev;
  prev = swap;

  for (i = x - 1; i >= -1; --i) { // backward search for MEMs
    c = i < 0
            ? -1
            : q[i] < 4 ? q[i] : -1; // c==-1 if i<0 or q[i] is an ambiguous base
    for (j = 0, curr->n = 0; static_cast<bwtint_t>(j) < prev->n; ++j) {
      bwtintv_t* p = &prev->a[j];
      if (c >= 0 && ik.x[2] >= max_intv)
        bwt_extend(bwt, p, ok, 1);
      if (c < 0 || ik.x[2] < max_intv ||
          ok[c].x[2] < static_cast<bwtint_t>(min_intv)) { // keep the hit if reaching the beginning or
                                   // an ambiguous base or the intv is small
                                   // enough
        if (curr->n ==
            0) { // test curr->n>0 to make sure there are no longer matches
          if (mem->n == 0 ||
              static_cast<bwtint_t>(i) + 1 < mem->a[mem->n - 1].info >> 32) { // skip contained matches
            ik = *p;
            ik.info |= (uint64_t)(i + 1) << 32;
            kv_push(bwtintv_t, *mem, ik);
          }
        } // otherwise the match is contained in another longer match
      } else if (curr->n == 0 || ok[c].x[2] != curr->a[curr->n - 1].x[2]) {
        ok[c].info = p->info;
        kv_push(bwtintv_t, *curr, ok[c]);
      }
    }
    if (curr->n == 0)
      break;
    swap = curr;
    curr = prev;
    prev = swap;
  }
  bwt_reverse_intvs(mem); // s.t. sorted by the start coordinate

  if (tmpvec == 0 || tmpvec[0] == 0)
    free(a[0].a);
  if (tmpvec == 0 || tmpvec[1] == 0)
    free(a[1].a);
  return ret;
}

int bwt_smem1_with_kmer(const bwt_t* bwt, int len, const uint8_t* q, int x,
                        int min_intv, bwtintv_t initial_interval,
                        bwtintv_v* mem, bwtintv_v* tmpvec[2]) {
  return bwt_smem1a_with_kmer(bwt, len, q, x, min_intv, 0, initial_interval,
                              mem, tmpvec);
}
} // namespace bwautils
