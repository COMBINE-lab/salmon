#include "MemCollector.hpp"
#include <cmath>
#include <limits>

// using spp:sparse_hash_map;

// polute the namespace --- put this in the functions that need it.
namespace kmers = combinelib::kmers;

template <typename PufferfishIndexT>
void MemCollector<PufferfishIndexT>::configureMemClusterer(uint32_t max) {
  mc.setMaxAllowedRefsPerHit(max);
}

template <typename PufferfishIndexT>
void MemCollector<PufferfishIndexT>::setHitFilterPolicy(pufferfish::util::HitFilterPolicy hfp) {
  mc.setHitFilterPolicy(hfp);
}

template <typename PufferfishIndexT>
pufferfish::util::HitFilterPolicy MemCollector<PufferfishIndexT>::getHitFilterPolicy() const {
  return mc.getHitFilterPolicy();
}


template <typename PufferfishIndexT>
size_t MemCollector<PufferfishIndexT>::expandHitEfficient(pufferfish::util::ProjectedHits& hit,
                       pufferfish::CanonicalKmerIterator& kit, 
                       ExpansionTerminationType& et) {

  auto& allContigs = pfi_->getSeq();
  // startPos points to the next kmer in contig (which can be the left or
  // right based on the orientation of match)
  size_t cStartPos =
    hit.globalPos_ - hit.contigPos_   ; // next kmer in the read
  size_t cEndPos = cStartPos + hit.contigLen_;
  size_t cCurrPos = hit.globalPos_; // start from next character if fw match
  if (hit.contigOrientation_) { // if match is fw, go to the next k-mer in the
                // contig
  cCurrPos += k;
  }
  int currReadStart = kit->second + 1;
  auto readSeqView = kit.seq();
  auto readSeqLen = readSeqView.size();
  auto readSeqStart = currReadStart;
  auto readSeqOffset = currReadStart + k - 1;
  int fastNextReadCode{0};
  bool stillMatch = true;
  bool foundTermCondition = false;

  while (stillMatch and
     (cCurrPos < cEndPos) and
     (cCurrPos > cStartPos) and
     readSeqOffset < readSeqLen) { 

    if (hit.contigOrientation_) { // if fw match, compare read last base with
                  // contig first base and move fw in the
                  // contig
      auto baseCnt = k < cEndPos - cCurrPos ? k : cEndPos - cCurrPos;
      uint64_t fk = allContigs.get_int(2*cCurrPos, 2*baseCnt);
      cCurrPos += baseCnt;
      for (size_t i = 0; i < baseCnt && readSeqOffset < readSeqLen; i++) {
        // be dirty and peek into the underlying read
        fastNextReadCode = kmers::codeForChar(readSeqView[readSeqOffset]);
        int contigCode = (fk >> (2 * i)) & 0x3;
        if (fastNextReadCode != contigCode) {
          stillMatch = false;
          et = ExpansionTerminationType::MISMATCH;
          foundTermCondition = true;
          break;
        }
        hit.k_++;
        readSeqOffset++;
        readSeqStart++;
      }
    } else { // if rc match, compare read last base with contig last base and
       // move backward in the contig
      auto baseCnt = k < cCurrPos - cStartPos ? k : cCurrPos - cStartPos;
      uint64_t fk = allContigs.get_int(2*(cCurrPos - baseCnt), 2*baseCnt);
      cCurrPos -= baseCnt;
      for (int i = baseCnt - 1; i >= 0 && readSeqOffset < readSeqLen; i--) {
        // be dirty and peek into the underlying read
        fastNextReadCode = kmers::codeForChar(kmers::complement(readSeqView[readSeqOffset]));
        int contigCode = (fk >> (2 * i)) & 0x3;
        if (fastNextReadCode != contigCode) {
          stillMatch = false;
          et = ExpansionTerminationType::MISMATCH;
          foundTermCondition = true;
          break;
        }
        hit.k_++;
        readSeqOffset++;
        readSeqStart++;
      }
    }
  }

  if (!foundTermCondition) {
    et = (cCurrPos >= cEndPos or cCurrPos <= cStartPos) ? 
    ExpansionTerminationType::CONTIG_END : ExpansionTerminationType::READ_END;
  }
  if (!hit.contigOrientation_) {
    hit.contigPos_ -= (hit.k_ - k);
    hit.globalPos_ -= (hit.k_ - k);
  }
  kit.jumpTo(readSeqStart);
  return currReadStart;
}

template <typename PufferfishIndexT>
void MemCollector<PufferfishIndexT>::setConsensusFraction(double cf) {
  mc.setConsensusFraction(cf);
}

template <typename PufferfishIndexT>
double MemCollector<PufferfishIndexT>::getConsensusFraction() const {
  return mc.getConsensusFraction();
}

template <typename PufferfishIndexT>
void MemCollector<PufferfishIndexT>::setAltSkip(uint32_t as) {
  altSkip = as;
}

  
enum class LastSkipType : uint8_t { NO_HIT=0, SKIP_READ=1, SKIP_UNI=2 };

struct FastHitInfo {
  // return true if it is currently valid to do 
  // a fast check, and false otherwise.
  inline bool valid() { return fast_check; }

  // set the state of this FastHitInfo instance
  // true if it is valid to check upon next query
  // and false otherwise.
  inline void valid(bool v) { fast_check = v; }

  bool fast_check{false};
  int32_t offset{0};
  uint64_t ref_kmer{0};
};

// Idea is to move the logic of the search into here.
// We should also consider the optimizations that can 
// be done here (like having small checks expected to)
// be on the same contig bypass a hash lookup.
template <typename PufferfishIndexT>
struct SkipContext {

  SkipContext(std::string& read, PufferfishIndexT* pfi_in, int32_t k_in) : 
    kit1(read), kit_tmp(read), pfi(pfi_in), read_len(static_cast<int32_t>(read.length())),
    read_target_pos(0), read_current_pos(0), read_prev_pos(0), safe_skip(1),
    k(k_in), expected_cid(invalid_cid), last_skip_type(LastSkipType::NO_HIT), miss_it(0),
    global_contig_pos(-1) { }
  
  inline bool is_exhausted() {
    return kit1 == kit_end;
  }

  inline CanonicalKmer& curr_kmer() {
    return kit1->first;
  }

  inline bool query_kmer(pufferfish::util::QueryCache& qc) {
    bool found_match = false;
    if (fast_hit.valid()) {
      auto keq = kit1->first.isEquivalent(fast_hit.ref_kmer);
      if (keq != KmerMatchType::NO_MATCH) {
        found_match = true;
        // how the k-mer hits the contig (true if k-mer in fwd orientation,
        // false otherwise)
        bool hit_fw = (keq == KmerMatchType::IDENTITY_MATCH);
        phits.contigOrientation_ = hit_fw;
        phits.globalPos_ += fast_hit.offset;
        phits.contigPos_ += fast_hit.offset;
      }
      fast_hit.valid(false);
    }

    if (!found_match) {
      phits = pfi->getRefPos(kit1->first, qc);
    }

    return !phits.empty();
  }

  // Returns true if the current hit occurred 
  // on a unitig other than that which was expected.
  // If there is no expecation about the unitig 
  // (i.e. if `clear_expectation()` was called before
  // obtaining this hit) or if the current unitig 
  // matches expectation, then it returns false.
  inline bool hit_is_unexpected() {
    return (expected_cid != invalid_cid and phits.contigIdx_ != expected_cid);
  }

  // Returns:
  // * True if there is a target position (end of unitig / read)
  //   jump that is _beyond_ the current hit's position on the read.
  // * False if the current hit is at the position to which 
  //   the jump was intended.
  // This value can be used to determine the appropriate courese
  // of action during safe skipping.
  inline bool hit_is_before_target_pos() {
    return (miss_it > 0 and read_target_pos > kit1->second);
  }

  // Returns the current position of the skip context 
  // on the read.
  inline int32_t read_pos() { return kit1->second; }

  // Returns the most recent projected hits object 
  // obtained by this skip context.
  inline pufferfish::util::ProjectedHits proj_hits() { return phits; }

  // Clears out any expectation we have about the unitig 
  // on which the next hit should occur.
  inline void clear_expectation() { expected_cid = invalid_cid; }

  // Computes the amount by which we will move along 
  // the read and query the reference until we hit the 
  // original target position.  This skip is used in the 
  // case that we query and find that we have arrived 
  // at an unexpected unitig (i.e. JE(HD) ).
  inline int32_t compute_safe_skip() {

    // read_target_pos = //hit_is_before_target_pos() ? read_target_pos : read_target_pos-1;

    // There are 2 cases we may encounter here.
    // Either:
    // (1) The current hit was found at the original
    // jump location.  In that case we got 
    // an unexpected hit. To avoid iterating to that
    // already queried hit again, the target pos
    // will be the position before that hit
    // (2) The current hit was found at some position
    // before the original jump location.  This can 
    // happen if we have an intervening miss (i.e. JE(M))
    // before we find the discordant hit (i.e. JE(HD)).
    // In that case, we leave the target position 
    // the same.

    // In either case, we have *already* checked the k-mer 
    // at read_target_pos, and it was either a hit (1) or a miss (2).
    // So, we would like to avoid doing that check again.
    // Thus, we traverse up to position (read_target_pos - 1)

    int32_t safe_skip_end = read_target_pos - 1;
    read_prev_pos = kit_tmp->second;
    int32_t dist = (safe_skip_end) - read_prev_pos;
    
    safe_skip = 1;
    // if the distance is large enough
    if (dist > 2) {
      safe_skip = (safe_skip_end - read_prev_pos) / 2;
      safe_skip = static_cast<int32_t>(safe_skip < 1 ? 1 : safe_skip);
    }
    return dist;
  }

  // Backup kit1 to the position of kit_tmp, the backup
  // position on the read.  Save the current location of 
  // kit1 in kit_swap so we can restore it later.
  inline void rewind() {
    // start from the next position on the read and 
    // walk until we hit the read end or the position 
    // that we wanted to skip to
    kit_swap = kit1;
    kit1 = kit_tmp;
    kit_tmp = kit_swap;
  }

  // Advance kit1 and kit_tmp to (at least) the position following
  // read_target_pos.  This function is intended to be called at 
  // the end of a `safe_skip` walk when the query at the original
  // read_target_pos yielded a miss.
  // Side effect : Resets the miss_it counter to 0.
  inline void advance_to_target_pos_successor() {
    int32_t diff = (read_target_pos - kit1->second) + 1;
    kit1 += diff;
    kit_tmp = kit1;
  }

  // Reset the miss counter to 0.
  inline void clear_miss_counter() { miss_it = 0; }

  // Reset kit1 to the place it was at when 
  // `rewind()` was called.
  inline void fast_forward() {
    // NOTE: normally we would assume that 
    // kit1->second >= kit_tmp->second
    // but, this be violated if we skipped farther than we asked in the iterator
    // because there were 'N's in the read.
    // in that case case, ignore the fact that fast forward 
    // could be a rewind.  We always want to skip back to 
    // where we were before.  However, in this case we 
    // ASSUME that no hit *after* kit_tmp has been added
    // to the set of collected hits.
    kit1 = kit_tmp;
  }

  /**
   * Advance by the safe skip amount calculated in 
   * `compute_safe_skip`.  
   * Returns true if we are still before or at the expected
   * end point, and false otherwise.
   */
  inline bool advance_safe() {
    int32_t safe_skip_end = read_target_pos - 1;
    int32_t remaining_len = safe_skip_end - kit1->second;
    safe_skip = (remaining_len <= safe_skip) ? remaining_len : safe_skip;
    safe_skip = (safe_skip < 1) ? 1 : safe_skip;
    kit1 += safe_skip;
    return (kit1->second <= safe_skip_end) and (kit1 != kit_end);
  }
  
  /**
   * Given that we have obtained a hit in response to a query,
   * and given that that hit was not unexpected (i.e. the jump 
   * that lead to it was not of type JE(HD)), then this function
   * will:
   * (1) calculate the next position to be checked.
   * (2) advance kit1 to that position
   * (3) preserve the current hit iterator state in kit_tmp
   * (4) set up the appropriate context for a fast hit check if appropriate
   */
  inline void advance_from_hit() {
      int32_t skip = 1;
      // the offset of the hit on the read
      int32_t read_offset = kit1->second;
      // the skip that would take us to the last base of the read
      int32_t read_skip = (read_len - k) - read_offset;

      // If we got here after a miss, and we have an expectation, and this hit 
      // matches it, then we want to avoid replaying this whole scenario again.
      // Because this hit is on the same unitig as the one that eventually led 
      // us here, we will compute the same attempted skip.  In this case, we 
      // found what we expected, but after 1 or more misses.  In that case, 
      // we're satisfied with the rest of this unitig so we simply move on 
      // to the position after our original skip position.
      if ( (miss_it > 0) and (expected_cid != invalid_cid) and (expected_cid == phits.contigIdx_) ) {
        skip = read_target_pos - read_offset + 1;
        skip = (skip < 1) ? 1 : skip;
        kit1 += skip;
        kit_tmp = kit1;
        clear_miss_counter();
        clear_expectation();
        return;
      }

      // any valid skip should always be at least 1 base
      read_skip = (read_skip < 1) ? 1 : read_skip;
      
      size_t cStartPos = phits.globalPos_ - phits.contigPos_; 
      size_t cEndPos = cStartPos + phits.contigLen_;
      size_t cCurrPos = phits.globalPos_; 
      global_contig_pos = static_cast<int64_t>(cCurrPos);
      int32_t ctg_skip = 1;
      // fw ori
      if (phits.contigOrientation_) {
        ctg_skip = static_cast<int64_t>(cEndPos) - (static_cast<int64_t>(cCurrPos + k));
      } else { // rc ori
        ctg_skip = static_cast<int32_t>(phits.contigPos_);
      }
      // we're already at the end of the contig
      bool at_contig_end = (ctg_skip == 0);
      // if we're at the end of the contig
      // we'll set the contig skip to be one 
      // (i.e. look for the next k-mer), but we won't
      // set an expectation on what that contig should be
      if (at_contig_end) {
        ctg_skip = 1;
        clear_expectation();
      } else { // otherwise, we expect the next contig to be the same
        expected_cid = phits.contigIdx_;
      }

      // remember where we are coming from
      kit_tmp = kit1;
      
      // remember the reason we are skipping
      last_skip_type = (read_skip <= ctg_skip) ? LastSkipType::SKIP_READ : LastSkipType::SKIP_UNI;
      
      // skip will be min of read and contig
      skip = (last_skip_type == LastSkipType::SKIP_READ) ? read_skip : ctg_skip;

      // skip must be at least 1
      skip = skip < 1 ? 1 : skip;
      
      // onward
      int32_t pos_before_skip = kit1->second;
      kit1 += skip;

      read_target_pos = (kit1->second > read_len - k) ? read_len - k : kit1->second;
      clear_miss_counter();

      // was the skip we got the one we expected?
      // that is, was the intervening part of the read 
      // free of `N` characters?
      bool expected_skip = (kit1->second - pos_before_skip) == skip;
      
      // if we didn't get the expected skip, then un-set our expectation
      // about where we will land.  
      // NOTE : This should be relatively infrequent, but check the effect
      if (!expected_skip) { clear_expectation(); }

      // if we got the skip we expected, then 
      // set ourselves up for a fast check in case we see 
      // what we expect to see.
      if (expected_skip and (expected_cid != invalid_cid) and (kit1 != kit_end)) {
        auto& allContigs = pfi->getSeq();
        if (phits.contigOrientation_) { 
          // if match is fw, go to the next k-mer in the contig
          cCurrPos += skip;
          /*
          if (cCurrPos + k > cEndPos) {
            std::cerr << cCurrPos << " + " << k << " > " << cEndPos << " : shouldn't happen!\n"; 
            fast_hit.valid(false);
          } else {
          */
          fast_hit.valid(true);
          fast_hit.offset = skip;
          fast_hit.ref_kmer = allContigs.get_int(2*cCurrPos, 2*k);
          //}
          return;
        } else {
          cCurrPos -= skip;
          /*
          if (cCurrPos < cStartPos ) {
            std::cerr << cCurrPos << " - " << k << " < " << cStartPos << " : shouldn't happen!\n"; 
            fast_hit.valid(false);
          } else {
          */
          fast_hit.valid(true);
          fast_hit.offset = -skip;
          fast_hit.ref_kmer = allContigs.get_int(2*cCurrPos, 2*k);
          //}
          return;
        }
      }
  }

  inline void advance_from_miss() {
      int32_t skip = 1;

      // distance from backup position 
      int32_t dist = read_target_pos - kit_tmp->second;

      switch (last_skip_type) {
        // we could have not yet seen a hit
        // should move alt_skip at a time
        case LastSkipType::NO_HIT : {
          //int32_t dist_to_end = (read_len - (kit1->second + k));
          kit1 += 5;
          return;
        }
        break;

        // we could have seen a hit, and tried to jump to the end of the read / uni
        // and the previous search was either that hit, or a miss
        case LastSkipType::SKIP_READ  :
        case LastSkipType::SKIP_UNI : {
          // if distance from backup is 1, that means we already tried
          // last position, so we should just move to the next position 
          int32_t pos_before_skip = kit_tmp->second;
          if (dist <= 1) {
            // if we're past the position we tried to jump to
            // then we're willing to skip a little faster
            skip = (dist < 0) ? 2 : 1;
            kit_tmp += skip;
            kit1 = kit_tmp;
          } else {
            // otherwise move the backup position toward us
            // and move the current point to the backup
            skip = (miss_it < 4) ? dist / 2 : 2;
            skip = (skip < 2) ? 2 : skip;
            kit_tmp += skip;
            kit1 = kit_tmp;
          }

          // If we can advance to the next query 
          // point without having to consult the index
          // then do it.
          if (// if we have an expecation
              (expected_cid != invalid_cid) and 
              // and the prev hit was on the expected unitig
              (phits.contigIdx_ == expected_cid) and 
              // and we are still before the final target position 
              (kit1->second < read_target_pos)) {
            
            // Here, we compute the actual amount we skipped by, since
            // there may have been intervening 'N's in the read and 
            // we could have skipped more than the expected amount.
            int actual_skip = (kit1->second - pos_before_skip);

            // If this was the first miss we encountered 
            // in the most recent round of search, then 
            // reset the fast_hit offset.
            if (miss_it == 0) { fast_hit.offset = 0; }
            
            auto& allContigs = pfi->getSeq();
            if (phits.contigOrientation_) {
              fast_hit.offset += actual_skip;
              global_contig_pos += actual_skip;
            } else {
              fast_hit.offset -= actual_skip;
              global_contig_pos -= actual_skip;
            }
            fast_hit.valid(true);
            fast_hit.ref_kmer = allContigs.get_int(2*global_contig_pos, 2*k);
          }
          // if we pass the read target position, then 
          // we no longer have an expectation of what 
          // we should see.
          if (kit1->second >= read_target_pos) { 
            //fast_hit.valid(false);
            clear_expectation(); 
          }
          // increment the miss iterator
          miss_it += 1;
          return;
        }
        break;
      }
      // we could have previously not seen a hit either
      // that doesn't matter since it is recursively one of the above cases.
  }

  pufferfish::CanonicalKmerIterator kit1;
  pufferfish::CanonicalKmerIterator kit_tmp;
  pufferfish::CanonicalKmerIterator kit_end;
  pufferfish::CanonicalKmerIterator kit_swap;
  PufferfishIndexT* pfi={nullptr};
  int32_t read_len;
  int32_t read_target_pos;
  int32_t read_current_pos;
  int32_t read_prev_pos;
  int32_t safe_skip;
  int32_t k;
  FastHitInfo fast_hit;
  uint32_t expected_cid;
  LastSkipType last_skip_type{LastSkipType::NO_HIT};
  int32_t miss_it;
  int64_t global_contig_pos;
  pufferfish::util::ProjectedHits phits;
  static constexpr uint32_t invalid_cid{std::numeric_limits<uint32_t>::max()};
};

// This method performs k-mer / hit collection 
// using a custom implementation of the corresponding 
// part of the pseudoalignment algorithm as described in (1).
// Specifically, it attempts to find a small set of 
// k-mers along the fragment (`read`) that are shared with
// a set of unitigs in the compacted colored de Bruijn graph, 
// using the structure of the graph to avoid queries that
// are unlikely to change the resulting set of unitigs 
// that are discovered.  This function only fills out the 
// set of hits, and does not itself implement any 
// consensus mechanism.  This hit collection mechanism 
// prioritizes speed compared to e.g. the uniMEM collection
// strategy implemented in the `operator()` method of this 
// class.  Currently, this strategy is only used in the 
// `--sketch` mode of alevin.  One may refer to the 
// [release notes](https://github.com/COMBINE-lab/salmon/releases/tag/v1.4.0)
// of the relevant version of salmon and links therein 
// for a more detailed discussion of the downstream effects of 
// different strategies.
//
// The function fills out the appropriate raw hits 
// member of this MemCollector instance.  
// 
// If the `isLeft` flag is set to true, then the left 
// raw hits are filled and can be retreived with 
// `get_left_hits()`.  Otherwise, the right raw hits are 
// filled and can be retrieved with `get_right_hits()`.
//
// This function returns `true` if at least one hit was 
// found shared between the read and reference and 
// `false` otherwise.
// 
// [1] Bray NL, Pimentel H, Melsted P, Pachter L. 
// Near-optimal probabilistic RNA-seq quantification. 
// Nat Biotechnol. 2016;34(5):525-527.
//
template <typename PufferfishIndexT>
bool MemCollector<PufferfishIndexT>::get_raw_hits_sketch(std::string &read,
                  pufferfish::util::QueryCache& qc,
                  bool isLeft,
                  bool verbose) {
  (void) verbose;
  pufferfish::util::ProjectedHits phits;
  auto& raw_hits = isLeft ? left_rawHits : right_rawHits;

  CanonicalKmer::k(k);
  int32_t k = static_cast<int32_t>(CanonicalKmer::k());
  SkipContext<PufferfishIndexT> skip_ctx(read, pfi_, k);
  
  // while it is possible to search further
  while (!skip_ctx.is_exhausted()) {
    
    // if we had a hit
    if (skip_ctx.query_kmer(qc)) {

      // record this hit
      int32_t read_pos = skip_ctx.read_pos();
      auto proj_hits = skip_ctx.proj_hits();
      
      // if the hit was not inline with what we were 
      // expecting. 
      if (skip_ctx.hit_is_unexpected()){
        // there are two scenarios now. Either: 
        // (1) the hit we found was at the inital place 
        // we tried to jump to, but just on the wrong unitig. 
        // OR
        // (2) we encountered >= 1 miss in between and so 
        // the hit we found that led us here is before the 
        // original jump point.
        //
        // In case (1), we want to add the current hit 
        // *after* doing safe advances to the end position, and 
        // then do a normal `advance_from_hit()`; in case (2) 
        // we want to add the current hit *before* doing advances 
        // to the end position, and then do a normal `advance_from_miss()`.
        
        bool hit_at_end = !skip_ctx.hit_is_before_target_pos();

        // we are in case (2)
        if (!hit_at_end) {
          raw_hits.push_back(std::make_pair(read_pos, proj_hits));
        }

        int32_t dist_to_target = skip_ctx.compute_safe_skip();

        // special case is the prev hit was the 
        // one right before this, in this case
        // no need to do this work, we already 
        // have both hits.  In that case, just
        // proceed as usual.
        if (dist_to_target > 0) {

            // start from the next position on the read 
            skip_ctx.rewind();

            // walk until we hit the read end or the position
            // that we wanted to skip to, collecting the 
            // hits we find along the way.
            while (skip_ctx.advance_safe()) {
              if (skip_ctx.query_kmer(qc)) {
                raw_hits.push_back(
                    std::make_pair(skip_ctx.read_pos(), skip_ctx.proj_hits()));
              }
            }
            
            // now we examined the positions in between.
            // if we were in case (1), jump back to the valid 
            // hit that we first encountered.  If we were in 
            // case (2), advance to at least the point right
            // after the initial jump.
            if (hit_at_end) {
              skip_ctx.fast_forward();
            } else {
              skip_ctx.advance_to_target_pos_successor();
              skip_ctx.clear_miss_counter();
            }
          }

          // if we got here, then either we are done, or we 
          // reached our target position so we don't
          // have an expectation of what should come next.
          skip_ctx.clear_expectation();

          // If we were in case (1)
          if (hit_at_end) {
            // set the phits for the skip_ctx so the 
            // next jump can be properly computed
            skip_ctx.phits = proj_hits;
            // add the hit here
            raw_hits.push_back(std::make_pair(read_pos, proj_hits));
            // advance as if this was a normal hit
            skip_ctx.advance_from_hit();
          } else {
            // We were in case 2
            // now advance as if this was a normal miss.
            skip_ctx.advance_from_miss();
          }
      } else {
        // We got a hit and either we had no expectation
        // or the hit was in accordance with our expectation.
        
        // push this hit and advance
        raw_hits.push_back(std::make_pair(read_pos, proj_hits));
        skip_ctx.advance_from_hit();
      }
    } else {
      // if we got here, we looked for a match and didn't find one
      skip_ctx.advance_from_miss();
    }
  }
  
  return raw_hits.size() != 0;
}

template <typename PufferfishIndexT>
bool MemCollector<PufferfishIndexT>::operator()(std::string &read,
                  pufferfish::util::QueryCache& qc,
                  bool isLeft,
                  bool verbose) {

  // currently unused:
  // uint32_t readLen = static_cast<uint32_t>(read.length()) ;
  pufferfish::util::ProjectedHits phits;
  auto& rawHits = isLeft ? left_rawHits : right_rawHits;

  CanonicalKmer::k(k);
  pufferfish::CanonicalKmerIterator kit_end;
  pufferfish::CanonicalKmerIterator kit1(read);
  if (verbose) {
    std::cerr << "ORIGINAL READ:\n";
    std::cerr << read << "\n";
  }

  /**
   *  Testing heuristic.  If we just succesfully matched a k-mer, and extended it to a uni-MEM, then
   *  the nucleotide after the end of the match is very likely to be a sequencing error (or a variant).
   *  In this case, the next "k" k-mers will still likely overlap that error and not match.  So, instead
   *  of looking at each of them (or skipping all of them, which might decrease sensitivity), we will
   *  skip (by a Mohsen number ;P) until we are at least k bases away from the end of the valid match.
   *  Then, so as to maintain high sensitivity, we will start skipping by only 1 k-mer at a time.
   **/

  // Start off pretending we are at least k bases away from the last hit
  uint32_t skip{1};
  uint32_t homoPolymerSkip{1};
  int32_t signedK = static_cast<int32_t>(k);
  int32_t basesSinceLastHit{signedK};
  ExpansionTerminationType et {ExpansionTerminationType::MISMATCH};

  while (kit1 != kit_end) {
    auto phits = pfi_->getRefPos(kit1->first, qc);
    skip = (basesSinceLastHit >= signedK) ? 1 : altSkip;
    if (!phits.empty()) {
      // kit1 gets updated inside expandHitEfficient function
      // stamping the readPos
      // NOTE: expandHitEfficient advances kit1 by *at least* 1 base
      size_t readPosOld = kit1->second;

      if (kit1->first.is_homopolymer()) {
        kit1 += homoPolymerSkip;
        continue;
      }

      expandHitEfficient(phits, kit1, et);
			if (verbose){
			  std::cerr<<"after expansion\n";
				std::cerr<<"readPosOld:"<<readPosOld<<" kmer:"<< kit1->first.to_str() <<"\n";
			}
      rawHits.push_back(std::make_pair(readPosOld, phits));

      basesSinceLastHit = 1;
      skip = (et == ExpansionTerminationType::MISMATCH) ? altSkip : 1;
      kit1 += (skip-1);
    } else {
     basesSinceLastHit += skip;
     kit1 += skip;
    }
  }

  // To consider references this end maps to for allowing hits on the other end
  /*auto& ref_ids = isLeft ? left_refs : right_refs;
  for (auto &hit : core::range<decltype(rawHits.begin())>(rawHits.begin(), rawHits.end())) {
    auto &projHits = hit.second;
    auto &refs = projHits.refRange;
    if (refs.size() < mc.getMaxAllowedRefsPerHit()) {
      for (auto &posIt : refs) {
        auto refid = posIt.transcript_id();
        ref_ids[refid] = true;
      }
    }
  }*/
  return rawHits.size() != 0;
}

template <typename PufferfishIndexT>
void MemCollector<PufferfishIndexT>::setChainSubOptThresh(double pre_merge_chain_sub_thresh) {
  pre_merge_chain_sub_thresh_ = pre_merge_chain_sub_thresh;
  inv_pre_merge_chain_sub_thresh_ = std::nexttoward(1.0 / pre_merge_chain_sub_thresh, std::numeric_limits<long double>::infinity());
  mc.set_chain_sub_opt_thresh_(pre_merge_chain_sub_thresh_, inv_pre_merge_chain_sub_thresh_);
}


template <typename PufferfishIndexT>
double MemCollector<PufferfishIndexT>::chainSubOptThresh() const {
  return pre_merge_chain_sub_thresh_;
}

template <typename PufferfishIndexT>
bool MemCollector<PufferfishIndexT>::findChains(std::string &read,
                                                pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>>& memClusters,
                                                //phmap::flat_hash_map<size_t, std::vector<pufferfish::util::MemCluster>>& memClusters,
                  uint32_t maxSpliceGap,
                  pufferfish::util::MateStatus mateStatus,
                  bool hChain,
                  bool isLeft,
                  bool verbose) {
  (void) maxSpliceGap;
  uint64_t firstDecoyIndex = pfi_->firstDecoyIndex();
  auto& rawHits = isLeft ? left_rawHits : right_rawHits;
  // if this is the right end of a paired-end read, use memCollectionRight,
  // otherwise (left end or single end) use memCollectionLeft.
  isSingleEnd = mateStatus == pufferfish::util::MateStatus::SINGLE_END;
  auto* memCollection = (mateStatus == pufferfish::util::MateStatus::PAIRED_END_RIGHT) ?
  &memCollectionRight : &memCollectionLeft;
  if (rawHits.size() > 0) {
    auto& other_end_refs = isLeft ? right_refs : left_refs;
    trMemMap.clear();
    mc.findOptChain(rawHits, memClusters, maxSpliceGap, *memCollection, read.length(), 
                    other_end_refs, hChain, trMemMap, firstDecoyIndex, verbose);
    /*
    if (verbose) {
      std::cerr << "lets see what we have\n";
      for (auto kv : trMemMap) {
        std::cerr <<"tid:" <<  kv.first.first << "\n";
        for (auto mem : kv.second) {
          std::cerr << "--- t" << mem.tpos << " r"
              << mem.rpos << " cid:"
              << mem.memInfo->cid << " cpos: "
              << mem.memInfo->cpos << " len:"
              << mem.extendedlen << " re:"
              << mem.memInfo->readEnd << " fw:"
              << mem.isFw << "\n";
          std::cerr << read.substr(mem.rpos,mem.extendedlen) << "\n";
        }
      }
    }
    */
    return true;
  }
  return false;
}


template <typename PufferfishIndexT>
void MemCollector<PufferfishIndexT>::clear() {
  memCollectionLeft.clear();
  memCollectionRight.clear();
  trMemMap.clear();
  left_refs.clear();
  right_refs.clear();
  left_rawHits.clear();
  right_rawHits.clear();
}

template class MemCollector<PufferfishIndex>;

template class MemCollector<PufferfishSparseIndex>;

template class MemCollector<PufferfishLossyIndex>;
