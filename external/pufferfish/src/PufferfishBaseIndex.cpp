#include "PufferfishBaseIndex.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "PufferfishLossyIndex.hpp"

// from : https://www.fluentcpp.com/2017/05/19/crtp-helper/
template <typename T>
T& PufferfishBaseIndex<T>::underlying() { return static_cast<T&>(*this); }
template <typename T>
T const& PufferfishBaseIndex<T>::underlying() const { return static_cast<T const&>(*this); }

//template <typename T>
//PufferfishBaseIndex<T>::PufferfishBaseIndex() {}
//PufferfishBaseIndex::PufferfishBaseIndex(const std::string& indexDir) {}

template <typename T>
pufferfish::types::EqClassID PufferfishBaseIndex<T>::getEqClassID(uint32_t contigID) {
  return underlying().eqClassIDs_[contigID];
}

template <typename T>
const pufferfish::types::EqClassLabel&
PufferfishBaseIndex<T>::getEqClassLabel(uint32_t contigID) {
  return underlying().eqLabels_[getEqClassID(contigID)];
}

/** Seems not to be used, ignore for now **/
/*
uint64_t PufferfishBaseIndex::getRawPos(CanonicalKmer& mer) {
  auto km = mer.getCanonicalWord();
  size_t res = hash_raw_->lookup(km);
  if (res < numKmers_) {
    uint64_t pos = pos_[res];
    uint64_t fk = seq_.get_int(2 * pos, twok_);
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      return pos;
    }
  }
  return std::numeric_limits<uint64_t>::max();
}
*/

/** Seems not to be used, ignore for now **/
/*
bool PufferfishBaseIndex::contains(CanonicalKmer& mer) {
  return isValidPos(getRawPos(mer));
}
*/

/** Seems not to be used, ignore for now **/
/*
uint32_t PufferfishIndex::contigID(CanonicalKmer& mer) {
  auto km = mer.getCanonicalWord();
  size_t res = hash_raw_->lookup(km);
  if (res < numKmers_) {
    uint64_t pos = pos_[res];
    uint64_t fk = seq_.get_int(2 * pos, twok_);
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      auto rank = contigRank_(pos);
      return rank;
    }
  }
  return std::numeric_limits<uint32_t>::max();
}
*/

template <typename T>
bool PufferfishBaseIndex<T>::isValidPos(uint64_t pos) {
  return pos != std::numeric_limits<uint64_t>::max();
}

template <typename T>
std::string PufferfishBaseIndex<T>::getRefSeqStr(size_t start, int64_t length) {
  if (length == 0) { return ""; }
  auto& rseq_ = underlying().refseq_;
  std::string seq(length, 'N');
  uint64_t c = 0;
  uint64_t bucket_offset = (start)*2;
  auto len_on_vector = length * 2;
  int32_t toFetch = len_on_vector;
  while (toFetch > 0) {
    uint32_t len = (toFetch >= 64) ? 64 : toFetch;
    toFetch -= len;
    uint64_t word = rseq_.get_int(bucket_offset, len);
    for (uint32_t i = 0; i < len; i += 2) {
      uint8_t next_bits = ((word >> i) & 0x03);
      seq[c++] = "ACGT"[next_bits];
    }
    bucket_offset += len;
  }
  return seq;
}

template <typename T>
std::string PufferfishBaseIndex<T>::getSeqStr(size_t start, int64_t length, bool isFw) {
  if (length <= 0) { return ""; }
  auto& seq_ = underlying().seq_;
  std::string outstr(length, 'N');
  // length - 1 must be >= 0 because we return above if length <= 0.
  uint64_t c = isFw ? 0 : static_cast<uint64_t>(length - 1);

  uint64_t bucket_offset = (start)*2;
  auto len_on_vector = length * 2;
  int32_t toFetch = len_on_vector;
  while (toFetch > 0) {
    uint32_t len = (toFetch >= 64) ? 64 : toFetch;
    toFetch -= len;
    uint64_t word = seq_.get_int(bucket_offset, len);
    if (isFw) {
      for (uint32_t i = 0; i < len; i += 2) {
        uint8_t next_bits = ((word >> i) & 0x03);
        outstr[c++] = "ACGT"[next_bits];
      }
    } else {
      for (uint32_t i = 0; i < len; i += 2) {
        uint8_t next_bits = ((word >> i) & 0x03);
        outstr[c++] = "TGCA"[next_bits];
      }
    }
    bucket_offset += len;
  }
  // we can simply fill this in in reverse, but I don't feel like 
  // checking the arithmetic (bounds) now, and, for the time being,
  // this function should never actually be called.
  if (!isFw) {
    std::reverse(outstr.begin(), outstr.end());
  }
  return outstr;
}


/**
 * Returns a ProjectedHits object containing all of the reference loci matching this
 * provided Canonical kmer (including the oritentation of the match).  The provided
 * QueryCache argument will be used to avoid redundant rank / select operations if feasible.
 */

template <typename T>
auto PufferfishBaseIndex<T>::getRefPos(CanonicalKmer& mer, pufferfish::util::QueryCache& qc)
    -> pufferfish::util::ProjectedHits {
    return underlying().getRefPos(mer, qc);
    }

template <typename T>
auto PufferfishBaseIndex<T>::getRefPos(CanonicalKmer& mer) -> pufferfish::util::ProjectedHits {
    return underlying().getRefPos(mer);
}

template <typename T>
uint32_t PufferfishBaseIndex<T>::k() { return underlying().k_; }


template <typename T>
CanonicalKmer PufferfishBaseIndex<T>::getStartKmer(uint64_t rank){
  T& derived = underlying();
  auto& rankSelDict = derived.rankSelDict;
  auto& seq_ = derived.seq_;
  auto k_ = derived.k_;
  CanonicalKmer::k(k_) ;
  CanonicalKmer kb ;
  uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict.select(rank)) + 1;
  uint64_t fk = seq_.get_int(2*sp, 2*k_) ;
  kb.fromNum(fk) ;
  return kb ;
}

template <typename T>
CanonicalKmer PufferfishBaseIndex<T>::getEndKmer(uint64_t rank){
  T& derived = underlying();
  auto& rankSelDict = derived.rankSelDict;
  auto& seq_ = derived.seq_;
  auto k_ = derived.k_;
  CanonicalKmer::k(k_) ;
  CanonicalKmer kb ;
  uint64_t contigEnd = rankSelDict.select(rank + 1);

  uint64_t fk = seq_.get_int(2*(contigEnd - k_ + 1), 2*k_) ;
  kb.fromNum(fk) ;
  return kb ;
}

/** Seems not to be used, so skip for now **/
/*
std::vector<CanonicalKmer> PufferfishIndex::getNextKmerOnGraph(uint64_t rank, pufferfish::util::Direction dir, bool isCurContigFwd){
  //get the edge vec
  std::vector<CanonicalKmer> nextKmers ;
  uint8_t edgeVec = edge_[rank] ;
  uint8_t mask = 1 ;
  std::vector<char> nuclmap = {'C','G','T','A','C','G','T','A'} ;
  std::map<char, char> cMap = {{'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}} ;

  if(dir == pufferfish::util::Direction::FORWARD){
    // We need to append so let's concentrate on the lower 4 bits
    auto ke = getEndKmer(rank) ;
    auto ktmp = ke ;
    for(uint8_t i=0; i < 4; ++i){
      ktmp = ke ;
      if(edgeVec & (mask << i)){
        char c = nuclmap[i] ;
        char charToAdd = (isCurContigFwd) ? c : cMap[c] ;
        ktmp.shiftFw(charToAdd) ;
        nextKmers.push_back(ktmp) ;
      }
    }
  }else{
    auto kb = getStartKmer(rank) ;
    auto ktmp = kb ;
    for(uint8_t i=4; i < 8; ++i){
      ktmp = kb ;
      if(edgeVec & (mask << i)){
        char c = nuclmap[i] ;
        char charToAdd = (isCurContigFwd) ? c : cMap[c] ;
        ktmp.shiftBw(charToAdd) ;
        nextKmers.push_back(ktmp) ;
      }
    }
  }
  return nextKmers ;
}
*/

template <typename T>
uint32_t PufferfishBaseIndex<T>::getContigLen(uint64_t rank){
  T& derived = underlying();
  auto& rankSelDict = derived.rankSelDict;
  uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict.select(rank)) + 1;
  uint64_t contigEnd = rankSelDict.select(rank + 1);
  return (static_cast<uint32_t>(contigEnd - sp + 1)) ;
}

template <typename T>
uint64_t PufferfishBaseIndex<T>::getGlobalPos(uint64_t rank){
  T& derived = underlying();
  auto& rankSelDict = derived.rankSelDict;
  uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict.select(rank)) + 1;
  return sp ;
}

template <typename T>
auto  PufferfishBaseIndex<T>::getContigBlock(uint64_t rank)->pufferfish::util::ContigBlock{
  T& derived = underlying();
  auto& rankSelDict = derived.rankSelDict;
  auto& seq_ = derived.seq_;
  auto  k_ = derived.k_;
  CanonicalKmer::k(k_) ;
  CanonicalKmer kb;
  CanonicalKmer ke;

  uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict.select(rank)) + 1;
  uint64_t contigEnd = rankSelDict.select(rank+1) ;

  uint32_t clen = static_cast<uint32_t>(contigEnd - sp + 1) ;
  uint64_t fk = seq_.get_int(2*sp, 2*k_) ;
  kb.fromNum(fk) ;

  fk = seq_.get_int(2*(contigEnd - k_ + 1), 2*k_) ;
  ke.fromNum(fk) ;

  std::string seq = getSeqStr(sp,clen) ;

  return {rank,sp,clen,seq};
}

// Get the number of contigs in the underlying indexed cDBG
template <typename T>
uint64_t PufferfishBaseIndex<T>::numContigs() const {
  return underlying().numContigs_;
}

/**
 * Return the position list (ref_id, pos) corresponding to a contig.
 */
template <typename T>
const core::range<std::vector<pufferfish::util::Position>::iterator>
PufferfishBaseIndex<T>::refList(uint64_t contigRank) {
  return contigRange(contigRank);
}


template <typename T>
bool PufferfishBaseIndex<T>::hasReferenceSequence() const {
  return underlying().haveRefSeq_;
}

template <typename T>
const std::vector<std::string>& PufferfishBaseIndex<T>::getFullRefNames() {
  return underlying().refNames_;
}

template <typename T>
const std::vector<uint32_t>& PufferfishBaseIndex<T>::getFullRefLengths() const {
  return underlying().refLengths_;
}

template <typename T>
const std::vector<uint32_t>& PufferfishBaseIndex<T>::getFullRefLengthsComplete() const {
  return underlying().completeRefLengths_;
}

template <typename T>
uint64_t PufferfishBaseIndex<T>::getIndexedRefCount() const {
  return underlying().refExt_.size();
}

template <typename T>
uint8_t PufferfishBaseIndex<T>::getEdgeEntry(uint64_t contigRank) const {
    return underlying().edge_[contigRank];
}

template <typename T>
typename PufferfishBaseIndex<T>::seq_vector_t& PufferfishBaseIndex<T>::getRefSeq() {
  return underlying().refseq_;
}

template <typename T>
typename PufferfishBaseIndex<T>::seq_vector_t& PufferfishBaseIndex<T>::getSeq() {
  return underlying().seq_;
}

template <typename T>
typename PufferfishBaseIndex<T>::edge_vector_t& PufferfishBaseIndex<T>::getEdge() {
  return underlying().edge_;
}

template <typename T>
bool PufferfishBaseIndex<T>::isDecoy(uint64_t rank) const {
  return (underlying().numDecoys_ > 0) ? (rank >= underlying().firstDecoyIndex_) : false;
}

template <typename T>
bool PufferfishBaseIndex<T>::isDecoyEncodedIndex(uint64_t rank) const {
  return (underlying().numDecoys_ > 0) ? (rank >= underlying().firstDecoyEncodedIndex_) : false;
}

template <typename T>
uint64_t PufferfishBaseIndex<T>::firstDecoyIndex() const { return underlying().firstDecoyIndex_ ;}

template <typename T>
uint64_t PufferfishBaseIndex<T>::firstDecoyEncodedIndex() const { return underlying().firstDecoyEncodedIndex_ ;}

template class  PufferfishBaseIndex<PufferfishIndex>;
template class  PufferfishBaseIndex<PufferfishSparseIndex>;
template class  PufferfishBaseIndex<PufferfishLossyIndex>;
