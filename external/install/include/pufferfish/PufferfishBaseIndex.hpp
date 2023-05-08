#ifndef _PUFFERFISHBASE_INDEX_HPP_
#define _PUFFERFISHBASE_INDEX_HPP_

#include <vector>
#include <iterator>

#include "core/range.hpp"
#include "cereal/archives/json.hpp"
#include "compact_vector/compact_vector.hpp"

#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "BooPHF.hpp"
#include "Util.hpp"
#include "PufferfishTypes.hpp"

template <typename T>
class PufferfishBaseIndex {
    private:
    // from : https://www.fluentcpp.com/2017/05/19/crtp-helper/
    T& underlying();
    T const& underlying() const;

protected:
  inline core::range<std::vector<pufferfish::util::Position>::iterator> contigRange(uint64_t contigRank) {
      auto spos = underlying().contigOffsets_[contigRank];
      auto epos = underlying().contigOffsets_[contigRank+1];
      std::vector<pufferfish::util::Position>::iterator startIt = underlying().contigTable_.begin() + spos;
      std::vector<pufferfish::util::Position>::iterator endIt = startIt + (epos - spos);
      return core::range<std::vector<pufferfish::util::Position>::iterator>(startIt, endIt);
    }

  using pos_vector_t = compact::vector<uint64_t>;
  using seq_vector_t = compact::vector<uint64_t, 2>;
  using edge_vector_t = compact::vector<uint64_t, 8>;
  using bit_vector_t = compact::vector<uint64_t, 1>;

  public:

  // Get the equivalence class ID (i.e., rank of the equivalence class)
  // for a given contig.
  pufferfish::types::EqClassID getEqClassID(uint32_t contigID);

  // Get the equivalence class label for a contig (i.e., the set of reference
  // sequences containing
  // the contig).
  const pufferfish::types::EqClassLabel& getEqClassLabel(uint32_t contigID);

  // Get the k value with which this index was built.
  uint32_t k();

  // Get the number of contigs in the underlying indexed cDBG
  uint64_t numContigs() const;

  // Get the list of reference sequences & positions corresponding to a contig
  const core::range<std::vector<pufferfish::util::Position>::iterator> refList(uint64_t contigRank);

  // Get the name of a given reference sequence
  inline const std::string& refName(uint64_t refRank) {
    return underlying().refNames_[refRank+underlying().refExt_[refRank]];
  }

  // Get the length of a reference sequence
  inline uint32_t refLength(uint64_t refRank) const {
    return underlying().refLengths_[refRank+underlying().refExt_[refRank]];
  }

  inline uint32_t completeRefLength(uint64_t refRank) const {
    return underlying().completeRefLengths_[refRank+underlying().refExt_[refRank]];
  }

  // Get the list of reference names
  const std::vector<std::string>& getFullRefNames() ;
  // and lengths
  const std::vector<uint32_t>& getFullRefLengths() const;
  // and lengths
  const std::vector<uint32_t>& getFullRefLengthsComplete () const;

  uint64_t getIndexedRefCount() const;

  inline uint64_t getRefId(uint64_t id) const {
    return id + underlying().refExt_[id];
  }

  // Returns true if pos is a valid position in the compacted sequence array
  // and false otherwise.
  bool isValidPos(uint64_t pos);

  // Returns a ProjectedHits object that contains all of the
  // projected reference hits for the given kmer.
  auto getRefPos(CanonicalKmer& mer) -> pufferfish::util::ProjectedHits;

  // Returns a copy of the string value of contig sequence vector starting from position `globalPos` with `length` bases
  // and reverse-complements the string if `isFw` is false
  std::string getSeqStr(size_t globalPos, int64_t length, bool isFw=true);

  // Returns a copy of the string value of the reference sequence vector starting from position `start` for `length` bases
  std::string getRefSeqStr(size_t start, int64_t length);

  // Returns a ProjectedHits object that contains all of the
  // projected reference hits for the given kmer.  Uses the results
  // of the previous contig info (start, end) from qc if the same
  // contig contains the match.  For correlated searches (e.g., from a read)
  // this can considerably speed up querying.
  auto getRefPos(CanonicalKmer& mer, pufferfish::util::QueryCache& qc) -> pufferfish::util::ProjectedHits;

  typename PufferfishBaseIndex<T>::seq_vector_t& getSeq(); 
  typename PufferfishBaseIndex<T>::seq_vector_t& getRefSeq(); 
  typename PufferfishBaseIndex<T>::edge_vector_t& getEdge(); 

  uint8_t getEdgeEntry(uint64_t contigRank) const;
  //uint8_t getRevEdgeEntry(uint64_t contigRank) {return revedge_[contigRank];}
  CanonicalKmer getStartKmer(uint64_t cid) ;
  CanonicalKmer getEndKmer(uint64_t cid) ;
  uint32_t getContigLen(uint64_t cid) ;
  uint64_t getGlobalPos(uint64_t cid) ;
  auto  getContigBlock(uint64_t rank) -> pufferfish::util::ContigBlock ;

  bool hasReferenceSequence() const; 

  // Returns true if the reference with the given
  // rank is a decoy, and false otherwise.
  bool isDecoy(uint64_t rank) const;
  bool isDecoyEncodedIndex(uint64_t rank) const;
  uint64_t firstDecoyIndex() const;
  uint64_t firstDecoyEncodedIndex() const;


  // Returns true if the given k-mer appears in the dBG, false otherwise
  /*
  bool contains(CanonicalKmer& mer);

  // Returns the contigID associated with the k-mer mer
  uint32_t contigID(CanonicalKmer& mer);

  // Returns the position in the compacted dBG sequence vector where the
  // given k-mer occurs, or std::numeric_limits<uint64_t>::max() otherwise.
  uint64_t getRawPos(CanonicalKmer& mer);
  */
  /*
  std::vector<CanonicalKmer> getNextKmerOnGraph(uint64_t cid, pufferfish::util::Direction dir, bool isCurContigFwd);
  */
  //void getRawSeq(pufferfish::util::ProjectedHits& phits, CanonicalKmerIterator& kit, std::string& contigStr, int readLen);
};

#endif // _PUFFERFISHBASE_INDEX_HPP_
