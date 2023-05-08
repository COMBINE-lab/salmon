#ifndef _PUFFERFISHSPARSE_INDEX_HPP_
#define _PUFFERFISHSPARSE_INDEX_HPP_

#include <vector>

#include "rank9b.hpp"
#include "core/range.hpp"
#include "cereal/archives/json.hpp"

#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "BooPHF.hpp"
#include "Util.hpp"
#include "PufferfishBaseIndex.hpp"
#include "rank9sel.hpp"

class PufferfishSparseIndex : public PufferfishBaseIndex<PufferfishSparseIndex> {
  friend PufferfishBaseIndex;
  using hasher_t = pufferfish::types::hasher_t;
  using boophf_t = pufferfish::types::boophf_t;
  using pos_vector_t = PufferfishBaseIndex<PufferfishSparseIndex>::pos_vector_t;
  using seq_vector_t = PufferfishBaseIndex<PufferfishSparseIndex>::seq_vector_t;
  using edge_vector_t = PufferfishBaseIndex<PufferfishSparseIndex>::edge_vector_t;
  using bit_vector_t = PufferfishBaseIndex<PufferfishSparseIndex>::bit_vector_t;

private:
  uint32_t k_{0};
  uint32_t twok_{0};
  int32_t extensionSize_{0};
  uint64_t numKmers_{0};
  uint64_t lastSeqPos_{0};
  uint64_t numSampledKmers_{0};
  bool haveEdges_{false};
  bool haveRefSeq_{false};
  bool haveEqClasses_{true};

  std::vector<uint32_t> eqClassIDs_;
  std::vector<std::vector<uint32_t>> eqLabels_;
  std::vector<std::string> refNames_;
  std::vector<uint32_t> refLengths_;
  std::vector<uint32_t> completeRefLengths_;
  std::vector<uint32_t> refExt_;
  std::vector<pufferfish::util::Position> contigTable_;
  compact::vector<uint64_t> contigOffsets_{16};

  uint64_t numContigs_{0};
  bit_vector_t contigBoundary_;
  rank9sel rankSelDict;

  seq_vector_t seq_;
  edge_vector_t edge_;
  //for sparse representation
  bit_vector_t presenceVec_;
  bit_vector_t canonicalNess_;
  bit_vector_t directionVec_ ;
  rank9b presenceRank_;
  compact::vector<uint64_t> extSize_{16};
  compact::vector<uint64_t> auxInfo_{16};
  pos_vector_t sampledPos_{16};

  std::unique_ptr<boophf_t> hash_{nullptr};
  uint64_t numDecoys_{0};
  uint64_t firstDecoyIndex_{0};
  uint64_t firstDecoyEncodedIndex_{0};

  static const constexpr uint64_t shiftTable_[] = {
    0x0, 0x7, 0x38, 0x1c0, 0xe00, 0x7000, 0x38000, 0x1c0000,
    0xe00000, 0x7000000, 0x38000000, 0x1c0000000, 0xe00000000,
    0x7000000000, 0x38000000000, 0x1c0000000000, 0xe00000000000,
    0x7000000000000, 0x38000000000000, 0x1c0000000000000,
    0xe00000000000000, 0x7000000000000000};

public:
  compact::vector<uint64_t, 2> refseq_;
  std::vector<uint64_t> refAccumLengths_;


  PufferfishSparseIndex();
  PufferfishSparseIndex(const std::string& indexPath, pufferfish::util::IndexLoadingOpts opts = pufferfish::util::IndexLoadingOpts());
  //Returns the position in the compacted bBG sequence from the sparse
  //index the above routine can be replaced by this code in
  //future versions
  // uint64_t getSparseRawPos(CanonicalKmer& mer);

  // Returns a ProjectedHits object that contains all of the
  // projected reference hits for the given kmer.
  auto getRefPos(CanonicalKmer mer) -> pufferfish::util::ProjectedHits;
  auto getRefPos(CanonicalKmer mer, pufferfish::util::QueryCache& qc) -> pufferfish::util::ProjectedHits;

private:
  auto getRefPosHelper_(CanonicalKmer& mer, uint64_t pos, bool didWalk = false) -> pufferfish::util::ProjectedHits;
  auto getRefPosHelper_(CanonicalKmer& mer, uint64_t pos, pufferfish::util::QueryCache& qc, bool didWalk = false) -> pufferfish::util::ProjectedHits;

};

#endif // _PUFFERFISH_INDEX_HPP_
