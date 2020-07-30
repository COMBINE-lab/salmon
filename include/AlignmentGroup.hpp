#ifndef ALIGNMENT_GROUP
#define ALIGNMENT_GROUP

extern "C" {
#include "io_lib/os.h"
#include "io_lib/scram.h"
#undef max
#undef min
}

// Cereal includes
#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"

#include "ReadPair.hpp"
#include "SalmonMath.hpp"
#include "UtilityFunctions.hpp"
#include <vector>

struct EmptyCellInfo {};

template <typename FragT, typename BCType = EmptyCellInfo, typename UMIType = EmptyCellInfo> class AlignmentGroup {
public:
  AlignmentGroup() : read_(nullptr), isUniquelyMapped_(true) {
    alignments_.reserve(10);
  }
  AlignmentGroup(AlignmentGroup& other) = delete;
  AlignmentGroup(AlignmentGroup&& other) = default;
  AlignmentGroup& operator=(AlignmentGroup& other) = delete;
  AlignmentGroup& operator=(AlignmentGroup&& other) = delete;

  void setRead(std::string* r) { read_ = r; }
  std::string* read() { return read_; }

   void setBarcode(BCType b ) {barcode_ = b;}
   void setUMI(UMIType b) {umi_ = b;}
   BCType barcode() {return barcode_;}
   UMIType umi() {return umi_;}

  inline std::vector<FragT>& alignments() { return alignments_; }
  void emplaceAlignment(FragT&& p) { alignments_.emplace_back(p); }
  void addAlignment(FragT& p) { alignments_.push_back(p); }

  void clearAlignments() {
    auto vecLen = alignments_.size();
    alignments_.clear();
    if (vecLen > 10) {
      alignments_.shrink_to_fit();
    }
    isUniquelyMapped_ = true;
  }

  inline bool& isUniquelyMapped() { return isUniquelyMapped_; }
  inline size_t numAlignments() const { return alignments_.size(); }
  inline size_t size() const { return numAlignments(); }

  template <typename Archive> void serialize(Archive& archive) {
    archive(alignments_);
  }

  /**
   *  Sort the alignments by their transcript ids
   */
  inline void sortHits() {
    std::sort(alignments_.begin(), alignments_.end(),
              [](const FragT& x, const FragT& y) -> bool {
                return x->transcriptID() < y->transcriptID();
              });
  }

  /**
   * Sort the alignments by their transcript ids and returns in
   * primaryIndex the index in the sorted vector of the primary
   * alignment for each alignment.
   */
  void sortHits(std::vector<int>& primaryIndex, bool debug = false);

private:
  std::vector<FragT> alignments_;
  std::string* read_;
  bool isUniquelyMapped_;
  BCType barcode_;
  UMIType umi_;
};


// Implementation

template <typename FragT, typename BCType, typename UMIType>
void AlignmentGroup<FragT, BCType, UMIType>::sortHits(std::vector<int>& primaryIndex, bool debug) {
  // Create initial primaryIndex
  if(primaryIndex.size() < alignments_.size())
    primaryIndex.resize(alignments_.size(), 0);

  //  int curPrimary = 0;
  for(int i = 0; i < alignments_.size(); ++i) {
    if(alignments_[i]->isSecondary())
      primaryIndex[i] = 0;
    else
      primaryIndex[i] = i;
  }

  // Create permuation array of sorted alignments
  chobo::small_vector<int> permutation(range{0}, range{alignments_.size()});
  std::sort(permutation.begin(), permutation.end(),
            [&](const int x, const int y) -> bool {
              return alignments_[x]->transcriptID() < alignments_[y]->transcriptID();
            });
  if(debug) {
    std::ostringstream os;
    os << "perm";
    for(auto x : permutation)
      os << ' ' << x;
    os << '\n';
    os << "primary";
    for(auto x : primaryIndex)
      os << ' ' << x;
    os << '\n';
    std::cerr << os.str();
  }

  // Update primaryIndex
  chobo::small_vector<int> permutation2(alignments_.size());
  for(size_t i = 0; i < permutation.size(); ++i)
    permutation2[permutation[i]] = i;
  for(size_t i = 0; i < alignments_.size(); ++i)
    primaryIndex[i] = permutation2[primaryIndex[i]];

  if(debug) {
    std::ostringstream os;
    os << "perm2";
    for(auto x : permutation)
      os << ' ' << x;
    os << "\npermutation2";
    for(auto x : permutation2)
      os << ' ' << x;
    os << "\nprimary2";
    for(auto x : primaryIndex)
      os << ' ' << x;
    os << '\n';
    std::cerr << os.str();
  }

  // Reorder alignments in place
  permute(permutation, alignments_, primaryIndex);
}


#endif // ALIGNMENT_GROUP
