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
  inline size_t numAlignments() { return alignments_.size(); }
  inline size_t size() { return numAlignments(); }

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

private:
  std::vector<FragT> alignments_;
  std::string* read_;
  bool isUniquelyMapped_;
  BCType barcode_;
  UMIType umi_;
};

#endif // ALIGNMENT_GROUP
