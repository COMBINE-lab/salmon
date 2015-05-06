#ifndef ALIGNMENT_GROUP
#define ALIGNMENT_GROUP

extern "C" {
#include "io_lib/scram.h"
#include "io_lib/os.h"
#undef max
#undef min
}


// Cereal includes
#include "cereal/types/vector.hpp"
#include "cereal/archives/binary.hpp"

#include <vector>
#include "SalmonMath.hpp"
#include "ReadPair.hpp"

template <typename FragT>
class AlignmentGroup {
    public:
        AlignmentGroup() : read_(nullptr), isUniquelyMapped_(true) { alignments_.reserve(10); }
        AlignmentGroup(AlignmentGroup& other) = delete;
        AlignmentGroup(AlignmentGroup&& other) = delete;
        AlignmentGroup& operator=(AlignmentGroup& other) = delete;
        AlignmentGroup& operator=(AlignmentGroup&& other) = delete;

        void setRead(std::string* r) { read_ = r; }
        std::string* read() { return read_; }

        inline std::vector<FragT>& alignments() { return alignments_; }
        void emplaceAlignment(FragT&& p) { alignments_.emplace_back(p); }
        void addAlignment(FragT& p) { alignments_.push_back(p); }

        void clearAlignments() {
            alignments_.clear();
            isUniquelyMapped_ = true;
        }

        inline bool& isUniquelyMapped() { return isUniquelyMapped_; }
        inline size_t numAlignments() { return alignments_.size(); }
        inline size_t size() { return numAlignments(); }

        template <typename Archive>
        void serialize(Archive& archive) {
            archive(alignments_);
        }

    private:
        std::vector<FragT> alignments_;
        std::string* read_;
        bool isUniquelyMapped_;
};
#endif // ALIGNMENT_GROUP
