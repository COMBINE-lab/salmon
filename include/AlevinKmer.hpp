#ifndef __ALEVIN_KMER_HPP__
#define __ALEVIN_KMER_HPP__

#include "rapmap/Kmer.hpp"

//code taken from
//https://raw.githubusercontent.com/COMBINE-lab/pufferfish/a51c41308142cdd0186724b0ca2bf773f5882072/include/CanonicalKmer.hpp

namespace alevin {
  namespace kmer {
    using Mer = combinelib::kmers::Kmer<32, 2>;

    // NO_MATCH => two k-mers k1, and k2 are distinct such that k1 != k2 and rc(k1)
    // != k2
    // IDENTITY_MATCH => k1 = k2
    enum class KmerMatchType : uint8_t { NO_MATCH = 0, IDENTITY_MATCH };


    /**
     * This class wraps a kmer class
     * for efficient umi handling
     */
    class AlvKmer {
    private:
      Mer umi_;

    public:
      AlvKmer() = default;
      AlvKmer(AlvKmer&& other) = default;
      AlvKmer(AlvKmer& other) = default;
      AlvKmer(const AlvKmer& other) = default;
      AlvKmer& operator=(AlvKmer& other) = default;

      AlvKmer(const size_t& umi_length){
        umi_.k(umi_length);
      }

      static inline void k(int kIn) { Mer::k(kIn); }
      static inline int k() { return Mer::k(); }

      inline bool fromStr(const std::string& s) {
        return umi_.fromCharsSafe(s);
        // NOTE : @Avi -- what happens below if s is of
        // the appropriate length, but contains an 'N'
        // character?
        /*
        auto k = Mer::k();
        if (s.length() < k) {
          return false;
        }
        for (size_t i = 0; i < k; ++i) {
          umi_.shift_right(s[i]);
        }
        return true;
        */
      }

      inline void fromNum(uint64_t w) {
        umi_.word__(0) = w;
      }

      inline const Mer& umi() const { return umi_; }

      inline uint64_t umiWord() const { return umi_.word(0); }

      inline KmerMatchType isEquivalent(const Mer& m) const {
        return m.word(0) == umiWord()
          ? KmerMatchType::IDENTITY_MATCH : KmerMatchType::NO_MATCH;
      }

      inline KmerMatchType isEquivalent(uint64_t m) const {
        return m == umiWord() ? KmerMatchType::IDENTITY_MATCH : KmerMatchType::NO_MATCH;
      }

      inline std::string to_str() const {
        //std::string s = umi_.to_str();
        std::string s = umi_.toStr();
        std::reverse(s.begin(), s.end());
        return s;
      }

      bool operator==(const AlvKmer& rhs) const {
        return this->umi_ == rhs.umi_;
      }
      bool operator!=(const AlvKmer& rhs) const {
        return !this->operator==(rhs);
      }
    };
  }
}

#endif // __ALEVIN_KMER_HPP__
