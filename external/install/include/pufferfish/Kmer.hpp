#ifndef __COMBINELIB_KMERS_HPP__
#define __COMBINELIB_KMERS_HPP__

#include <cassert>
#include <climits>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <type_traits>

namespace combinelib {
namespace kmers {

#ifndef __DEFINE_LIKELY_MACRO__
#define __DEFINE_LIKELY_MACRO__
#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x), 1)
#define UNLIKELY(x) __builtin_expect((x), 0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif
#endif

/**
 *
 * The following lookup tables and reverse complement code is taken from
 *Jellyfish
 * https://github.com/gmarcais/Jellyfish/blob/master/include/jellyfish/mer_dna.hpp
 *
 **/

#define R -1
#define I -2
#define O -3
#define A 0
#define C 1
#define G 2
#define T 3
static constexpr int codes[256] = {
    O, O, O, O, O, O, O, O, O, O, I, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, R, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, A, R, C, R, O, O, G,
    R, O, O, R, O, R, R, O, O, O, R, R, T, O, R, R, R, R, O, O, O, O, O, O,
    O, A, R, C, R, O, O, G, R, O, O, R, O, R, R, O, O, O, R, R, T, O, R, R,
    R, R, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O};

static constexpr char complements[256] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N'};

#undef R
#undef I
#undef O
#undef A
#undef C
#undef G
#undef T

// Checkered mask. cmask<uint16_t, 1> is every other bit on
// (0x55). cmask<uint16_t,2> is two bits one, two bits off (0x33). Etc.
template <typename U, int len, int l = sizeof(U) * 8 / (2 * len)> struct cmask {
  static const U v =
      (cmask<U, len, l - 1>::v << (2 * len)) | ((static_cast<U>(1) << len) - 1);
};
template <typename U, int len> struct cmask<U, len, 0> {
  static const U v = 0;
};

// Fast reverse complement of one word through bit tweedling.
static inline uint64_t word_reverse_complement(uint64_t w, uint16_t k_) {
  typedef uint64_t U;
  w = ((w >> 2) & cmask<U, 2>::v) | ((w & cmask<U, 2>::v) << 2);
  w = ((w >> 4) & cmask<U, 4>::v) | ((w & cmask<U, 4>::v) << 4);
  w = ((w >> 8) & cmask<U, 8>::v) | ((w & cmask<U, 8>::v) << 8);
  w = ((w >> 16) & cmask<U, 16>::v) | ((w & cmask<U, 16>::v) << 16);
  w = (w >> 32) | (w << 32);
  return ((static_cast<U>(-1)) - w) >> (2 * (32 - k_));
}
static constexpr char revCodes[4] = {'A', 'C', 'G', 'T'};
/**
 * The above from Jellyfish (mer_dna.hpp)
 */

static constexpr int8_t rc_table[128] = {
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, // 15
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, // 31
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, // 787
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, // 63
    78, 84, 78, 71, 78, 78, 78, 67, 78, 78, 78, 78, 78, 78, 78, 78, // 79
    78, 78, 78, 78, 65, 65, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, // 95
    78, 84, 78, 71, 78, 78, 78, 67, 78, 78, 78, 78, 78, 78, 78, 78, // 101
    78, 78, 78, 78, 65, 65, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78  // 127
};

/**
 * Since we define these implementations at file scope in a header, we mark them
 *constant to
 * avoid duplicate symbol errors due to external linkage.
 **/

static decltype(codes[0]) codeForChar(char c) {
  return codes[static_cast<uint8_t>(c)];
}
static char charForCode(int i) { return revCodes[i]; }
static decltype(complements[0]) complement(char c) {
  return complements[static_cast<uint8_t>(c)];
}
static int complement(int i) { return 0x3 - i; }
static bool isValidNuc(int i) { return i >= 0; }
static bool isValidNuc(char c) { return isValidNuc(codeForChar(c)); }
static bool notValidNuc(int i) { return !isValidNuc(i); }
static bool notValidNuc(char c) { return !isValidNuc(c); }

// from :
// https://stackoverflow.com/questions/1392059/algorithm-to-generate-bit-mask
template <typename R> static constexpr R bitmask(unsigned int const onecount) {
  return (onecount == 0)
             ? 0
             : (static_cast<R>(-(onecount != 0)) &
                (static_cast<R>(-1) >> ((sizeof(R) * CHAR_BIT) - onecount)));
}

// table that contains bit patterns to mask out the top bits of a word.
// The table is such that maskTable[k] will mask out the top (64 - 2*k) bits of
// the word.
static const constexpr uint64_t maskTable[] = {
    bitmask<uint64_t>(0),  bitmask<uint64_t>(2),  bitmask<uint64_t>(4),
    bitmask<uint64_t>(6),  bitmask<uint64_t>(8),  bitmask<uint64_t>(10),
    bitmask<uint64_t>(12), bitmask<uint64_t>(14), bitmask<uint64_t>(16),
    bitmask<uint64_t>(18), bitmask<uint64_t>(20), bitmask<uint64_t>(22),
    bitmask<uint64_t>(24), bitmask<uint64_t>(26), bitmask<uint64_t>(28),
    bitmask<uint64_t>(30), bitmask<uint64_t>(32), bitmask<uint64_t>(34),
    bitmask<uint64_t>(36), bitmask<uint64_t>(38), bitmask<uint64_t>(40),
    bitmask<uint64_t>(42), bitmask<uint64_t>(44), bitmask<uint64_t>(46),
    bitmask<uint64_t>(48), bitmask<uint64_t>(50), bitmask<uint64_t>(52),
    bitmask<uint64_t>(54), bitmask<uint64_t>(56), bitmask<uint64_t>(58),
    bitmask<uint64_t>(60), bitmask<uint64_t>(62)};

constexpr const uint64_t nucleotidesPerByte = 4;

// from :
// https://stackoverflow.com/questions/31952237/looking-for-a-constexpr-ceil-function
constexpr uint64_t ceil(double num) {
  return (static_cast<double>(static_cast<uint64_t>(num)) == num)
             ? static_cast<uint64_t>(num)
             : static_cast<uint64_t>(num) + ((num > 0) ? 1 : 0);
}

constexpr uint64_t numWordsRequired(uint64_t K) {
  return ceil(K / (1.0 * nucleotidesPerByte * (sizeof(uint64_t))));
}

/**
 * @returns the binary encoding for character c
 **/
static int64_t doEncodeBinary(char c) { return codes[static_cast<uint8_t>(c)]; }

/**
 * @returns true of the character `c` was a valid nucleotide and false
 *otherwise. The corresponding
 * code for this character is placed in the parameter `code`.
 **/
static bool encodeBinary(char c, int64_t& code) {
  code = codes[static_cast<uint8_t>(c)];
  return code >= 0;
}

static char decodeBinary(uint64_t n) { return revCodes[n]; }

/**
 * Convert an ascii character to the corresponding 2-bit encoding
 *
 * Following the encoding suggested [here](https://www.biostars.org/p/113640/),
 *originally
 * suggested by G. Rizk:
 * A : 0
 * C : 1
 * G : 3
 * T : 2
 * N : 4
 *
 * This function will work with both lower and upper case nucleotides.
 * @ASSUMPTION : c is in {A,C,G,T,N,a,c,g,t,n}
 **/
static uint64_t charToBitsGATB(char c) {
  // Convert to uppercase
  // https://stackoverflow.com/questions/10688831/fastest-way-to-capitalize-words
  return static_cast<uint64_t>((((c & ~0x20) >> 1) & 0x03) + ((c & 0x08) >> 3));
}

// Adapted from
// https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/8c9933a1685e0ab50c7d8b7926c9068bc0c9d7d2/src/main.c#L36
static void reverseComplement(const std::string& seq, std::string& readWork) {
  readWork.resize(seq.length(), 'A');
  int32_t end = seq.length() - 1, start = 0;
  while (LIKELY(start < end)) {
    readWork[start] = (char)rc_table[(int8_t)seq[end]];
    readWork[end] = (char)rc_table[(int8_t)seq[start]];
    ++start;
    --end;
  }
  // If odd # of bases, we still have to complement the middle
  if (start == end) {
    readWork[start] = (char)rc_table[(int8_t)seq[start]];
  }
}

static std::string reverseComplement(const std::string& seq) {
  std::string work;
  reverseComplement(seq, work);
  return work;
}

static std::string stringRevComp(const std::string& seq) {
  return reverseComplement(seq);
}

/**
* From https://www.biostars.org/p/113640/.  This only works for a given word
*right now;
* will determine how to best generalize later.
**/
static uint64_t word_reverse_complement_gatb(uint64_t x, size_t k_) {
  uint64_t res = x;
  res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
  res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
  res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
  res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
  res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
  res = res ^ 0xAAAAAAAAAAAAAAAA;
  return (res >> (2 * (32 - k_)));
}

// Some template magic to detect if a template type has a ``length()'' function.
template <typename...> using combinelib_void_t = void;

template <typename, typename = void>
struct has_length : public std::false_type {};

template <typename T>
struct has_length<T, combinelib_void_t<decltype(T().length())>>
    : public std::true_type {};

/**
 * The first template parameter, K, is the maximum length (in nucleotides)
 * of the k-mer that can be represented with this class.
 *
 * The second template parameter, CID, is a class-type specific tag that
 * will allow all instances of this particular class to share a value of
 * their k.  This idea is used in Jellyfish, which inspired the use here.
 **/
template <uint64_t K, uint64_t CID = 0> class Kmer {
  static_assert(
      K <= 32,
      "Currently, the Kmer class can only represent k-mers of size <= 32");

public:
  using base_type = uint64_t;

  explicit Kmer() {}

  template <
      typename ViewT,
      typename = typename std::enable_if<has_length<ViewT>::value, void>::type>
  Kmer(ViewT& v) {
    fromChars(v);
  }

  template <
      typename IterT,
      typename = typename std::enable_if<!has_length<IterT>::value, void>::type>
  Kmer(IterT v) {
    fromChars(v);
  }

  Kmer(const Kmer& other) = default;
  Kmer(Kmer&& other) = default;
  Kmer(Kmer& other) = default;
  Kmer& operator=(Kmer& other) = default;
  Kmer& operator=(const Kmer& other) = default;

  template <
      typename IterT,
       typename = typename std::enable_if<!has_length<IterT>::value, void>::type>
  Kmer& operator=(IterT iter) {
    fromChars(iter);
    return *this;
  }

  template <
      typename ViewT,
      typename = typename std::enable_if<has_length<ViewT>::value, void>::type>
  Kmer& operator=(ViewT& v) {
    fromChars(v);
    return *this;
  }

  /**
   * FROM JELLYFISH
   **/
  // Get bits [start, start+len). start must be < 2k, len <=
  // sizeof(base_type) and start+len < 2k. No checks
  // performed. start and len are in bits, not bases.
  base_type get_bits(unsigned int start, unsigned int len) const {
    unsigned int q = start / wbits;
    unsigned int r = start % wbits;

    base_type res = data_[q] >> r;
    /* this case should never happen with k <= 32
     * and valid input arguments.
    if(len > wbits - r)
      res |= data_[q + 1] << (wbits - r);
    */
    return len < (unsigned int)wbits ? res & (((base_type)1 << len) - 1) : res;
  }

  /**
   * Populate this kmer by consuming characters pointed to by iter.
   *
   * @ASSUMPTIONS:
   *  There are at least k_ characters to consume or
   *  (2) we will encounter a non-nucleotide (i.e. \0) character
   **/
  template <
      typename IterT,
      typename = typename std::enable_if<!has_length<IterT>::value, void>::type>
  bool fromChars(IterT iter) {
    // std::memset(&data_[0], 0, sizeof(data_));
    data_[0] = 0;
    auto toConsume = 1; // numWordsRequired(k_);
    int64_t code{0};
    bool success = true;
    int32_t remK = static_cast<int32_t>(k_);
    for (int32_t w = 0; w < toConsume; ++w) {
      int32_t shift = std::min((2 * remK) - 2, 62);
      auto& currWord = data_[w];
      for (; remK > 0 and shift >= 0; ++iter, --remK, shift -= 2) {
        // success &= encodeBinary(*iter, code);
        if (!encodeBinary(*iter, code))
          return false;
        currWord |= (code << shift);
      }
    }
    return success;
  }

  /**
   *  This is the same as the above function, but it will be called if the
   *argument type
   *  `ViewT` has a "length()" member.  In that case, the function will
   *additionally check
   *  that the length of `v` is >= k_.
   **/
  template <
      typename ViewT,
      typename = typename std::enable_if<has_length<ViewT>::value, void>::type>
  bool fromCharsSafe(ViewT& v) {
    return (v.length() >= k_) ? fromChars(v.begin()) : false;
  }

  /**
   *  This is a convenience function taht lets us call fromChars on a string, or
   *string_vew (or similar object);
   **/
  template <
      typename ViewT,
      typename = typename std::enable_if<has_length<ViewT>::value, void>::type>
  bool fromChars(ViewT& v) {
    return fromChars(v.begin());
  }
  bool fromChars(Kmer& k) {
      data_[0] = k.data_[0];
      return true;
  }
  /**
   * Append the character `c` to the end of the k-mer
   **/
  uint64_t append(char c) {
    auto r = (data_[0] >> (2 * k_ - 2)) & 0x03;
    data_[0] = maskTable[k_] & ((data_[0] << 2) | doEncodeBinary(c));
    return r;
  }

  /**
   * Prepend the character `c` to the beginning of the k-mer
   **/
  uint64_t prepend(char c) {
    auto r = (data_[0] & 0x03);
    data_[0] = (data_[0] >> 2) | (doEncodeBinary(c) << (2 * k_ - 2));
    return r;
  }

   /**
   * Append the character `c` to the end of the k-mer
   **/
  uint64_t append(int i) {
    auto r = (data_[0] >> (2 * k_ - 2)) & 0x03;
    data_[0] = maskTable[k_] & ((data_[0] << 2) | static_cast<base_type>(i));
    return r;
  }

  /**
   * Prepend the character `c` to the beginning of the k-mer
   **/
  uint64_t prepend(int i) {
    auto r = (data_[0] & 0x03);
    data_[0] = (data_[0] >> 2) | (static_cast<base_type>(i) << (2 * k_ - 2));
    return r;
  }

  /**
   * @returns a `uint64_t` that represents the encoded `idx`-th word of this
   *k-mer
   **/
  uint64_t word(uint32_t idx) const { return data_[idx]; }

  /**
   * @returns a reference to the `uint64_t` that represents the encoded `idx`-th
   *word of this k-mer
   **/
  uint64_t& word__(uint32_t idx) { return data_[idx]; }

  const base_type* data() const { return &data_[0]; }

  /**
   *  @returns the number of bytes required by this k-mer
   **/
  uint64_t sizeInBytes() const { return sizeof(data_); }

  /**
   *  @returns the number of words required by this k-mer
   **/
  uint64_t sizeInWords() const { return sizeof(data_) / sizeof(base_type); }

  /**
   *  @returns the number of words required by this k-mer
   **/
  uint64_t nb_words() const { return sizeInWords(); }

  /**
   * Set the dynamic length of this k-mer class to be kIn nucleotides.
   * @returns the value of k for this class prior to this update.
   **/
  static uint16_t k(uint16_t kIn) {
    assert(kIn < K);
    std::swap(k_, kIn);
    return kIn;
  }

  /**
   * @returns the value of k used for this k-mer class
   **/
  static uint16_t k() { return k_; }

  std::string toStr() const {
    std::string s(k_, 'X');
    auto& d = data_[0];
    int32_t offset = (2 * k_) - 2;
    for (int32_t idx = 0; offset >= 0; offset -= 2, ++idx) {
      s[idx] = decodeBinary((d >> offset & 0x03));
    }
    return s;
  }

  bool isHomoPolymer() const {
    auto nuc = data_[0] & 0x3;
    return (data_[0] == (maskTable[k_] & ((data_[0] << 2) | nuc)));
  }
  bool is_homopolymer() const { return isHomoPolymer(); }

  void rc() { data_[0] = word_reverse_complement(data_[0], k_); }

  Kmer<K, CID> getRC() const {
    Kmer<K, CID> nk;
    nk.data_[0] = word_reverse_complement(data_[0], k_);
    return nk;
  }

  void canonicalize() {
    auto wrc = word_reverse_complement(data_[0], k_);
    data_[0] = (wrc < data_[0]) ? wrc : data_[0];
  }

  Kmer<K, CID> getCanonical() {
    Kmer<K, CID> rck = getRC();
    return (rck < *this) ? rck : *this;
  }

  template <uint64_t KP, uint64_t CIDP>
  friend std::ostream& operator<<(std::ostream& os, const Kmer<KP, CIDP>& k);

  template <uint64_t KP, uint64_t CIDP>
  friend bool operator==(const Kmer<KP, CIDP>& lhs, const Kmer<KP, CIDP>& rhs);

  template <uint64_t KP, uint64_t CIDP>
  friend bool operator!=(const Kmer<KP, CIDP>& lhs, const Kmer<KP, CIDP>& rhs);

  template <uint64_t KP, uint64_t CIDP>
  friend bool operator<(const Kmer<KP, CIDP>& lhs, const Kmer<KP, CIDP>& rhs);

  template <uint64_t KP, uint64_t CIDP>
  friend bool operator>(const Kmer<KP, CIDP>& lhs, const Kmer<KP, CIDP>& rhs);

private:
  base_type data_[numWordsRequired(K)] = {};
  static uint16_t k_;
  static constexpr const int32_t wshift = sizeof(base_type) * 8 - 2; // left shift in 1 word
  static constexpr const int32_t wbases = 4 * sizeof(base_type); // bases in a word
  static constexpr const int32_t wbits  = 8 * sizeof(base_type); // bits in a word
};

template <uint64_t K, uint64_t CID> uint16_t Kmer<K, CID>::k_ = 0;

template <uint64_t K, uint64_t CID>
std::ostream& operator<<(std::ostream& os, const Kmer<K, CID>& k) {
  os << k.toStr();
  return os;
}

template <uint64_t K, uint64_t CID>
bool operator==(const Kmer<K, CID>& lhs, const Kmer<K, CID>& rhs) {
  return lhs.data_[0] == rhs.data_[0];
}

template <uint64_t K, uint64_t CID>
bool operator!=(const Kmer<K, CID>& lhs, const Kmer<K, CID>& rhs) {
  return !(lhs == rhs);
}

template <uint64_t K, uint64_t CID>
bool operator<(const Kmer<K, CID>& lhs, const Kmer<K, CID>& rhs) {
  return (lhs.data_[0] < rhs.data_[0]);
}

template <uint64_t K, uint64_t CID>
bool operator>(const Kmer<K, CID>& lhs, const Kmer<K, CID>& rhs) {
  return (lhs.data_[0] > rhs.data_[0]);
}

} // namespace kmers
} // namespace combinelib

#endif // __COMBINELIB_KMERS_HPP__ 