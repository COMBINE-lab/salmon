// This code is derived and modified from Pall Melsted's
// KmerIterator class (part of bfgraph) :
// https://github.com/pmelsted/bfgraph/blob/master/src/KmerIterator.cpp

#ifndef MER_ITERATOR_HPP
#define MER_ITERATOR_HPP

#include "CanonicalKmer.hpp"
#include "string_view.hpp"
#include <iterator>

namespace pufferfish {
  namespace kmers = combinelib::kmers;
// class CanonicalKmerIterator : public std::iterator<std::input_iterator_tag,
// std::pair<CanonicalKmer, int>, int> {
class CanonicalKmerIterator
  : public std::iterator<std::input_iterator_tag, std::pair<CanonicalKmer, int>, int> {
  stx::string_view s_;
  std::pair<CanonicalKmer, int> p_;
  //CanonicalKmer km_;
  //int pos_;
  bool invalid_;
  int lastinvalid_;
  int k_;

public:
  typedef std::pair<CanonicalKmer,int> value_type;
  typedef value_type& reference;
  typedef value_type* pointer;
  typedef std::input_iterator_tag iterator_category;
  typedef int64_t difference_type;
  CanonicalKmerIterator()
    : s_(), p_(), /*km_(), pos_(),*/ invalid_(true), lastinvalid_(-1),
        k_(CanonicalKmer::k()) {}
  CanonicalKmerIterator(const std::string& s)
    : s_(s), p_(), /*km_(), pos_(),*/ invalid_(false), lastinvalid_(-1),
        k_(CanonicalKmer::k()) {
    find_next(-1, -1);
  }
  CanonicalKmerIterator(const CanonicalKmerIterator& o)
    : s_(o.s_), p_(o.p_), /*km_(o.km_), pos_(o.pos_),*/ invalid_(o.invalid_),
        lastinvalid_(o.lastinvalid_), k_(o.k_) {}

  CanonicalKmerIterator& operator=(const CanonicalKmerIterator& o) {
    if (this != &o) {
      s_ = o.s_;
      p_ = o.p_;
      invalid_ = o.invalid_;
      lastinvalid_ = o.lastinvalid_;
      k_ = o.k_;
    }
    return *this;
  }

private:
  inline void find_next(int i, int j) {
    ++i;
    ++j;
    // j is the last nucleotide in the k-mer we're building
    for (; j < static_cast<int>(s_.length()); ++j) {
      // get the code for the last nucleotide, save it as c
      int c = kmers::codeForChar(s_[j]);
      // c is a valid code if != -1
      if (c != -1) {
        p_.first.shiftFw(c);
        if (j - lastinvalid_ >= k_) {
          p_.second = i;
          return;
        }
      } else {
        // this k-mer is clearly not valid
        // if c is not a valid code, then j is the last invalid position
        lastinvalid_ = j;
        i = j + 1;
      }
    }
    invalid_ = true;
  }

public:
  inline stx::string_view seq() { return s_; }
  // use:  ++iter;
  // pre:
  // post: *iter is now exhausted
  //       OR *iter is the next valid pair of kmer and location
  inline CanonicalKmerIterator& operator++() {
    auto lpos = p_.second + k_;
    invalid_ = invalid_ || lpos >= static_cast<int>(s_.length());
    if (!invalid_) {
      find_next(p_.second, lpos - 1);
      /** --- implementation that doesn't skip non-{ACGT}
      int c = kmers::codeForChar(s_[lpos]);
      if (c!=-1) { km_.shiftFw(c); } else { lastinvalid_ = pos_ + k_; }
      ++pos_;
      */
    }
    return *this;
  }

  // use:  iter++;
  // pre:
  // post: iter has been incremented by one
  inline CanonicalKmerIterator operator++(int) {
    CanonicalKmerIterator tmp(*this);
    operator++();
    return tmp;
  }

  // use:  iter += constant int;
  // pre:
  // post: *iter is now exhausted
  //       OR *iter is the next valid pair of kmer and location after advancing
  inline CanonicalKmerIterator& operator+=(int advance) {
    //CanonicalKmerIterator tmp(*this) ;
    while(advance > 0){
        operator++() ;
        advance-- ;
    }
    return *this;
  }


  // use:  val = (a == b);
  // pre:
  // post: (val == true) if a and b are both exhausted
  //       OR a and b are in the same location of the same string.
  //       (val == false) otherwise.
  inline bool operator==(const CanonicalKmerIterator& o) {
    return (invalid_ || o.invalid_) ? invalid_ && o.invalid_
                                    : ((p_.first == o.p_.first) && (p_.second == o.p_.second));
  }

  inline bool operator!=(const CanonicalKmerIterator& o) {
    return !this->operator==(o);
  }

  inline bool kmerIsValid() { return (p_.second + k_ - lastinvalid_ > k_); }

  // use:  p = *iter;
  // pre:
  // post: p is NULL or a pair of Kmer and int
  inline reference operator*() { return p_; }

  // use:  example 1: km = iter->first;
  //       example 2:  i = iter->second;
  // pre:  *iter is not NULL
  // post: km will be (*iter).first, i will be (*iter).second
  inline pointer operator->() { return &(operator*()); }

  void jumpTo(int pos) {
    lastinvalid_ = pos-1;
    find_next(pos-1,(pos-1));
  }


private:
  /*
// use:  find_next(i,j, last_valid);
// pre:
// post: *iter is either invalid or is a pair of:
//       1) the next valid kmer in the string that does not have any 'N'
//       2) the location of that kmer in the string
inline void find_next(int i, int j){//}, bool last_valid) {
++i;
++j;
bool valid{false};
// j is the last nucleotide in the k-mer we're building
while (j < s_.length()) {
  // get the code for the last nucleotide, save it as c
  int c = kmers::codeForChar(s_[j]);
  // c is a valid code if != -1
  if (c != -1) {
    km_.shiftFw(c);
    valid = (j - lastinvalid_ >= k_);
  } else {
    // if c is not a valid code, then j is the last invalid position
    lastinvalid_ = j;
    // the start position is the next (potentially) valid position
    i = j+1;
    // this k-mer is clearly not valid
    valid = false;
  }
  if (valid) {
    //p_.second = i;
    pos_ = i;
    return;
  }
  ++j;
}
invalid_ = true;
}
*/

private:
};
}

#endif // MER_ITERATOR_HPP
