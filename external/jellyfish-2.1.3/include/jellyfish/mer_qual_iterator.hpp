/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef __MER_QUAL_ITERATOR_HPP__
#define __MER_QUAL_ITERATOR_HPP__

#include <iterator>
#include <jellyfish/mer_dna.hpp>

namespace jellyfish {
template<typename SequencePool, typename MerType>
class mer_qual_iterator : public std::iterator<std::input_iterator_tag, MerType> {
  typename SequencePool::job* job_;
  std::string::const_iterator cseq_, eseq_;
  std::string::const_iterator cqual_, equal_;
  MerType                     m_; // mer
  MerType                     rcm_; // reverse complement mer
  unsigned int                filled_;
  const char                  min_qual_;
  const bool                  canonical_;
  size_t                      index_;

public:
  typedef MerType      mer_type;
  typedef SequencePool sequence_parser_type;

  mer_qual_iterator(SequencePool& seq, char min_qual, bool canonical = false) :
    job_(new typename SequencePool::job(seq)),
    filled_(0), min_qual_(min_qual), canonical_(canonical), index_(0)
  {
    if(job_->is_empty()) {
      delete job_;
      job_ = 0;
    } else {
      init_from_job();
      this->operator++();
    }
  }
  mer_qual_iterator() : job_(0), filled_(0), canonical_(false), index_(0) { }
  //  mer_iterator(const mer_iterator& rhs) : job_(rhs.job_), cseq_(rhs.cseq_), m_(rhs.m_), filled_(rhs.filled_) { }
  ~mer_qual_iterator() {
    delete job_;
  }

  bool operator==(const mer_qual_iterator& rhs) const { return job_ == rhs.job_; }
  bool operator!=(const mer_qual_iterator& rhs) const { return job_ != rhs.job_; }

  operator void*() const { return (void*)job_; }
  const mer_type& operator*() const { return !canonical_ || m_ < rcm_ ? m_ : rcm_; }
  const mer_type* operator->() const { return &this->operator*(); }
  mer_qual_iterator& operator++() {
    while(true) {
      while(cseq_ == eseq_) {
        ++index_;
        while(index_ >= (*job_)->nb_filled) {
          index_ = 0;
          job_->next();
          if(job_->is_empty()) {
            delete job_;
            job_  = 0;
            return *this;
          }
        }
        init_from_job();
        filled_ = 0;
        }

      do {
        const int  code = m_.code(*cseq_++);
        const char qual = cqual_ < equal_ ? *cqual_++ : std::numeric_limits<char>::min();
        if(code >= 0 && qual >= min_qual_) {
          m_.shift_left(code);
          if(canonical_)
            rcm_.shift_right(rcm_.complement(code));
          filled_ = std::min(filled_ + 1, mer_dna::k());
        } else
          filled_ = 0;
      } while(filled_ < m_.k() && cseq_ < eseq_);
      if(filled_ >= m_.k())
        break;
    }
    return *this;
  }

  mer_qual_iterator operator++(int) {
    mer_qual_iterator res(*this);
    ++*this;
    return res;
  }

private:
  void init_from_job() {
    std::string& seq = (*job_)->data[index_].seq;
    cseq_  = seq.begin();
    eseq_  = seq.end();
    std::string& quals = (*job_)->data[index_].qual;
    cqual_ = quals.begin();
    equal_ = quals.end();
  }
};

} // namespace jellyfish {

#endif /* __MER_QUAL_ITERATOR_HPP__ */
