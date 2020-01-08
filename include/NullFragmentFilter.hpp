#ifndef __NULL_FRAGMENT_FILTER_HPP__
#define __NULL_FRAGMENT_FILTER_HPP__

template <typename FragT> class NullFragmentFilter {
public:
  inline void processFrag(FragT* /*f*/) {}
};

#endif // __NULL_FRAGMENT_FILTER_HPP__
