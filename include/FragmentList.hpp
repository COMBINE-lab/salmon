#ifndef FRAGMENT_LIST_HPP
#define FRAGMENT_LIST_HPP

extern "C" {
#include "container.h"
}

class FragmentList {
public:
  FragmentList();

  // Free all of the structures allocated for this fragment list
  ~FragmentList();

  void freeFragList(Container* frags);

  bool computeBestChain_(Container* frags, double& maxScore, uint32_t& bestPos);

  void computeBestChain();

  void addFragMatch(uint32_t refStart, uint32_t queryStart, uint32_t len);

  void addFragMatch(uint32_t refStart, uint32_t refEnd, uint32_t queryStart,
                    uint32_t queryEnd);

  void addFragMatchRC(uint32_t refStart, uint32_t queryStart, uint32_t len);

  void addFragMatchRC(uint32_t refStart, uint32_t refEnd, uint32_t queryStart,
                      uint32_t queryEnd);

  bool isForward();

  Container* fragments;
  Container* fragmentsRC;
  size_t numFrag;
  double bestHitScore;
  uint32_t bestHitPos;
  bool isForward_{true};
};

#endif // FRAGMENT_LIST_HPP
