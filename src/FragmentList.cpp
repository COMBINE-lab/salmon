#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>

extern "C" {
#include "container.h"
#include "debug.h"
#include "fileio.h"
#include "info.h"
#include "manopt.h"
#include "rangetree.h"
#include "slchain.h"
#include "sltypes.h"
//#include "clasp.h"
}

#include "FragmentList.hpp"

#define SOP ((unsigned char)(0 << 0))
#define LIN ((unsigned char)(1 << 0))

unsigned char mute;

FragmentList::FragmentList() {
  fragments = (Container*)malloc(sizeof(Container));
  bl_containerInit(fragments, 20, sizeof(slmatch_t));
  fragmentsRC = (Container*)malloc(sizeof(Container));
  bl_containerInit(fragmentsRC, 20, sizeof(slmatch_t));
}

// Free all of the structures allocated for this fragment list
FragmentList::~FragmentList() {
  freeFragList(fragments);
  freeFragList(fragmentsRC);
}

void FragmentList::freeFragList(Container* frags) {
  for (uint32_t i = 0; i < bl_containerSize(frags); i++) {
    slmatch_t* sl = (slmatch_t*)bl_containerGet(frags, i);
    if (sl->chain != nullptr) {
      std::cerr << "still at least one chain not freed before end: " << i
                << "\n";
      exit(-1);
      bl_slchainDestruct(sl->chain);
      free(sl->chain);
    }
  }
  bl_containerDestruct(frags, bl_slmatchDestruct);
  free(frags);
}

bool FragmentList::computeBestChain_(Container* frags, double& maxScore,
                                     uint32_t& bestPos) {
  double epsilon = 0.0;
  double lambda = 1.0;
  unsigned char chainMode = LIN;
  // Did we find a chain with a score higher than the original
  // maxScore parameter?
  bool updatedMaxScore = false;

  /* sort fragments */
  qsort(frags->contspace, bl_containerSize(frags), sizeof(slmatch_t),
        cmp_slmatch_qsort);
  uint32_t begin = 0;
  for (uint32_t i = 1; i <= bl_containerSize(frags); i++) {
    /*
     * end of fragments list or different database sequence
     * --> process fragment[begin]...fragment[i-1], write output
     *     and free chains (less memory consumption with large input files)
     */
    if (i == bl_containerSize(frags) ||
        ((slmatch_t*)bl_containerGet(frags, begin))->subject !=
            ((slmatch_t*)bl_containerGet(frags, i))->subject) {
      // fprintf(info.dev, "%d\t%d\n", begin, i-begin);
      if (chainMode == SOP) {
        /* only use chaining without clustering if no ids are specified */
        bl_slChainSop((slmatch_t*)frags->contspace + begin, i - begin, epsilon,
                      lambda, 10);
        /*
           bl_slClusterSop((slmatch_t *) info.fragments->contspace + begin, i -
           begin, info.epsilon, info.lambda, info.maxgap);
           */
      } else {
        bl_slChainLin((slmatch_t*)frags->contspace + begin, i - begin, epsilon,
                      lambda, 10);
        /*
           bl_slClusterLin((slmatch_t *) info.fragments->contspace + begin, i -
           begin, info.epsilon, info.lambda, info.maxgap);
           */
      }

      for (uint32_t j = begin; j < i; j++) {
        slmatch_t* match = (slmatch_t*)bl_containerGet(frags, j);

        if (match->chain) {
          slchain_t* chain = (slchain_t*)match->chain;
          if (chain->scr >= maxScore) {
            maxScore = chain->scr;
            bestPos = chain->p;
            updatedMaxScore = true;
          }

          /*
          // output matches (if desired)
          if (info.outputm){
          fprintf(info.dev, "M\t");
          if (!info.outputorig){
          if (info.idcol != NULL){
          fprintf(info.dev, "%s\t",
           *(char **) bl_containerGet(info.subject, match->subject));
           }
           fprintf(info.dev, "%d\t%d\t%d\t%d\t%.3f\n", match->i,
           match->i + match->j - 1, match->p,
           match->p + match->q - 1, match->scr);
           }
          // output in original format as input
          else {
          fprintf(info.dev, "%s\n", *(char **) bl_containerGet(info.lines, j));
          }
          }
          if (match->chain){
          slchain_t *chain = (slchain_t *) match->chain;
          if (info.outputc && chain->scr >= info.minscore &&
          bl_containerSize(chain->matches) >= info.minfrag){
          fprintf(info.dev, "C\t");
          if (info.idcol != NULL){
          fprintf(info.dev, "%s\t", *(char **) bl_containerGet(info.subject,
          chain->subject));
          }
          fprintf(info.dev, "%d\t%d\t%d\t%d\t%.3f\n", chain->i,
          chain->i + chain->j - 1, chain->p,
          chain->p + chain->q - 1, chain->scr);
          }
          // output chains and fragments (if requested)
          if (info.outputf && chain->scr >= info.minscore &&
          bl_containerSize(chain->matches) >= info.minfrag){
          for (k = 0; k < bl_containerSize(chain->matches); k++){
          slmatch_t *frag = *(slmatch_t **)
          bl_containerGet(chain->matches, k);
          fprintf(info.dev, "F\t");
          if (!info.outputorig){
          if (info.idcol != NULL){
          fprintf(info.dev, "%s\t",
           *(char **) bl_containerGet(info.subject, frag->subject));
           }
           fprintf(info.dev, "%d\t%d\t%d\t%d\t%.3f\n", frag->i,
           frag->i + frag->j - 1, frag->p, frag->p + frag->q - 1,
           frag->scr);
           }
          // output in original format as input
          else {
          fprintf(info.dev, "%s\n",
           *(char **) bl_containerGet(info.lines, frag->idx));
           }
           }
           }
           */
          bl_slchainDestruct(chain);
          free(chain);
          match->chain = nullptr;
        } /* END OF if (frag->chain) */
      }   /* END OF for (j = begin; j < i; j++) */
      begin = i;
    } /* END OF  if (i == bl_containerSize(info.fragments) ||
         ((slmatch_t *) bl_containerGet(info.fragments, begin))->subject !=
         ((slmatch_t *) bl_containerGet(info.fragments, i))->subject) */
  }   /* END OF for (i = 1; i <= bl_containerSize(info.fragments); i++) */

  return updatedMaxScore;
}

bool FragmentList::isForward() { return isForward_; }

void FragmentList::computeBestChain() {
  double maxScore = 0.0;
  uint32_t bestPos = 0;

  // we don't need the return value from the first call
  static_cast<void>(computeBestChain_(fragments, maxScore, bestPos));
  bool reverseIsBest = computeBestChain_(fragmentsRC, maxScore, bestPos);
  isForward_ = not reverseIsBest;

  bestHitScore = maxScore;
  bestHitPos = bestPos;
}

void FragmentList::addFragMatch(uint32_t refStart, uint32_t queryStart,
                                uint32_t len) {
  addFragMatch(refStart, refStart + len, queryStart, queryStart + len);
}

void FragmentList::addFragMatch(uint32_t refStart, uint32_t refEnd,
                                uint32_t queryStart, uint32_t queryEnd) {
  slmatch_t frag;
  bl_slmatchInit(&frag, 0);
  frag.p = refStart;
  frag.q = refEnd;
  frag.i = queryStart;
  frag.j = queryEnd;
  frag.scr = static_cast<double>(queryEnd) - queryStart;
  bl_containerAdd(fragments, &frag);
}

void FragmentList::addFragMatchRC(uint32_t refStart, uint32_t queryStart,
                                  uint32_t len) {
  addFragMatchRC(refStart, refStart + len, queryStart, queryStart + len);
}

void FragmentList::addFragMatchRC(uint32_t refStart, uint32_t refEnd,
                                  uint32_t queryStart, uint32_t queryEnd) {
  slmatch_t frag;
  bl_slmatchInit(&frag, 0);
  frag.p = refStart;
  frag.q = refEnd;
  frag.i = queryStart;
  frag.j = queryEnd;
  frag.scr = static_cast<double>(queryEnd) - queryStart;
  bl_containerAdd(fragmentsRC, &frag);
}

/*
int main(int argc, char* argv[]) {
    FragmentList fl;
    fl.addFragment(0, 0, 15);
    fl.addFragment(25, 25, 30);
    fl.addFragment(10, 10, 20);
    fl.computeBestChain();
    std::cerr << "maxChain score = " << fl.maxChainScore << "\n";
}
*/
