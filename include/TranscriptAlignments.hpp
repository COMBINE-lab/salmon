#ifndef __TRANSCRIPT_ALIGNMENTS_HPP__
#define __TRANSCRIPT_ALIGNMENTS_HPP__

#include <vector>
#include "SalmonMath.hpp"

class SMEMAlignment;

class TranscriptAlignments {
    public:
        TranscriptAlignments() : alignments(std::vector<SMEMAlignment*>()),
                                 totalProb(salmon::math::LOG_0) {}

        std::vector<SMEMAlignment*> alignments;
        double totalProb;
        double logMassPrior;
        double logMassPosterior;
};

#endif //__TRANSCRIPT_ALIGNMENTS_HPP__
