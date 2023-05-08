#include <string>
#include <algorithm> // std::sort
#include <iostream>

#include "cedar/Taxa.h"

void TaxaNode::addInterval(uint64_t begin, uint64_t len, ReadEnd readEnd) {
    if (readEnd == ReadEnd::LEFT)
        lintervals.emplace_back(begin, begin+len);
    else
        rintervals.emplace_back(begin, begin+len);
}

void TaxaNode::updateScore() {
    score = 0;
    for (auto& it : lintervals) {
        score += it.end - it.begin;
    }
    for (auto& it : rintervals) {
        score += it.end - it.begin;
    }
    //std::cout << lintervals.size() << " " << rintervals.size() << " score: " << score << "\n";
}
void TaxaNode::updateIntervals(TaxaNode* child, ReadEnd readEnd) {

    std::vector<Interval>* intervals;
    
    // merge two sorted interval lists into parent
    // update parent score
    auto& childIntervals = child->getIntervals(readEnd);
    if (readEnd == ReadEnd::LEFT)
        intervals = &lintervals;
    else
        intervals = &rintervals;
    std::vector<Interval> parentIntervals(intervals->size());
    std::copy(intervals->begin(), intervals->end(), parentIntervals.begin());
    intervals->clear();
    intervals->reserve(parentIntervals.size()+childIntervals.size());
    
    std::vector<Interval>::iterator pit = parentIntervals.begin();
    std::vector<Interval>::iterator cit = childIntervals.begin();
    std::vector<Interval>::iterator fit = intervals->begin();

    // add the smallest interval as the first interval
    if (cit != childIntervals.end() && (pit == parentIntervals.end() || cit->begin < pit->begin)) {
        intervals->emplace_back(cit->begin, cit->end);
        cit++;
    }
    else if (pit != parentIntervals.end()) {
        intervals->emplace_back(pit->begin, pit->end);
        pit++;
    }
    else {
        std::cerr << "ERROR!! Both parent an child intervals were empty.\n";
        std::exit(1);
    }
    std::vector<Interval>::iterator cur;
    while (pit != parentIntervals.end() || cit != childIntervals.end()) {
        // find the smallest interval between the heads of the two lists
        if (cit != childIntervals.end() && (pit == parentIntervals.end() || cit->begin < pit->begin) )  {
            cur = cit;
            cit++;
        }
        else if (pit != parentIntervals.end()) {
            cur = pit;
            pit++;  
        }
        else {
            std::cerr << "ERROR!! Shouldn't even enter the loop.\n";
            std::exit(1);
        }
        // merge the new interval
        // Note: since both lists are sorted
        // the new interval's begin is always >= the last inserted interval's
        if (fit->end >= cur->begin) { // if the new interval has an overlap with the last inserted one
            fit->end = std::max(cur->end, fit->end); // merge them
        } else { // insert the interval as a separate one and move fit forward
            intervals->emplace_back(cur->begin, cur->end);
            fit++;
        }
    }
}
/**
 * Sorts intervals
 * Merge intervals if possible
 * Calculates score
**/
void TaxaNode::cleanIntervals(ReadEnd readEnd) {
    // if we were writing intervals wrt read position we wouldn't need this part
    std::vector<Interval>* intervals;
    if (readEnd == ReadEnd::LEFT)
        intervals = &lintervals;
    else
        intervals = &rintervals;
    
    std::sort(intervals->begin(), intervals->end(), 
    [](Interval& i1, Interval& i2){
        return i1.begin != i2.begin?i1.begin < i2.begin:i1.end < i2.end;
    });
    // merge intervals if necessary
    // calculate score / coverage !! this whole process is repetition of mapping coverage calc!!
    for (auto it=intervals->begin(); it != intervals->end();) {
        // start from next item and merge (and erase) as much as possible
        bool merged = true;
        auto next = it+1;
        while (next != intervals->end() && merged) {
            // if they overlap, merge them
            if (it->end >= next->begin) {
                if (it->end < next->end)
                    it->end = next->end;
                intervals->erase(next); // erase next after merging with it
                next = it+1; // it pointer and its next are always valid
            }
            else {// no overlap anymore and no merging. Update score and leave the loop
                merged = false;
            }
        }
        // there is nothing to merge with the current it, so increase it
        it++;
    }

    isCleaned = true;
}

bool TaxaNode::addChild(TaxaNode* child) { 
    // if this child was not already added
    if (activeChildren.find(child->getId()) == activeChildren.end()) {
        activeChildren.insert(child->getId()); // add it
        notIncorporatedChildrenCounter++; // increase not incorprated children counter
        return true;
    }
    return false;
}

void TaxaNode::reset() {
     lintervals.clear();
     rintervals.clear();
     activeChildren.clear();
     notIncorporatedChildrenCounter = 0;
     score = 0;
}

bool TaxaNode::compareIntervals(TaxaNode& other) {
    if (!isCleaned) {
        cleanIntervals(ReadEnd::LEFT);
        cleanIntervals(ReadEnd::RIGHT);
    }
    if (!other.isCleaned) {
        other.cleanIntervals(ReadEnd::LEFT);
        other.cleanIntervals(ReadEnd::RIGHT);
    }

    return lintervals == other.lintervals and rintervals == other.rintervals;
}

