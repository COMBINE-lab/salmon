#include "RefSeqConstructor.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "PufferfishLossyIndex.hpp"

#include <sparsepp/spp.h>

//#define verbose true

#define suffixIfFw true
#define prefixIfFw false


/* story 
NOTE: Since we are ALWAYS traversing the mems --and hence
the unmapped sequences-- in the forward direction wrt the transcript,
whenever we are moving backward in a contig that means that the contig
is in reverse orientation wrt the transcript, so the sequence we fetch from the contig
should be reverse-complemented!!!
EXCEPTION ::::: constructing the string before the first unimem :(((
TODO: make sure about the fact above startp : relative position in curContig
that toBeAligned string starts endp : relative position in endContig that toBeAligned string ends
txpDist : number of bases that should be fetched before we stop spreading the path through BFS
seq : will be appended through out the traversal
*/


template <typename PufferfishIndexT>
RefSeqConstructor<PufferfishIndexT>::RefSeqConstructor(PufferfishIndexT* pfi,
                                                       spp::sparse_hash_map<uint32_t,
                                                       pufferfish::util::ContigBlock>* contigSeqCache)
                                                       : pfi_(pfi), contigSeqCache_(contigSeqCache) { k = pfi_->k(); }
template <typename PufferfishIndexT>
Task RefSeqConstructor<PufferfishIndexT>::fillSeq(size_t tid,
                                             size_t tpos,
                                             bool isCurContigFw,
                                             pufferfish::util::ContigBlock& curContig,
                                             uint32_t startp,
                                             pufferfish::util::ContigBlock& endContig,
                                             uint32_t endp,
                                             bool isEndContigFw,
                                             uint32_t txpDist,
                                                  std::string& seq,
                                                  bool verbose) {
  
  //std::cerr << curContig.isDummy()<<endContig.isDummy() << " ";

  if (curContig.contigIdx_ == endContig.contigIdx_ && endp-startp-1 == txpDist) {
    append(seq, curContig, startp, endp, isCurContigFw);
    //std::cerr << txpDist << "-" << seq.length() <<  " ";
    return Task::SUCCESS;
  }
  size_t endTpos = tpos + txpDist + 1;
  
  if (!endContig.isDummy()) {
    auto endRemLen = remainingLen(endContig, endp, isEndContigFw, prefixIfFw);
    if (endRemLen >= txpDist) {
      appendByLen(seq, endContig, endp, txpDist, isEndContigFw, prefixIfFw);
      return Task::SUCCESS;
    } else if (endRemLen > 0) {
      appendByLen(seq, endContig, endp, endRemLen, isEndContigFw, prefixIfFw);
      txpDist-=endRemLen;
      if (isCurContigFw) endp -= endRemLen; // as it is prefix, in case of forward contig, we update the endp by reducing it
      else endp += endRemLen;
      tpos += endRemLen;
    }
  }
  if (curContig.isDummy()) {
    std::string tmp;
    doDFS(tid, endTpos, isEndContigFw, endContig, endp, curContig, isCurContigFw, txpDist, false, tmp);
    seq = rc(tmp) + seq;
    return Task::SUCCESS;
  }
  if (verbose) std::cerr << std::this_thread::get_id() << " "  << "\n\nWOOOOOT!! Got to bfs\n";
  return doDFS(tid, tpos, isCurContigFw, curContig, startp, endContig, isEndContigFw, txpDist, true, seq, verbose);
}

template <typename PufferfishIndexT>
Task RefSeqConstructor<PufferfishIndexT>::doDFS(size_t tid,
                                                size_t tpos,
                                                bool isCurContigFw,
                                                pufferfish::util::ContigBlock& curContig,
                                                uint32_t startp,
                                                pufferfish::util::ContigBlock& endContig,
                                                bool isEndContigFw,
                                                uint32_t txpDist,
                                                bool walkForward,
                                                std::string& seq,
                                                bool verbose) {



  if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\n[doBFS]\n" << "called for txp " << tid << " with pos " << tpos << " with curr contig: "
                        << curContig.contigIdx_ << " of length " << curContig.contigLen_ << " ori " << isCurContigFw
                        << " start " << startp << "\n"
                        << "end contig index "<< endContig.contigIdx_ << " is dummy: " << endContig.isDummy() << "\ntxpDist: " << txpDist << "\n";

  if (startp >= curContig.contigLen_) {
      std::cerr << std::this_thread::get_id() << " "  << "ERROR!!! shouldn't happen ---> startp >= curContig.contigLen_ : " << startp << ">" << curContig.contigLen_ << "\n";
      std::cerr << std::this_thread::get_id() << " "  << "called for txp " << tid << " with pos " << tpos << " with curr contig: "
                << curContig.contigIdx_ << " of length " << curContig.contigLen_ << " ori " << isCurContigFw
                << " start " << startp << "\n"
                << "end contig index "<< endContig.contigIdx_
                << "\nis end dummy: " << endContig.isDummy()
                << " is start dummy: " << curContig.isDummy() << "\ntxpDist: " << txpDist << "\n";
      std::exit(1);
  }
    // used in all the following terminal conditions
    auto remLen = remainingLen(curContig, startp, walkForward?isCurContigFw:!isCurContigFw, suffixIfFw);
    if (remLen >= txpDist) {
      if (endContig.isDummy()){
        appendByLen(seq, curContig, startp, txpDist, walkForward?isCurContigFw:!isCurContigFw, suffixIfFw);
        if (verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[After append] " << seq << "\n";
        return Task::SUCCESS;
      }
      // DON'T GET STUCK IN INFINITE LOOPS
      // if we have not reached the last contigId
      // and also end of the path is NOT open
      // then if the remaining of the contig from this position is less than txpDist it should be counted as a failure
      // because we couldn't find the path between start and end that is shorter than some txpDist
      else {
        if (remLen > txpDist) {

          if (remLen-txpDist < endContig.contigLen_ &&
              getRemSeq(curContig, remLen-txpDist, isCurContigFw, suffixIfFw) == getRemSeq(endContig, remLen-txpDist, isEndContigFw, prefixIfFw)) {
            appendByLen(seq, curContig, startp, txpDist, isCurContigFw, suffixIfFw);
            return Task::SUCCESS;
          }
          else {
            if(verbose) std::cerr << std::this_thread::get_id() << " "  << "[doBFS] returning failure\n";
            return Task::FAILURE;
          }
        }
        // terminal condition
        // called even when txpDist == 0
        // remLen == txpDist
        auto tmp =  fetchSuccessors(curContig, isCurContigFw, tid, tpos, txpDist);
        if (tmp.size() > 4) {
          std::cerr << "Failure for this f*cking size: " << tmp.size();
          std::exit(1);
        }
        for (auto& c : tmp) {
          pufferfish::util::ContigBlock& cb = (*contigSeqCache_)[c.cid];
          if (cb.contigLen_ - (k - 1) <= endContig.contigLen_ &&
              getRemSeq(cb, cb.contigLen_ - (k - 1), c.isCurContigFw,
                        suffixIfFw) == getRemSeq(endContig,
                                                 cb.contigLen_ - (k - 1),
                                                 isEndContigFw, prefixIfFw)) {
            appendByLen(seq, curContig, startp, txpDist, isCurContigFw,
                        suffixIfFw);
          }
          return Task::SUCCESS;
        }
        return Task::FAILURE; // I'm in the middle of no where!! lost!!
      }
    }

    // if we are here, this means remLen < txpDist

    // If we didn't meet any terminal conditions, we need to dig deeper into the tree through BFS
    // The approach is the same for both valid and dummy end nodes
    //if(verbose) std::cerr << std::this_thread::get_id() << " "  << "[doBFS] remLen < txpDist : " << remLen << " < " << txpDist << "\n" ;
    appendByLen(seq, curContig, startp, remLen, walkForward?isCurContigFw:!isCurContigFw, suffixIfFw);
    txpDist-=remLen;
    if (walkForward) {
      if (isCurContigFw) startp += remLen;
      else startp -= remLen;
      tpos += remLen;
    } else {
      if (isCurContigFw) startp -= remLen;
      else startp += remLen;
      tpos -= remLen;
    }

    std::vector<nextCompatibleStruct> children;
    if (walkForward) {
      children = fetchSuccessors(curContig, isCurContigFw, tid, tpos, txpDist);
    }
    else {
      children = fetchPredecessors(curContig, isCurContigFw, tid, tpos, txpDist);
    }
    for (auto& c : children) {
      // act greedily and return with the first successfully constructed sequence.
      pufferfish::util::ContigBlock& cb = (*contigSeqCache_)[c.cid];
      if (doDFS(tid, c.tpos, c.isCurContigFw, cb, c.cpos, endContig, isEndContigFw, txpDist, walkForward, seq, verbose) == Task::SUCCESS) {
        return Task::SUCCESS;
      }
    }

    // If couldn't find any valid/compatible successors,
    //then this path was a dead end. Revert your append and return with failure
    if(verbose) std::cerr << std::this_thread::get_id() << " "  <<"[doBFS] failed!!\n";
    //revert back
    txpDist+=remLen;
    if (walkForward) {
      if (isCurContigFw) startp -= remLen;
      else startp += remLen;
      tpos -= remLen;
    } else {
      if (isCurContigFw) startp += remLen;
      else startp -= remLen;
      tpos += remLen;
    }

    cutoff(seq, remLen);
    return Task::FAILURE;
}


template <typename PufferfishIndexT>
size_t RefSeqConstructor<PufferfishIndexT>::remainingLen(pufferfish::util::ContigBlock& contig, size_t startp, bool isCurContigFw, bool fromTheEnd, bool verbose) {
  (void) verbose;
      if ( isCurContigFw == fromTheEnd)
      return contig.contigLen_ - startp - 1;
    else
      return startp ;
}


template <typename PufferfishIndexT>
void RefSeqConstructor<PufferfishIndexT>::append(std::string& seq,
                                                 pufferfish::util::ContigBlock& contig,
                                                 size_t startp, size_t endp,
                                                 bool isCurContigFw,
                                                 bool verbose) {

  if(isCurContigFw) {
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[append] 1 " << seq << " clipping by pos: from "<<startp+1<< " to " << endp << " in a contig with len "<<contig.contigLen_
                          << " str len: " << contig.seq.length() << "\n" ;
      seq += contig.substrSeq(startp+1, endp-startp-1);
  }
    else {
      if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[append] 2 rc " << seq << " clipping by pos: from "<<endp+1<< " to " << startp << " in a contig with len "<<contig.contigLen_<<"\n" ;
      // we are always building the seq by moving forward in transcript, so we always append (& never prepend) any substring that we construct
      seq += rc(contig.substrSeq(endp+1, startp-endp-1)); 
    }
}


template <typename PufferfishIndexT>
void RefSeqConstructor<PufferfishIndexT>::appendByLen(std::string& seq, pufferfish::util::ContigBlock& contig, size_t startp, size_t len, bool isCurContigFw, bool appendSuffix, bool verbose) {
  if (len == 0)
    return;
  if (isCurContigFw && appendSuffix) { // append suffix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[appendByLen] 1 from " << startp+1 << " to " << startp+1+len << " total length " << contig.contigLen_ << "\n";
    seq += contig.substrSeq(startp+1, len);
  }
  else if (isCurContigFw && !appendSuffix) {// append prefix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[appendByLen] 2 from " << startp-len << " to " << startp << " total length " << contig.contigLen_ << "\n";
    seq = contig.substrSeq(startp-len, len) + seq;
  }
  else if (!isCurContigFw && appendSuffix) {// append rc of the seq from the other end as a suffix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[appendByLen] 3 rc from " << startp-len << " to " << startp << " total length " << contig.contigLen_ << "\n";
    seq += rc(contig.substrSeq(startp-len, len));
  }
  else if (!isCurContigFw && !appendSuffix) {// append rc of the seq as prefix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[appendByLen] 4 rc from " << startp+1 << " to " << startp+1+len << " total length " << contig.contigLen_ << "\n";
    seq = rc(contig.substrSeq(startp+1, len)) + seq;
  }
}


template <typename PufferfishIndexT>
std::string RefSeqConstructor<PufferfishIndexT>::getRemSeq(pufferfish::util::ContigBlock& contig, size_t len, bool isCurContigFw, bool appendSuffix, bool verbose) {
  std::string seq = "";
  if (len == 0)
    return "";
  if (isCurContigFw && appendSuffix) {// append suffix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[getRemSeq] 1 from " << contig.contigLen_-len << " to " << contig.contigLen_ << " total length " << contig.contigLen_ << "\n";
    seq += contig.substrSeq(contig.contigLen_-len, len);
  }
  else if (isCurContigFw && !appendSuffix) {// append prefix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[getRemSeq] 2 from " << 0 << " to " << len << " total length " << contig.contigLen_ << "\n";
    seq = contig.substrSeq(0, len) + seq;
  }
  else if (!isCurContigFw && appendSuffix) {// append rc of the seq from the other end as a suffix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[getRemSeq] 3 rc from " << 0 << " to " << len << " total length " << contig.contigLen_ << "\n";
    seq += rc(contig.substrSeq(0, len));
  }
  else if (!isCurContigFw && !appendSuffix) {// append rc of the seq as prefix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[getRemSeq] 4 rc from " << contig.contigLen_-len << " to " << contig.contigLen_ << " total length " << contig.contigLen_ << "\n";
    seq = rc(contig.substrSeq(contig.contigLen_-len, len)) + seq;
  }
  return seq;
}

template <typename PufferfishIndexT>
void RefSeqConstructor<PufferfishIndexT>::cutoff(std::string& seq, size_t len, bool verbose) {
  (void) verbose;
    seq = seq.substr(0, seq.length()-len);
}


//TODO couldn't make it reference bc of the section seq += rc(str)
// now eachtime calling this we are copying a string twice which is not good
template <typename PufferfishIndexT>
std::string RefSeqConstructor<PufferfishIndexT>::rc(std::string str, bool verbose) {
  (void) verbose;
    for (uint32_t i = 0; i < str.length()/2; i++) {
      char tmp = str[i];
      str[i] = rev(str[str.length()-1-i]);
      str[str.length()-1-i] = rev(tmp);
    }
    return str;
}


template <typename PufferfishIndexT>
char RefSeqConstructor<PufferfishIndexT>::rev(const char& c) {
    switch(c){
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    }
    return 'N';
}

template <typename PufferfishIndexT>
std::vector<nextCompatibleStruct> RefSeqConstructor<PufferfishIndexT>::fetchSuccessors(pufferfish::util::ContigBlock& contig,
                                                                                       bool isCurContigFw,
                                                                                       size_t tid,
                                                                                       size_t tpos,
                                                                                       uint32_t txpDist,
                                                                                       bool verbose) {

  //if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[Fetch successors] \n" ;
    std::vector<nextCompatibleStruct> successors ;
    CanonicalKmer::k(k) ;

    auto& edges = pfi_->getEdge() ;
    pufferfish::util::Direction dir = isCurContigFw?pufferfish::util::Direction::FORWARD:pufferfish::util::Direction::BACKWORD ;

    uint8_t edgeVec = edges[contig.contigIdx_] ;
    std::vector<pufferfish::util::extension> ext = pufferfish::util::getExts(edgeVec) ;

    if(!ext.empty()){
      CanonicalKmer kb ;
      kb.fromStr(contig.seq.substr(0,k)) ;
      CanonicalKmer ke ;
      ke.fromStr(contig.seq.substr(contig.contigLen_-k,k)) ;
      CanonicalKmer kt ;

      for(auto& ed : ext){
        if(ed.dir == dir){
          (dir == pufferfish::util::Direction::FORWARD)?ke.shiftFw(ed.c):kb.shiftBw(ed.c) ;
          (dir == pufferfish::util::Direction::FORWARD)?kt.fromNum(ke.getCanonicalWord()):kt.fromNum(kb.getCanonicalWord()) ;

          auto nextHit = pfi_->getRefPos(kt) ;

          if(contigSeqCache_->find(nextHit.contigIdx_) == contigSeqCache_->end()){
            std::string contigseq = pfi_->getSeqStr(nextHit.globalPos_,nextHit.contigLen_);
            (*contigSeqCache_)[nextHit.contigIdx_] = {nextHit.contigIdx_, nextHit.globalPos_, nextHit.contigLen_, contigseq} ;
          }

          pufferfish::util::ContigBlock cb = (*contigSeqCache_)[nextHit.contigIdx_] ;
          nextCompatibleStruct theBest {cb.contigIdx_, tpos, 0, true};
          bool isBestValid = false;
          for(auto& posIt: nextHit.refRange){
            size_t succLastBaseTpos = posIt.pos() + cb.contigLen_ - 1;
            int overlap = tpos - posIt.pos() + 1;
            if (posIt.transcript_id() == tid and tpos < succLastBaseTpos and overlap > 0) {
              //heuristic : TODO make sure it's right
              // As a successor, if the txp distance is not less than txpDist, then they are not successors in this txp
              // that 1 base in tpos + txpDist + 1 was kept as an overlap with end contig to have sanity check
              if (theBest.tpos < succLastBaseTpos) {
                isBestValid = true;
                theBest.tpos = succLastBaseTpos; 
                theBest.isCurContigFw = posIt.orientation();
                theBest.cpos=theBest.isCurContigFw?(overlap-1):(cb.contigLen_-overlap);
                if (txpDist < succLastBaseTpos - tpos)
                  break;
              }
            }
          }
          if (isBestValid) {
            theBest.tpos = tpos;
            successors.push_back(theBest) ;
          }
        }
      }
    }

    if (verbose) std::cerr << " RETURNING SUCCESSOR!!\n";
    return successors;
}


template <typename PufferfishIndexT>
std::vector<nextCompatibleStruct> RefSeqConstructor<PufferfishIndexT>::fetchPredecessors(pufferfish::util::ContigBlock& contig,
                                                                                       bool isCurContigFw,
                                                                                       size_t tid,
                                                                                       size_t tpos,
                                                                                         uint32_t txpDist,
                                                                                         bool verbose) {

  //if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[Fetch successors] \n" ;
    std::vector<nextCompatibleStruct> predecessors ;
    CanonicalKmer::k(k) ;

    // [DEV BUG] TODO : Moving to having only one "edge set" breaks this.  FIX THIS!!
    // The below line used to be this:
    //auto& edges = pfi_->getRevEdge() ;
    auto& edges = pfi_->getEdge() ;

    pufferfish::util::Direction dir = isCurContigFw?pufferfish::util::Direction::BACKWORD:pufferfish::util::Direction::FORWARD ;

    uint8_t edgeVec = edges[contig.contigIdx_] ;
    std::vector<pufferfish::util::extension> ext = pufferfish::util::getExts(edgeVec) ;

    if(!ext.empty()){
      CanonicalKmer kb ;
      kb.fromStr(contig.seq.substr(0,k)) ;
      CanonicalKmer ke ;
      ke.fromStr(contig.seq.substr(contig.contigLen_-k,k)) ;
      CanonicalKmer kt ;

      for(auto& ed : ext){
        if(ed.dir == dir){
          (dir == pufferfish::util::Direction::FORWARD)?ke.shiftFw(ed.c):kb.shiftBw(ed.c) ;
          (dir == pufferfish::util::Direction::FORWARD)?kt.fromNum(ke.getCanonicalWord()):kt.fromNum(kb.getCanonicalWord()) ;

          auto nextHit = pfi_->getRefPos(kt) ;

          if(contigSeqCache_->find(nextHit.contigIdx_) == contigSeqCache_->end()){
            std::string contigseq = pfi_->getSeqStr(nextHit.globalPos_,nextHit.contigLen_);
            (*contigSeqCache_)[nextHit.contigIdx_] = {nextHit.contigIdx_, nextHit.globalPos_, nextHit.contigLen_, contigseq} ;
          }

          pufferfish::util::ContigBlock cb = (*contigSeqCache_)[nextHit.contigIdx_] ;
          nextCompatibleStruct theBest {cb.contigIdx_, tpos, 0, true};
          bool isBestValid = false;
          for(auto& posIt: nextHit.refRange){
            size_t predLastBaseTpos = posIt.pos() + cb.contigLen_ - 1;
            int overlap = predLastBaseTpos - tpos + 1;
            if (posIt.transcript_id() == tid and tpos > predLastBaseTpos and overlap > 0) {
              //heuristic : TODO make sure it's right
              // As a successor, if the txp distance is not less than txpDist, then they are not successors in this txp
              // that 1 base in tpos + txpDist + 1 was kept as an overlap with end contig to have sanity check
              if (theBest.tpos > predLastBaseTpos) {
                isBestValid = true;
                theBest.tpos = predLastBaseTpos; // NOTE: careful
                theBest.isCurContigFw = posIt.orientation();
                theBest.cpos=theBest.isCurContigFw?(cb.contigLen_-overlap):(overlap-1);
                if (txpDist < tpos - predLastBaseTpos)
                  break;
              }
            }
          }
          if (isBestValid) {
            theBest.tpos = tpos;
            predecessors.push_back(theBest);
          }
        }
      }
    }

    if (verbose) std::cerr << " RETURNING Predecessors!!\n";
    return predecessors;
}

template class RefSeqConstructor<PufferfishIndex>;
template class RefSeqConstructor<PufferfishSparseIndex>;
template class RefSeqConstructor<PufferfishLossyIndex>;