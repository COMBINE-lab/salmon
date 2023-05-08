#ifndef REF_SEQ_CONSTRUCTOR_HPP
#define REF_SEQ_CONSTRUCTOR_HPP


#include <sparsepp/spp.h>

#include "Util.hpp"
#include "CanonicalKmer.hpp"
#include <stack>

enum Task {
           SUCCESS,
           FAILURE
};


enum SearchType {
  PREDECESSOR,
  SUCCESSOR
};

struct nextCompatibleStruct {
  size_t cid;
  size_t tpos ;
  uint32_t cpos ;
  bool isCurContigFw ;

  nextCompatibleStruct(size_t cntgIn, size_t tposIn, uint32_t cposIn, bool mFw) : cid(cntgIn), tpos(tposIn), cpos(cposIn), isCurContigFw(mFw) {} 
} ;

struct strInfo {
  size_t pid;
  std::string str;
  bool shouldAppend;
};


template <typename PufferfishIndexT>
class RefSeqConstructor {

public:
  RefSeqConstructor(PufferfishIndexT* pfi, spp::sparse_hash_map<uint32_t, pufferfish::util::ContigBlock>* contigSeqCache);
  Task fillSeq(size_t tid,
               size_t tpos,
                                               bool isCurContigFw,
                                               pufferfish::util::ContigBlock& curContig,
                                               uint32_t startp,
                                               pufferfish::util::ContigBlock& endContig,
                                               uint32_t endp,
                                               bool isEndContigFw,
               uint32_t txpDist,
               std::string& refSeq,
               bool verbose=false);
  /* Task doBFS(size_t tid,
             size_t tpos,
             bool isCurContigFw,
             pufferfish::util::ContigBlock& curContig,
             size_t startp,
             pufferfish::util::ContigBlock& endContig,
             bool isEndContigFw,
             uint32_t threshold,
             std::string& seq);
  */
  Task doDFS(size_t tid,
             size_t tpos,
             bool isCurContigFw,
             pufferfish::util::ContigBlock& curContig,
             uint32_t startp,
             pufferfish::util::ContigBlock& endContig,
             bool isEndContigFw,
             uint32_t threshold,
             bool walkForward,
             std::string& refSeq,
             bool verbose=false);


  /*  //search predecessors
  Task fillSeqLeft(size_t tid,
                   size_t tpos,
                   pufferfish::util::ContigBlock& curContig,
                   bool isCurContigFw,
                   uint32_t cstart,
                   std::string& seq);

  Task doRevBFS(size_t tid,
                size_t tpos,
                pufferfish::util::ContigBlock& curContig,
                bool isCurContigFw,
                uint32_t cstart,
                uint32_t threshold,
                std::string& seq);
  */
private:
  PufferfishIndexT* pfi_ ;
  size_t k ;
  spp::sparse_hash_map<uint32_t, pufferfish::util::ContigBlock>* contigSeqCache_;



  size_t remainingLen(pufferfish::util::ContigBlock& contig, size_t startp, bool isCurContigFw, bool fromTheEnd, bool verbose=false);
  void append(std::string& seq, pufferfish::util::ContigBlock& contig, size_t startp, size_t endp, bool isCurContigFw, bool verbose=false);
  void appendByLen(std::string& refSeq, pufferfish::util::ContigBlock& contig, size_t startp, size_t len, bool isCurContigFw, bool appendSuffix, bool verbose=false);
  //TODO 
  //void prependByLen(std::string& seq, pufferfish::util::ContigBlock& contig, size_t startp, size_t len, bool isCurContigFw, bool appendSuffix);
  std::string getRemSeq(pufferfish::util::ContigBlock& contig, size_t len, bool isCurContigFw, bool appendSuffix, bool verbose=false);
  void cutoff(std::string& seq, size_t len, bool verbose=false);
  std::string rc(std::string str, bool verbose=false);
  char rev(const char& c);
  std::vector<nextCompatibleStruct> fetchSuccessors(pufferfish::util::ContigBlock& contig,
                                                 bool isCurContigFw,
                                                 size_t tid,
                                                    size_t tpos,
                                                    uint32_t txpDist,
                                                    bool verbose=false);

  std::vector<nextCompatibleStruct> fetchPredecessors(pufferfish::util::ContigBlock& contig,
                                                    bool isCurContigFw,
                                                    size_t tid,
                                                        size_t tpos,
                                                      uint32_t txpDist,
                                                      bool verbose=false);
};

#endif
