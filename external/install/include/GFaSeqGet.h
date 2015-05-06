#ifndef GFASEQGET_H
#define GFASEQGET_H
#include "GList.hh"

#define MAX_FASUBSEQ 0x20000000
//max 512MB sequence data held in memory at a time

class GSubSeq {
 public:
  uint sqstart; //1-based coord of subseq start on sequence
  uint sqlen;   //length of subseq loaded
  char* sq; //actual subsequence data will be stored here
                // (with end-of-line characters removed)

  /*char* xseq; //the exposed pointer to the last requested subsequence start
  off_t xstart; //the coordinate start for the last requested subseq
  off_t xlen; //the last requested subseq len*/
  GSubSeq() {
     sqstart=0;
     sqlen=0;
     sq=NULL;
     /* xseq=NULL;
     xstart=0;
     xlen=0;*/
     }
  ~GSubSeq() {
     GFREE(sq);
     }
  // genomic, 1-based coordinates:
  void setup(uint sstart, int slen, int sovl=0, int qfrom=0, int qto=0, uint maxseqlen=0);
    //check for overlap with previous window and realloc/extend appropriately
    //returns offset from seq that corresponds to sstart
    // the window will keep extending until MAX_FASUBSEQ is reached
};

class GFaSeqGet {
  char* fname;
  FILE* fh;
  //raw offset in the file where the sequence actually starts:
  off_t fseqstart;
  uint seq_len; //total sequence length, if known (when created from GFastaIndex)
  int line_len; //length of each line of text
  int line_blen; //binary length of each line
                 // = line_len + number of EOL character(s)
  GSubSeq* lastsub;
  void initialParse(off_t fofs=0, bool checkall=true);
  const char* loadsubseq(uint cstart, int& clen);
  void finit(const char* fn, off_t fofs, bool validate);
 public:
  GFaSeqGet() {
    fh=NULL;
    fseqstart=0;
    seq_len=0;
    line_len=0;
    line_blen=0;
    fname=NULL;
    lastsub=NULL;
    }
  GFaSeqGet(const char* fn, off_t fofs, bool validate=false) {
     seq_len=0;
     finit(fn,fofs,validate); 
     }
  GFaSeqGet(const char* fn, bool validate=false) {
     seq_len=0;
     finit(fn,0,validate);
     }

  GFaSeqGet(const char* faname, uint seqlen, off_t fseqofs, int l_len, int l_blen);
  //constructor from GFastaIndex record

  GFaSeqGet(FILE* f, off_t fofs=0, bool validate=false);

  ~GFaSeqGet() {
    if (fname!=NULL) {
       GFREE(fname);
       fclose(fh);
       }
    delete lastsub;
    }
  const char* subseq(uint cstart, int& clen);
  const char* getRange(uint cstart=1, uint cend=0) {
      if (cend==0) cend=(seq_len>0)?seq_len : MAX_FASUBSEQ;
      if (cstart>cend) { Gswap(cstart, cend); }
      int clen=cend-cstart+1;
      //int rdlen=clen;
      return subseq(cstart, clen);
      }

  char* copyRange(uint cstart, uint cend, bool revCmpl=false, bool upCase=false);
  //caller is responsible for deallocating the return string

  void loadall(uint32 max_len=0) {
    //TODO: better read the whole sequence differently here - line by line
    //so when EOF or another '>' line is found, the reading stops!
    int clen=(seq_len>0) ? seq_len : ((max_len>0) ? max_len : MAX_FASUBSEQ);
    subseq(1, clen);
    }
  void load(uint cstart, uint cend) {
     //cache as much as possible
      if (seq_len>0 && cend>seq_len) cend=seq_len; //correct a bad request
      int clen=cend-cstart+1;
      subseq(cstart, clen);
     }
  int getsublen() { return lastsub!=NULL ? lastsub->sqlen : 0 ; }
  off_t getseqofs() { return fseqstart; }
  int getLineLen() { return line_len; }
  int getLineBLen() { return line_blen; }
  //reads a subsequence starting at genomic coordinate cstart (1-based)
 };


#endif
