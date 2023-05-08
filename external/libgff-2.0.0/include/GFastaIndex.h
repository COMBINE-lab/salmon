/*
 * GFaIdx.h
 *
 *  Created on: Aug 25, 2010
 *      Author: gpertea
 */

#ifndef GFAIDX_H_
#define GFAIDX_H_

#include "GHash.hh"
#include "GList.hh"

class GFastaRec {
 public:
  char* seqname;
  uint seqlen;
  off_t fpos;
  int line_len; //effective line length (without EoL)
  int line_blen; //length of line including EoL characters
  GFastaRec(uint slen=0, off_t fp=0, int llen=0, int llenb=0) {
    seqname=NULL; //only a pointer copy
    seqlen=slen;
    fpos=fp;
    line_len=llen;
    line_blen=llenb;
    }
  bool operator==(GFastaRec& d){
      return (fpos==d.fpos);
      }
  bool operator>(GFastaRec& d){
     return (fpos>d.fpos);
     }
  bool operator<(GFastaRec& d){
    return (fpos<d.fpos);
    }

};

class GFastaIndex {
  char* fa_name;
  char* fai_name;
  bool haveFai;
 public:
  GHash<GFastaRec> records;
  void addRecord(const char* seqname, uint seqlen,
                    off_t foffs, int llen, int llen_full);

  GFastaRec* getRecord(const char* seqname) {
    return records.Find(seqname);
    }
  bool hasIndex() { return haveFai; }
  int loadIndex(const char* finame);
  int buildIndex(); //build index in memory by parsing the whole fasta file
  int storeIndex(const char* finame);
  int storeIndex(FILE* fai);
  int getCount() { return records.Count(); }
  GFastaIndex(const char* fname, const char* finame=NULL):records() {
    if (fileExists(fname)!=2) GError("Error: fasta file %s not found!\n",fname);
    if (fileSize(fname)<=0) GError("Error: invalid fasta file %s !\n",fname);
    fa_name=Gstrdup(fname);
    fai_name=finame!=NULL ? Gstrdup(finame) : NULL;
    if (fileSize(fa_name)==0) {
      GError("Error creating GFastaIndex(%s): invalid fasta file!\n",fa_name);
      }
    haveFai=false;
    if (fai_name!=NULL && fileSize(fai_name)>0) {
       //try to load the index file if it exists
       loadIndex(fai_name);
       haveFai=(records.Count()>0);
       }
    }
  ~GFastaIndex() {
    GFREE(fa_name);
    GFREE(fai_name);
    }
};

#endif /* GFAIDX_H_ */
