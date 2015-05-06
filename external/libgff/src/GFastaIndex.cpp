/*
 * GFastaIndex.cpp
 *
 *  Created on: Aug 25, 2010
 *      Author: gpertea
 */

#include "GFastaIndex.h"
#define ERR_FAIDXLINE "Error parsing fasta index line: \n%s\n"
#define ERR_FALINELEN "Error: sequence lines in a FASTA record must have the same length!\n"
void GFastaIndex::addRecord(const char* seqname, uint seqlen, off_t foffs, int llen, int llen_full) {
     GFastaRec* farec=records.Find(seqname);
     if (farec!=NULL) {
          GMessage("Warning: duplicate sequence ID (%s) added to the fasta index! Only last entry data will be kept.\n");
          farec->seqlen=seqlen;
          farec->fpos=foffs;
          farec->line_len=llen;
          farec->line_blen=llen_full;
          }
     else {
         farec=new GFastaRec(seqlen,foffs,llen,llen_full);
         records.Add(seqname,farec);
         farec->seqname=records.getLastKey();
         }
}

int GFastaIndex::loadIndex(const char* finame) { //load record info from existing fasta index
    if (finame==NULL) finame=fai_name;
    if (finame!=fai_name) {
      fai_name=Gstrdup(finame);
      }
    if (fai_name==NULL) GError("Error: GFastaIndex::loadIndex() called with no file name!\n");
    records.Clear();
    haveFai=false;
    FILE* fi=fopen(fai_name,"rb");
    if (fi==NULL) {
       GMessage("Warning: cannot open fasta index file: %s!\n",fai_name);
       return 0;
       }
    GLineReader fl(fi);
    char* s=NULL;
    while ((s=fl.nextLine())!=NULL) {
      if (*s=='#') continue;
      char* p=strchrs(s,"\t ");
      if (p==NULL) GError(ERR_FAIDXLINE,s);
      *p=0; //s now holds the genomic sequence name
      p++;
      uint len=0;
      int line_len=0, line_blen=0;
#ifdef __WIN32__
         long offset=-1;
         sscanf(p, "%d%ld%d%d", &len, &offset, &line_len, &line_blen);
#else
         long long offset=-1;
         sscanf(p, "%d%lld%d%d", &len, &offset, &line_len, &line_blen);
#endif
      if (len==0 || line_len==0 || line_blen==0 || line_blen<line_len)
          GError(ERR_FAIDXLINE,p);
      addRecord(s,len,offset,line_len, line_blen);
      }
    fclose(fi);
    haveFai=(records.Count()>0);
    return records.Count();
}

int GFastaIndex::buildIndex() {
    //this parses the whole fasta file, so it could be slow
    if (fa_name==NULL)
       GError("Error: GFastaIndex::buildIndex() called with no fasta file!\n");
    FILE* fa=fopen(fa_name,"rb");
    if (fa==NULL) {
       GMessage("Warning: cannot open fasta index file: %s!\n",fa_name);
       return 0;
       }
    records.Clear();
    GLineReader fl(fa);
    char* s=NULL;
    uint seqlen=0;
    int line_len=0,line_blen=0;
    bool newSeq=false; //set to true after defline
    off_t newSeqOffset=0;
    int prevOffset=0;
    char* seqname=NULL;
    int last_len=0;
    bool mustbeLastLine=false; //true if the line length decreases
    while ((s=fl.nextLine())!=NULL) {
     if (s[0]=='>') {
        if (seqname!=NULL) {
         if (seqlen==0)
            GError("Warning: empty FASTA record skipped (%s)!\n",seqname);
         else { //seqlen!=0
           addRecord(seqname, seqlen,newSeqOffset, line_len, line_blen);
           }
         }
        char *p=s;
        while (*p > 32) p++;
        *p=0;
        GFREE(seqname);
        seqname=Gstrdup(&s[1]);
        newSeq=true;
        newSeqOffset=fl.getfpos();
        last_len=0;
        line_len=0;
        line_blen=0;
        seqlen=0;
        mustbeLastLine=false;
        } //defline parsing
     else { //sequence line
       int llen=fl.length();
       int lblen=fl.getFpos()-prevOffset;
        if (newSeq) { //first sequence line after defline
          line_len=llen;
          line_blen=lblen;
          }
        else {//next seq lines after first
          if (mustbeLastLine || llen>last_len)
             GError(ERR_FALINELEN);
          if (llen<last_len) mustbeLastLine=true;
          }
        seqlen+=llen;
        last_len=llen;
        newSeq=false;
        } //sequence line
     prevOffset=fl.getfpos();
     }//for each line of the fasta file
    if (seqlen>0)
       addRecord(seqname, seqlen, newSeqOffset, line_len, line_blen);
    GFREE(seqname);
    fclose(fa);
    return records.Count();
}


int GFastaIndex::storeIndex(const char* finame) { //write the hash to a file
    if (records.Count()==0)
       GError("Error at GFastaIndex:storeIndex(): no records found!\n");
    FILE* fai=fopen(finame, "w");
    if (fai==NULL) GError("Error creating fasta index file: %s\n",finame);
    int rcount=storeIndex(fai);
    GFREE(fai_name);
    fai_name=Gstrdup(finame);
    return rcount;
}

int GFastaIndex::storeIndex(FILE* fai) {
  int rcount=0;
  GList<GFastaRec> reclist(true,false,true); //sorted, don't free members, unique
  records.startIterate();
  GFastaRec* rec=NULL;
  while ((rec=records.NextData())!=NULL) {
    reclist.Add(rec);
    }
  //reclist has records sorted by file offset
  for (int i=0;i<reclist.Count();i++) {
#ifdef __WIN32__
    int written=fprintf(fai, "%s\t%d\t%ld\t%d\t%d\n",
            reclist[i]->seqname,reclist[i]->seqlen,(long)reclist[i]->fpos,
              reclist[i]->line_len, reclist[i]->line_blen);
#else
    int written=fprintf(fai, "%s\t%d\t%lld\t%d\t%d\n",
            reclist[i]->seqname, reclist[i]->seqlen, (long long)(reclist[i]->fpos),
              reclist[i]->line_len, reclist[i]->line_blen);
#endif
    if (written>0) rcount++;
       else break; //couldn't write anymore
    }
  fclose(fai);
  haveFai=(rcount>0);
  return rcount;
}
