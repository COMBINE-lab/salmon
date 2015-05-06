#include "GFaSeqGet.h"
#include "gdna.h"
#include <ctype.h>

void GSubSeq::setup(uint sstart, int slen, int sovl, int qfrom, int qto, uint maxseqlen) {
     if (sovl==0) {
       GFREE(sq);
       sqstart=sstart;
       uint max_len=(maxseqlen>0) ? maxseqlen : MAX_FASUBSEQ;
       sqlen = (slen==0 ? max_len : slen);
       GMALLOC(sq, sqlen);
       return;
       }
  //overlap -- copy the overlapping region
  char* newsq=NULL;
  GMALLOC(newsq, slen);
  memcpy((void*)&newsq[qto], (void*)&sq[qfrom], sovl);
  GFREE(sq);
  sq=newsq;
  sqstart=sstart;
  sqlen=slen;
}

void GFaSeqGet::finit(const char* fn, off_t fofs, bool validate) {
 fh=fopen(fn,"rb");
 if (fh==NULL) {
   GError("Error (GFaSeqGet) opening file '%s'\n",fn);
   }
 fname=Gstrdup(fn);
 initialParse(fofs, validate);
 lastsub=new GSubSeq();
}

GFaSeqGet::GFaSeqGet(const char* faname, uint seqlen, off_t fseqofs, int l_len, int l_blen) {
//for GFastaIndex use mostly -- the important difference is that
//the file offset is to the sequence, not to the defline
  fh=fopen(faname,"rb");
  if (fh==NULL) {
    GError("Error (GFaSeqGet) opening file '%s'\n",faname);
    }
  fname=Gstrdup(faname);
  line_len=l_len;
  line_blen=l_blen;
  seq_len=seqlen;
  if (line_blen<line_len)
       GError("Error (GFaSeqGet): invalid line length info (len=%d, blen=%d)\n",
              line_len, line_blen);
  fseqstart=fseqofs;
  lastsub=new GSubSeq();
}

GFaSeqGet::GFaSeqGet(FILE* f, off_t fofs, bool validate) {
  fname=NULL;
  fseqstart=0;
  if (f==NULL) GError("Error (GFaSeqGet) : null file handle!\n");
  seq_len=0;
  fh=f;
  initialParse(fofs, validate);
  lastsub=new GSubSeq();
}

void GFaSeqGet::initialParse(off_t fofs, bool checkall) {
 static const char gfa_ERRPARSE[]="Error (GFaSeqGet): invalid FASTA file format.\n";
 if (fofs!=0) { fseeko(fh,fofs,SEEK_SET); } //e.g. for offsets provided by cdbyank
 //read the first two lines to determine fasta parameters
 fseqstart=fofs;
 int c=getc(fh);
 fseqstart++;
 if (c!='>') GError("Error (GFaSeqGet): not a fasta header?\n");
 while ((c=getc(fh))!=EOF) {
   fseqstart++;
   if (c=='\n' || c=='\r') { break; } //end of defline
   }

 if (c==EOF) GError(gfa_ERRPARSE);
 line_len=0;
 int lendlen=0;
 while ((c=getc(fh))!=EOF) {
  if (c=='\n' || c=='\r') { //end of line encountered
     if (line_len>0) { //end of the first "sequence" line
        lendlen++;
        break;
        }
      else {// another EoL char at the end of defline
        fseqstart++;
        continue;
        }
     }// end-of-line characters
  line_len++;
  }
 //we are at the end of first sequence line
 while ((c=getc(fh))!=EOF) {
   if (c=='\n' || c=='\r') lendlen++;
      else {
       ungetc(c,fh);
       break;
       }
   }
 line_blen=line_len+lendlen;
 if (c==EOF) return;
 // -- you don't need to check it all if you're sure it's safe
 if (checkall) { //validate the rest of the FASTA record
   int llen=0; //last line length
   int elen=0; //length of last line ending
   bool waseol=true;
   while ((c=getc(fh))!=EOF) {
     if (c=='>' && waseol) { ungetc(c,fh); break; }
     if (c=='\n' ||  c=='\r') {
        // eol char
        elen++;
        if (waseol) continue; //2nd eol char
        waseol=true;
        elen=1;
        continue;
        }
     if (c<=32) GError(gfa_ERRPARSE); //invalid character encountered
     //--- on a seq char here:
     if (waseol) {//beginning of a seq line
       if (elen && (llen!=line_len || elen!=lendlen))
           //GError(gfa_ERRPARSE);
         GError("Error: invalid FASTA format for GFaSeqGet; make sure that\n\
  the sequence lines have the same length (except for the last line)");
       waseol=false;
       llen=0;
       elen=0;
       }
     llen++;
     } //while reading chars
   }// FASTA checking was requested
 fseeko(fh,fseqstart,SEEK_SET);
}

const char* GFaSeqGet::subseq(uint cstart, int& clen) {
  //cstart is 1-based genomic coordinate within current fasta sequence
   int maxlen=(seq_len>0)?seq_len : MAX_FASUBSEQ;
   //GMessage("--> call: subseq(%u, %d)\n", cstart, clen);
  if (clen>maxlen) {
    GMessage("Error (GFaSeqGet): subsequence cannot be larger than %d\n", maxlen);
    return NULL;
    }
  if (seq_len>0 && clen+cstart-1>seq_len) {
     GMessage("Error (GFaSeqGet): end coordinate (%d) cannot be larger than sequence length %d\n", clen+cstart-1, seq_len);
     }
  if (lastsub->sq==NULL || lastsub->sqlen==0) {
    lastsub->setup(cstart, clen, 0,0,0,seq_len);
    loadsubseq(cstart, clen);
    lastsub->sqlen=clen;
    return (const char*)lastsub->sq;
    }
  //allow extension up to MAX_FASUBSEQ
  uint bstart=lastsub->sqstart;
  uint bend=lastsub->sqstart+lastsub->sqlen-1;
  uint cend=cstart+clen-1;
  int qlen=0; //only the extra len to be allocated/appended/prepended
  uint qstart=cstart; //start coordinate of the new seq block of length qlen to be read from file
  int newlen=0; //the new total length of the buffered sequence lastsub->sq
  int kovl=0;
  int czfrom=0;//0-based offsets for copying a previously read sequence chunk
  int czto=0;
  uint newstart=cstart;
  if (cstart>=bstart && cend<=bend) { //new reg contained within existing buffer
     return (const char*) &(lastsub->sq[cstart-bstart]) ;
    }
  //extend downward
  uint newend=GMAX(cend, bend);
  if (cstart<bstart) { //requested start < old buffer start
    newstart=cstart;
    newlen=(newend-newstart+1);
    if (newlen>MAX_FASUBSEQ) {
       newlen=MAX_FASUBSEQ;
       newend=cstart+newlen-1; //keep newstart, set newend
       }
    qlen=bstart-cstart;
    if (newend>bstart) { //overlap
       if (newend>bend) {// new region is larger & around the old one - so we have two regions to update
         kovl=bend-bstart+1;
         czfrom=0;
         czto=bstart-cstart;
         lastsub->setup(newstart, newlen, kovl, czfrom, czto, seq_len); //this should realloc and copy the kovl subseq
         qlen=bstart-cstart;
         loadsubseq(newstart, qlen);
         qlen=newend-bend;
         int toread=qlen;
         loadsubseq(bend+1, qlen);
         clen-=(toread-qlen);
         lastsub->sqlen=clen;
         return (const char*)lastsub->sq;
         }
        //newend<=bend
       kovl=newend-bstart+1;
       }
     else { //no overlap with previous buffer
       if (newend>bend) kovl=bend-bstart+1;
                   else kovl=newend-bstart+1;
       }
     qlen=bstart-cstart;
     czfrom=0;
     czto=qlen;
    } //cstart<bstart
   else { //cstart>=bstart, possibly extend upwards
    newstart=bstart;
    newlen=(newend-newstart+1);
    if (newlen>MAX_FASUBSEQ) {
       newstart=bstart+(newlen-MAX_FASUBSEQ);//keep newend, assign newstart
       newlen=MAX_FASUBSEQ;
       if (newstart<=bend) { //overlap with old buffer
          kovl=bend-newstart+1;
          czfrom=newstart-bstart;
          czto=0;
          }
       else { //not overlapping old buffer
         kovl=0;
         }
       } //newstart reassigned
    else { //we can extend the buffer to include the old one
      qlen=newend-bend; //how much to read from file
      qstart=bend+1;
      kovl=bend-bstart+1;
      czfrom=0;
      czto=0;
      }
    }
  lastsub->setup(newstart, newlen, kovl, czfrom, czto, seq_len); //this should realloc but copy any overlapping region
  lastsub->sqlen-=qlen; //appending may result in a premature eof
  int toread=qlen;
  loadsubseq(qstart, qlen); //read the missing chunk, if any
  clen-=(toread-qlen);
  lastsub->sqlen+=qlen;
  return (const char*)(lastsub->sq+(cstart-newstart));
}

char* GFaSeqGet::copyRange(uint cstart, uint cend, bool revCmpl, bool upCase) {
  if (cstart>cend) { Gswap(cstart, cend); }
  int clen=cend-cstart+1;
  const char* gs=subseq(cstart, clen);
  if (gs==NULL) return NULL;
  char* r=NULL;
  GMALLOC(r,clen+1);
  r[clen]=0;
  memcpy((void*)r,(void*)gs, clen);
  if (revCmpl) reverseComplement(r,clen);
  if (upCase) {
       for (int i=0;i<clen;i++)
            r[i]=toupper(r[i]);
       }
  return r;
 }

const char* GFaSeqGet::loadsubseq(uint cstart, int& clen) {
  //assumes enough lastsub->sq space allocated previously
  //only loads the requested clen chars from file, at offset &lastsub->sq[cstart-lastsub->sqstart]
  int sofs=cstart-lastsub->sqstart;
  int lendlen=line_blen-line_len;
  char* seqp=lastsub->sq+sofs;
  //find the proper file offset and read the appropriate lines
  uint seqofs=cstart-1;
  uint startlno = seqofs/line_len;
  int lineofs = seqofs % line_len;
  off_t fstart=fseqstart + (startlno*line_blen);
  fstart+=lineofs;

  fseeko(fh, fstart, SEEK_SET);
  int toread=clen;
  int maxlen=(seq_len>0)? seq_len-cstart+1 : MAX_FASUBSEQ ;
  if (toread==0) toread=maxlen; //read max allowed, or to the end of file
  int actualrlen=0;
  int sublen=0;
  if (lineofs>0) { //read the partial first line
    int reqrlen=line_len-lineofs;
    if (reqrlen>toread) reqrlen=toread; //in case we need to read just a few chars
    actualrlen=fread((void*)seqp, 1, reqrlen, fh);
    if (actualrlen<reqrlen) { //eof reached prematurely
      while (seqp[actualrlen-1]=='\n' || seqp[actualrlen-1]=='\r') actualrlen--;
      //check for new sequences in between
      clen=actualrlen;
      sublen+=actualrlen;
      return (const char*)seqp;
      }
    toread-=reqrlen;
    sublen+=reqrlen;
    fseeko(fh, lendlen, SEEK_CUR);
    }
  //read the rest of the lines
  while (toread>=line_len) {
    char* rseqp=&(seqp[sublen]);
    actualrlen=fread((void*)rseqp, 1, line_len, fh);
    /*
    char dbuf[256];dbuf[255]=0;
    strncpy(dbuf,rseqp, actualrlen);
    dbuf[actualrlen]=0;
    GMessage("<<<read line: %s\n",dbuf);
    */
    if (actualrlen<line_len) {
      while (rseqp[actualrlen-1]=='\n' || rseqp[actualrlen-1]=='\r') actualrlen--;
      sublen+=actualrlen;
      clen=sublen;
      return (const char*)seqp;
      }
    toread-=actualrlen;
    sublen+=actualrlen;
    fseeko(fh, lendlen, SEEK_CUR);
    }
  // read the last partial line, if any
  if (toread>0) {
    char* rseqp=&(seqp[sublen]);
    actualrlen=fread((void*)rseqp, 1, toread, fh);
    if (actualrlen<toread) {
      while (rseqp[actualrlen-1]=='\n' || rseqp[actualrlen-1]=='\r')
          actualrlen--;
      }
    sublen+=actualrlen;
    }
  //lastsub->sqlen+=sublen;
  clen=sublen;

  return (const char*)seqp;
  }


