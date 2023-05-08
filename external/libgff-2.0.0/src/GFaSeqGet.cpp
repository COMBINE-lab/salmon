#include "GFaSeqGet.h"
#include "gdna.h"
#include <ctype.h>

GFaSeqGet* fastaSeqGet(GFastaDb& gfasta, const char* seqid) {
  if (gfasta.fastaPath==NULL) return NULL;
  return gfasta.fetch(seqid);
}


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

GFaSeqGet::GFaSeqGet(const char* faname, uint seqlen, off_t fseqofs, int l_len, int l_blen):fname(NULL),
		fh(NULL), fseqstart(0), seq_len(0), line_len(0),
		line_blen(0), lastsub(NULL), seqname(NULL) {
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

GFaSeqGet::GFaSeqGet(FILE* f, off_t fofs, bool validate):fname(NULL), fh(NULL),
	    fseqstart(0), seq_len(0), line_len(0), line_blen(0),
		lastsub(NULL), seqname(NULL) {
  if (f==NULL) GError("Error (GFaSeqGet) : null file handle!\n");
  fh=f;
  initialParse(fofs, validate);
  lastsub=new GSubSeq();
}

void GFaSeqGet::initialParse(off_t fofs, bool checkall) {
 static const char gfa_ERRPARSE[]="Error (GFaSeqGet): invalid FASTA file format.\n";
 if (fofs!=0) { fseeko(fh,fofs,SEEK_SET); } //e.g. for offsets provided by fasta indexing
 //read the first two lines to determine fasta parameters
 if (seqname) GFREE(seqname);
 GDynArray<char> fseqname(64);
 fseqname.DetachPtr(); //will not free the allocated memory
 fseqstart=fofs;
 int c=getc(fh);
 fseqstart++;
 if (c!='>') //fofs must be at the beginning of a FASTA record!
	 GError("Error (GFaSeqGet): not a FASTA record?\n");

 bool getName=true;
 while ((c=getc(fh))!=EOF) {
   fseqstart++;
   if (getName) {
	   if (c<=32) getName=false;
	   else //seqname.append((char)c);
		   fseqname.Add((char)c);
   }
   if (c=='\n' || c=='\r') { break; } //end of defline
 }
 fseqname.Add('\0'); //terminate the string
 seqname=fseqname(); //takeover the string pointer
 if (c==EOF) GError(gfa_ERRPARSE);
 line_len=0;
 uint lendlen=0;
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
 if (checkall) { //validate the rest of the FASTA records
   uint llen=0; //last line length
   uint elen=0; //length of last line ending
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
     //GMessage("Error (GFaSeqGet): end coordinate (%d) cannot be larger than sequence length %d\n", clen+cstart-1, seq_len);
     //Adjust request:
     clen=seq_len-cstart+1;
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
  if (cstart>seq_len || lastsub->sqstart>cstart) {
	   clen=0; //invalid request
	   return NULL;
  }
  int eol_size=line_blen-line_len;
  char* seqp=lastsub->sq+(int)(cstart-lastsub->sqstart); //should be positive offset?
  //find the proper file offset and read the appropriate lines
  cstart--; //seq start offset, 0-based
  int lineofs = cstart % line_len;
  //file offset, relative to the first letter of the sequence in the file
  off_t f_start= ((int)(cstart/line_len))*line_blen + lineofs;
  uint letters_toread=clen; //actual sequence letters to read
  int maxlen=(seq_len>0)? seq_len-cstart : MAX_FASUBSEQ ;
  if (clen==0) letters_toread=maxlen; //read max allowed, or to the end of file
  uint c_end=cstart+letters_toread; //cstart+clen
  off_t f_end= ((int)(c_end/line_len))*line_blen + c_end % line_len;
  int bytes_toRead=f_end-f_start;
  f_start+=fseqstart; // file offset from the beginning of the file
  fseeko(fh, f_start, SEEK_SET);
  size_t actual_read=0;
  char* smem=NULL;
  GMALLOC(smem, bytes_toRead);
  actual_read=fread((void*)smem, 1, bytes_toRead, fh);
  if (actual_read==0) {
	  //error reading any bytes from the file, or invalid request
	  clen=0;
	  return (const char*)seqp;
  }
  uint mp=0; //current read offset in smem
  uint sublen=0; //current sequence letter storage offset in seqp
  //copySeqOnly(seqp, smem, actualrlen);
  bool rdone=false;
  if (lineofs>0) { //read the partial first line
    uint reqrlen=line_len-lineofs;
    if (reqrlen>letters_toread) {
    	reqrlen=letters_toread; //in case we need to read just a few chars
    	rdone=true;
    }
    if (reqrlen>actual_read) {
    	reqrlen=actual_read; //incomplete file read?
    	rdone=true;
    }
    memcpy((void*)seqp, (void*)smem, reqrlen);
    if (rdone) { //eof reached prematurely
      GFREE(smem);
      clen=reqrlen;
      return (const char*)seqp;
    }
    letters_toread-=reqrlen;
    sublen+=reqrlen;
    mp+=reqrlen+eol_size;
    if (mp>actual_read) {
        GFREE(smem);
        clen=reqrlen;
        return (const char*)seqp;
    }
  }//loading first line
  //read the rest of the lines
  while (letters_toread>=line_len && mp+line_len<actual_read) {
    //char* rseqp=&(seqp[sublen]);
    memcpy((void*)(&seqp[sublen]), (void*)(&smem[mp]), line_len);
    sublen+=line_len;
    letters_toread-=line_len;
    mp+=line_blen;
  }
  if (mp>=actual_read) {
	GFREE(smem);
	clen=sublen;
	return (const char*)seqp;
  }
  // read the last partial line, if any
  if (letters_toread>0) {
    if (mp+letters_toread>actual_read)
    	 letters_toread=actual_read-mp;
    if (letters_toread>0) {
       memcpy((void*)(&seqp[sublen]), (void*)(&smem[mp]), letters_toread);
       sublen+=letters_toread;
    }
  }
  //lastsub->sqlen+=sublen;
  GFREE(smem);
  clen=sublen;
  return (const char*)seqp;
}


