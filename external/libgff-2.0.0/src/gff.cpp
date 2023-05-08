#include "gff.h"

GffNames* GffObj::names=NULL;
//global set of feature names, attribute names etc.
// -- common for all GffObjs in current application!

const uint GFF_MAX_LOCUS = 7000000; //longest known gene in human is ~2.2M, UCSC claims a gene for mouse of ~ 3.1 M
const uint GFF_MAX_EXON  =   30000; //longest known exon in human is ~11K
const uint GFF_MAX_INTRON= 6000000; //Ensembl shows a >5MB mouse intron
const int  GFF_MIN_INTRON = 4; //for mergeCloseExons option
//bool gff_show_warnings = false; //global setting, set by GffReader->showWarnings()
int gff_fid_mRNA=0; //mRNA (has CDS)
int gff_fid_transcript=1; // generic "transcript" feature
int gff_fid_exon=2; // generic "exon"-like feature (exon,CDS,UTR,start/stop codon)
int gff_fid_CDS=3; // CDS feature (CDS, start/stop codon)
const char* exonTypes[]={ "None", "StartCodon", "StopCodon",
		  "CDS", "UTR", "CDS+UTR", "exon" };

const GffScore GFFSCORE_NONE;
//const uint gfo_flag_LEVEL_MSK        = 0x00FF0000;
//const byte gfo_flagShift_LEVEL           = 16;

void gffnames_ref(GffNames* &n) {
  if (n==NULL) n=new GffNames();
  n->numrefs++;
}

void gffnames_unref(GffNames* &n) {
  if (n==NULL) GError("Error: attempt to remove reference to null GffNames object!\n");
  n->numrefs--;
  if (n->numrefs==0) { delete n; n=NULL; }
}

const int CLASSCODE_OVL_RANK = 15;
int classcode_rank(char c) {
	switch (c) {
		case '=': return 0; //intron chain match or full exon chain match if strict matching is enabled
		case '~': return 1; //intron chain match when strict matching is enabled
		case 'c': return 2; //containment, perfect partial match (transfrag < reference)
		case 'k': return 6; // reverse containment (reference < transfrag)
		case 'm': return 6; // full span overlap with all reference introns either matching or retained
		case 'n': return 6; // partial overlap transfrag with at least one intron retention
		case 'j': return 6; // multi-exon transfrag with at least one junction match
		case 'e': return 12; // single exon transfrag partially overlapping an intron of reference (possible pre-mRNA fragment)
		case 'o': return 14; // other generic exon overlap
	//****  >15 = no-overlaps (not on the same strand) from here on *****
		case 's': return 16; //"shadow" - an intron overlaps with a ref intron on the opposite strand (wrong strand mapping?)
		case 'x': return 18; // generic overlap on opposite strand (usually wrong strand mapping)
		case 'i': return 20; // intra-intron (transfrag fully contained within a reference intron)
		case 'y': return 30; // no exon overlap: ref exons fall within transfrag introns!
		case 'p': return 90; //polymerase run
		case 'r': return 92; //repeats
		case 'u': return 94; //intergenic
		case  0 : return 100;
		 default: return 96;
		}
}

const char* strExonType(char xtype) {
	static const char* extbl[7]={"None", "start_codon", "stop_codon", "CDS", "UTR", "CDS_UTR", "exon"};
	if (xtype>0 && xtype<7)
	   return extbl[(int)xtype];
	else return "NULL";
}

int gfo_cmpByLoc(const pointer p1, const pointer p2) {
 GffObj& g1=*((GffObj*)p1);
 GffObj& g2=*((GffObj*)p2);
 if (g1.gseq_id==g2.gseq_id) {
             if (g1.start!=g2.start)
                    return (int)(g1.start-g2.start);
               else if (g1.getLevel()!=g2.getLevel())
                        return (int)(g1.getLevel()-g2.getLevel());
                    else
                        if (g1.end!=g2.end)
                              return (int)(g1.end-g2.end);
                        else return strcmp(g1.getID(), g2.getID());
             }
             else //return (int)(g1.gseq_id-g2.gseq_id); // input order !
            	 return strcmp(g1.getGSeqName(), g2.getGSeqName()); //lexicographic !
}

//comparator for ordering by reference sequence (chromosome) index
int gfo_cmpRefByID(const pointer p1, const pointer p2) {
 GffObj& g1=*((GffObj*)p1);
 GffObj& g2=*((GffObj*)p2);
 if (g1.gseq_id==g2.gseq_id) {
             if (g1.start!=g2.start)
                    return (int)(g1.start-g2.start);
               else if (g1.getLevel()!=g2.getLevel())
                        return (int)(g1.getLevel()-g2.getLevel());
                    else
                        if (g1.end!=g2.end)
                              return (int)(g1.end-g2.end);
                        else return strcmp(g1.getID(), g2.getID());
             }
             else return (g1.gseq_id-g2.gseq_id); // sort refs by their id# order
}

char* GffLine::extractGFFAttr(char* & infostr, const char* oline, const char* attr, bool caseStrict,
		   bool enforce_GTF2, int* rlen, bool deleteAttr) {
 //parse a key attribute and remove it from the info string
 //(only works for attributes that have values following them after ' ' or '=')
 static const char GTF2_ERR[]="Error parsing attribute %s ('\"' required for GTF) at line:\n%s\n";
 int attrlen=strlen(attr);
 char cend=attr[attrlen-1];
 //char* pos = (caseStrict) ? strstr(info, attr) : strifind(info, attr);
 //must make sure attr is not found in quoted text
 char* pos=infostr;
 char prevch=0;
 bool in_str=false;
 bool notfound=true;
 int (*strcmpfn)(const char*, const char*, int) = caseStrict ? Gstrcmp : Gstricmp;
 while (notfound && *pos) {
   char ch=*pos;
   if (ch=='"') {
     in_str=!in_str;
     pos++;
     prevch=ch;
     continue;
     }
   if (!in_str && (prevch==0 || prevch==' ' || prevch == ';')
          && strcmpfn(attr, pos, attrlen)==0) {
      //attr match found
      //check for word boundary on right
      char* epos=pos+attrlen;
      if (cend=='=' || cend==' ' || *epos==0 || *epos==' ') {
        notfound=false;
        break;
        }
      //not a perfect match, move on
      pos=epos;
      prevch=*(pos-1);
      continue;
      }
   //not a match or in_str
   prevch=ch;
   pos++;
   }
 if (notfound) return NULL;
 char* vp=pos+attrlen;
 while (*vp==' ') vp++;
 if (*vp==';' || *vp==0) {
      GMessage("Warning: cannot parse value of GFF attribute \"%s\" at line:\n%s\n", attr, oline);
      return NULL;
 }
 bool dq_enclosed=false; //value string enclosed by double quotes
 if (*vp=='"') {
     dq_enclosed=true;
     vp++;
     }
 if (enforce_GTF2 && !dq_enclosed)
      GError(GTF2_ERR, attr, oline);
 char* vend=vp;
 if (dq_enclosed) {
    while (*vend!='"' && *vend!=';' && *vend!=0) vend++;
    }
 else {
    while (*vend!=';' && *vend!=0) vend++;
    }
 if (enforce_GTF2 && *vend!='"')
     GError(GTF2_ERR, attr, oline);
 char *r=Gstrdup(vp, vend-1);
 if (rlen) *rlen = vend-vp;
 if (deleteAttr) {//-- remove this attribute from infostr
	 while (*vend!=0 && (*vend=='"' || *vend==';' || *vend==' ')) vend++;
	 if (*vend==0) vend--;
	 for (char *src=vend, *dest=pos;;src++,dest++) {
	   *dest=*src; //shift the rest of infostr (copy over)
	   if (*src==0) break;
	 }
 }
 return r;
}
BEDLine::BEDLine(GffReader* reader, const char* l): skip(true), dupline(NULL), line(NULL), llen(0),
		gseqname(NULL), fstart(0), fend(0), strand(0), ID(NULL), info(NULL),
		cds_start(0), cds_end(0), cds_phase(0), exons(1) {
  if (reader==NULL || l==NULL) return;
  llen=strlen(l);
  GMALLOC(line,llen+1);
  memcpy(line, l, llen+1);
  GMALLOC(dupline, llen+1);
  memcpy(dupline, l, llen+1);
  char* t[14];
  int i=0;
  int tidx=1;
  t[0]=line;
  if (startsWith(line, "browser ") || startsWith(line, "track "))
	  return;
  while (line[i]!=0) {
   if (line[i]=='\t') {
    line[i]=0;
    t[tidx]=line+i+1;
    tidx++;
    //our custom BED-13+ format, with GFF3 attributes in 13th column
    if (tidx>12) { info=t[12]; break; }
   }
   i++;
  }
  /* if (tidx<6) { // require BED-6+ lines
   GMessage("Warning: 6+ BED columns expected, instead found:\n%s\n", l);
   return;
   }
  */
  gseqname=t[0];
  char* p=t[1];
  if (!parseUInt(p,fstart)) {
    GMessage("Warning: invalid BED start coordinate at line:\n%s\n",l);
    return;
    }
  ++fstart; //BED start is 0 based
  p=t[2];
  if (!parseUInt(p,fend)) {
    GMessage("Warning: invalid BED end coordinate at line:\n%s\n",l);
    return;
    }
  if (fend<fstart) Gswap(fend,fstart); //make sure fstart<=fend, always
  if (tidx>5) {
	  strand=*t[5];
	  if (strand!='-' && strand !='.' && strand !='+') {
		  GMessage("Warning: unrecognized BED strand at line:\n%s\n",l);
		  return;
	  }
  }
  else strand='.';
  //if (tidx>12) ID=t[12];
  //        else ID=t[3];
  ID=t[3];
  //now parse the exons, if any
  if (tidx>11) {
	  int numexons=0;
	  p=t[9];
	  if (!parseInt(p, numexons) || numexons<=0) {
	      GMessage("Warning: invalid BED block count at line:\n%s\n",l);
	      return;
	  }
	  char** blen;
	  char** bstart;
	  GMALLOC(blen, numexons * sizeof(char*));
	  GMALLOC(bstart, numexons * sizeof(char*));
	  i=0;
	  int b=1;
	  blen[0]=t[10];
	  while (t[10][i]!=0 && b<=numexons) {
		if (t[10][i]==',') {
			t[10][i]=0;
			if (b<numexons)
			  blen[b]=t[10]+i+1;
			b++;
		}
		i++;
	  }
	  b=1;i=0;
	  bstart[0]=t[11];
	  while (t[11][i]!=0 && b<=numexons) {
		if (t[11][i]==',') {
			t[11][i]=0;
			if (b<numexons)
			  bstart[b]=t[11]+i+1;
			b++;
		}
		i++;
	  }
	  GSeg ex;
	  for (i=0;i<numexons;++i) {
		  int exonlen;
		  if (!strToInt(blen[i], exonlen) || exonlen<=0) {
		      GMessage("Warning: invalid BED block size %s at line:\n%s\n",blen[i], l);
		      return;
		  }
		  int exonstart;
		  if (!strToInt(bstart[i], exonstart) || exonstart<0) {
		      GMessage("Warning: invalid BED block start %s at line:\n%s\n",bstart[i], l);
		      return;
		  }
		  if (i==0 && exonstart!=0) {
			  GMessage("Warning: first BED block start is %d>0 at line:\n%s\n",exonstart, l);
			  return;
		  }
		  exonstart+=fstart;
		  uint exonend=exonstart+exonlen-1;
		  if ((uint)exonstart>fend || exonend>fend) {
			  GMessage("Warning: BED exon %d-%d is outside record boundary at line:\n%s\n",exonstart,exonend, l);
			  return;
		  }
		  ex.start=exonstart;ex.end=exonend;
		  exons.Add(ex);
	  }
	  GFREE(blen);
	  GFREE(bstart);
  }
  else { //take it as single-exon transcript
	  GSeg v(fstart, fend);
	  exons.Add(v);
  }
  if (info!=NULL) {
	  char* cdstr=GffLine::extractGFFAttr(info, dupline, "CDS=");
	  if (cdstr) {
		 char* p=strchr(cdstr, ':');
		 if (p!=NULL) {
			*p='\0'; ++p;
		 }
		if (strToUInt(cdstr, cds_start) && cds_start>=fstart-1) {
			++cds_start;
			if (!strToUInt(p, cds_end) || cds_end>fend) {
				GMessage("Warning: invalid CDS (%d-%d) discarded for line:\n%s\n",
						    cds_start, cds_end, dupline);
				cds_start=0;
				cds_end=0; //invalid CDS coordinates
			}
		}
		char* cdstr_phase=NULL;
		if (cds_start>0 && (cdstr_phase=GffLine::extractGFFAttr(info, dupline, "CDSphase="))!=NULL) {
			cds_phase=cdstr_phase[0];
			GFREE(cdstr_phase);
		}
		GFREE(cdstr);
	  }
  }
  if (cds_start==0 && cds_end==0 && tidx>7) {
	//check if columns 7,8 can be reasonably assumed to be CDS start-end coordinates
	if (strToUInt(t[6], cds_start) && strToUInt(t[7], cds_end) && cds_end>cds_start) {
       if (cds_start>=fstart-1 && cds_end<=fend)
    	    cds_start++;
       else { cds_start=0; cds_end=0; }
	}
  }
  skip=false;
}

bool GffLine::parseSegmentList(GVec<GSeg>& segs, char* str) {
	bool segs_valid=false;
	char* p=strchr(str, '-');
	if (p!=NULL && p>str) {
		GDynArray<char*> ss;
		strsplit(str, ss, ',');
		GSeg seg;
		segs_valid=true;
		for (uint i=0;i<ss.Count();++i) {
			char* p=strchr(ss[i], '-');
			if (p==NULL) {
				segs_valid=false;
				break;
			}
			*p='\0'; ++p;
			int xstart=0, xend=0;
			if (!strToInt(ss[i], xstart) || xstart<(int)fstart || xstart>(int)fend){
				segs_valid=false;
				break;
			}
			if (!strToInt(p, xend) || xend<(int)fstart || xend>(int)fend) {
				segs_valid=false;
				break;
			}
			if (xstart>xend) { seg.start=(uint)xend;seg.end=(uint)xstart; }
			else    { seg.start=(uint)xstart;seg.end=(uint)xend; }
			segs.Add(seg);
		} //parse all CDS segments
		if (segs_valid) {
			if (segs.Count()>1) segs.Sort();
		} else
			segs.Clear();
	}
	return segs_valid;
}

GffLine::GffLine(GffReader* reader, const char* l): _parents(NULL), _parents_len(0),
		dupline(NULL), line(NULL), llen(0), gseqname(NULL), track(NULL),
		ftype(NULL), ftype_id(-1), info(NULL), fstart(0), fend(0), //qstart(0), qend(0), qlen(0),
		score(0), score_decimals(-1), strand(0), flags(0), exontype(exgffNone), phase(0), cds_start(0), cds_end(0),
		exons(), cdss(), gene_name(NULL), gene_id(NULL), parents(NULL), num_parents(0), ID(NULL) {
 llen=strlen(l);
 GMALLOC(line,llen+1);
 memcpy(line, l, llen+1);
 GMALLOC(dupline, llen+1);
 memcpy(dupline, l, llen+1);
 skipLine=true; //clear only if we make it to the end of this function
 char* t[9];
 int i=0;
 int tidx=1;
 t[0]=line;
 char fnamelc[128];
 while (line[i]!=0) {
  if (line[i]=='\t') {
   line[i]=0;
   t[tidx]=line+i+1;
   tidx++;
   if (tidx>8) break;
   }
  i++;
  }
 if (tidx<8) { // ignore non-GFF lines
  return;
 }
 gffWarnings=reader->gff_warns;
 gseqname=t[0];
 track=t[1];
 ftype=t[2];
 info=t[8];
 char* p=t[3];
 if (!parseUInt(p,fstart)) {
   //chromosome_band entries in Flybase
   GMessage("Warning: invalid start coordinate at line:\n%s\n",l);
   return;
   }
 p=t[4];
 if (!parseUInt(p,fend)) {
   GMessage("Warning: invalid end coordinate at line:\n%s\n",l);
   return;
   }
 if (fend<fstart) {
	 GMessage("Error: invalid feature coordinates (end<start!) at line:\n%s\n",l);
	 Gswap(fend,fstart); //make sure fstart<=fend, always
	 //return;
 }
 p=t[5];
 if (p[0]=='.' && p[1]==0) {
  score=0;
  score_decimals=-1;
 }
 else {
  score_decimals=0;
  //count decimals
  char* pd=strchr(p,'.');
  if (pd) {
	  ++pd;
	  char* pde=pd;
	  while ((*pde)!=0) ++pde;
	  score_decimals=pde-pd;
  }
  if (!parseFloat(p, score))
       GError("Error parsing feature score from GFF line:\n%s\n",l);

  }
 strand=*t[6];
 if (strand!='+' && strand!='-' && strand!='.')
     GError("Error parsing strand (%c) from GFF line:\n%s\n",strand,l);
 phase=*t[7]; // must be '.', '0', '1' or '2'
 // exon/CDS/mrna filter
 strncpy(fnamelc, ftype, 127);
 fnamelc[127]=0;
 strlower(fnamelc); //convert to lower case
 if (endsWith(fnamelc, "match")) {
   //TODO: do not discard generic cDNA_match/protein_match features, convert them internally
   // (when a hit chain has multiple _match features with the same ID, e.g. multiple HSPs with the same subject)
   // and set GffObj::flag_DISCONTINUOUS
   return;
 }
 bool is_t_data=false;
 bool someRNA=false;
 if (strstr(fnamelc, "utr")!=NULL) {
	 exontype=exgffUTR;
	 is_exon=true;
	 is_t_data=true;
 }
 else if (endsWith(fnamelc, "exon")) {
	 exontype=exgffExon;
	 is_exon=true;
	 is_t_data=true;
 }
 else if (strstr(fnamelc, "stop")!=NULL &&
		 (strstr(fnamelc, "codon")!=NULL || strstr(fnamelc, "cds")!=NULL) &&
		 strstr(fnamelc, "redefined")==NULL && strstr(fnamelc, "selenocysteine")==NULL) {
	 exontype=exgffStopCodon;
	 is_exon=true;
	 is_cds=true; //though some place it outside the last CDS segment
	 is_t_data=true;
 }
 else if (strstr(fnamelc, "start") &&
		 ((strstr(fnamelc, "codon")!=NULL) || strstr(fnamelc, "cds")!=NULL)){
	 exontype=exgffStartCodon;
	 is_exon=true;
	 is_cds=true;
	 is_t_data=true;
 }
 else if (strcmp(fnamelc, "cds")==0) {
	 exontype=exgffCDS;
	 is_exon=true;
	 is_cds=true;
	 is_t_data=true;
 }
 else if (startsWith(fnamelc, "intron") || endsWith(fnamelc, "intron")) {
	 exontype=exgffIntron;
 }
 else if ((someRNA=endsWith(fnamelc,"rna")) || endsWith(fnamelc,"transcript")) { // || startsWith(fnamelc+1, "rna")) {
	 is_transcript=true;
	 is_t_data=true;
	 if (someRNA) ftype_id=GffObj::names->feats.addName(ftype);
 }
 else if (endsWith(fnamelc, "_gene_segment")) {
	 is_transcript=true;
	 is_t_data=true;
	 is_gene_segment=true;
 }
 else if (endsWith(fnamelc, "gene") || startsWith(fnamelc, "gene")) {
	 is_gene=true;
	 is_t_data=true; //because its name will be attached to parented transcripts
 }
 char* Parent=NULL;
 /*
  Rejecting non-transcript lines early if only transcripts are requested ?!
  It would be faster to do this here but there are GFF cases when we reject an
  unusual parent feature here  (e.g. protein with CDS children) and then
  their exon/CDS children show up and get assigned to an implicit parent mRNA
  The solution is to still load this parent as GffObj for now and BAN it later
  so its children get dismissed/discarded as well.
 */
 if (reader->ignoreLocus) {
	 if (strcmp(ftype, "locus")==0) return;
	 if (is_transcript || is_gene) {
		 char* locus=NULL;
        if (reader->is_gff3 || reader->gff_type==0)
        	locus=extractAttr("locus=");
		else locus=extractAttr("locus");
        if (locus!=NULL) { GFREE(locus); }
	 }
 }
 char *gtf_tid=NULL;
 char *gtf_gid=NULL;
 if (reader->is_gff3 || reader->gff_type==0) {
	ID=extractAttr("ID=",true);
	Parent=extractAttr("Parent=",true);
	if (reader->gff_type==0) {
		if (ID!=NULL || Parent!=NULL) reader->is_gff3=true;
			else { //check if it looks like a GTF
				gtf_tid=extractAttr("transcript_id", true, true);
				if (gtf_tid==NULL) {
					gtf_gid=extractAttr("gene_id", true, true);
					if (gtf_gid==NULL) return; //cannot determine file type yet
				}
				reader->is_gtf=true;
			}
	}
 }

 if (reader->is_gff3) {
	 //parse as GFF3
	 //if (ID==NULL && Parent==NULL) return; //silently ignore unidentified/unlinked features
	 if (ID!=NULL) {
		 //has ID attr so it's likely to be a parent feature

		 //look for explicit gene name
		 gene_name=getAttrValue("gene_name=");
		 if (gene_name==NULL) {
			 gene_name=getAttrValue("geneName=");
			 if (gene_name==NULL) {
				 gene_name=getAttrValue("gene_sym=");
				 if (gene_name==NULL) {
					 gene_name=getAttrValue("gene=");
				 }
			 }
		 }
		 gene_id=getAttrValue("geneID=");
		 if (gene_id==NULL) {
			 gene_id=getAttrValue("gene_id=");
		 }
		 /*
		 if (is_gene) { //--WARNING: this might be mislabeled (e.g. TAIR: "mRNA_TE_gene")
			 //---special case: keep the Name and ID attributes of the gene feature
			 //if (gene_name==NULL)
			 //  gene_name=extractAttr("Name=");
			 if (gene_id==NULL) //the ID is also gene_id in this case
				 gene_id=Gstrdup(ID);
			 //skip=false;
			 //return;
			 //-- we don't care about gene parents.. unless it's a mislabeled "gene" feature
		 } //gene feature (probably)
		*/
		 //--parse exons for TLF
		 char* segstr=extractAttr("exons=");
		 bool exons_valid=false;
		 if (segstr) {
			 exons_valid=parseSegmentList(exons, segstr);
			 char* exoncountstr=extractAttr("exonCount=");
			 if (exoncountstr) {
				 int exoncount=0;
				 if (!strToInt(exoncountstr, exoncount) || exoncount!=(int)exons.Count())
					 GMessage("Warning: exonCount attribute value doesn't match the exons attribute!\n");
				 GFREE(exoncountstr);
			 }
			 GFREE(segstr);
		 }
		 if (exons_valid) {
			 bool validCDS=false;
			 segstr=extractAttr("CDS=");
			 if (segstr) {
				 char* p=strchr(segstr, ':');
				 if (p!=NULL) { // CDS=start:end format
					 *p='\0'; ++p;
					 validCDS=true;
					 if (validCDS && strToUInt(segstr, cds_start) && cds_start>=fstart) {
						 if (!strToUInt(p, cds_end) || cds_end>fend) {
							 validCDS=false;
						 }
					 }
					 if (!validCDS || (int)cds_start<=0 || (int)cds_end<=0) {
						 GMessage("Warning: invalid CDS (%d-%d) discarded for line:\n%s\n",
								 cds_start, cds_end, dupline);
						 cds_start=0;
						 cds_end=0;
					 }
				 }  //CDS=start:end format
				 else { //CDS = list of start-end segments, just like the exons
					 validCDS=parseSegmentList(cdss, segstr);
					 if (validCDS && cdss.Count()>0) {
						 if (cds_start==0) cds_start=cdss.First().start;
						 if (cds_end==0) cds_end=cdss.Last().end;
					 }
				 }
				 GFREE(segstr);
			 }
			 if (validCDS) {
				 char* cds_phase=NULL;
				 if ((cds_phase=extractAttr("CDSphase="))!=NULL) {
					 phase=cds_phase[0];
					 GFREE(cds_phase);
				 }
			 } //CDS found
		 }//has valid exons
	 }// has GFF3 ID
	 if (Parent!=NULL) {
		 //keep Parent attr
		 //parse multiple parents
		 num_parents=1;
		 p=Parent;
		 int last_delim_pos=-1;
		 while (*p!=';' && *p!=0) {
			 if (*p==',' && *(p+1)!=0 && *(p+1)!=';') {
				 num_parents++;
				 last_delim_pos=(p-Parent);
			 }
			 p++;
		 }
		 _parents_len=p-Parent+1;
		 _parents=Parent;
		 GMALLOC(parents, num_parents*sizeof(char*));
		 parents[0]=_parents;
		 int i=1;
		 if (last_delim_pos>0) {
			 for (p=_parents+1;p<=_parents+last_delim_pos;p++) {
				 if (*p==',') {
					 char* ep=p-1;
					 while (*ep==' ' && ep>_parents) ep--;
					 *(ep+1)=0; //end the string there
					 parents[i]=p+1;
					 i++;
				 }
			 }
		 }
	 } //has Parent field
	 //special case for gene_id: for genes, this is the ID
	 if (is_gene && gene_id==NULL && ID!=NULL) {
	    	 gene_id=Gstrdup(ID);
	 }
	 //parse other potentially useful GFF3 attributes
	 /*
	 if ((p=strstr(info,"Target="))!=NULL) { //has Target attr
		 p+=7;
		 while (*p!=';' && *p!=0 && *p!=' ') p++;
		 if (*p!=' ') {
			 GError("Error parsing target coordinates from GFF line:\n%s\n",l);
		 }
		 if (!parseUInt(p,qstart))
			 GError("Error parsing target start coordinate from GFF line:\n%s\n",l);
		 if (*p!=' ') {
			 GError("Error parsing next target coordinate from GFF line:\n%s\n",l);
		 }
		 p++;
		 if (!parseUInt(p,qend))
			 GError("Error parsing target end coordinate from GFF line:\n%s\n",l);
	 }
	 if ((p=strifind(info,"Qreg="))!=NULL) { //has Qreg attr
		 p+=5;
		 if (!parseUInt(p,qstart))
			 GError("Error parsing target start coordinate from GFF line:\n%s\n",l);
		 if (*p!='-') {
			 GError("Error parsing next target coordinate from GFF line:\n%s\n",l);
		 }
		 p++;
		 if (!parseUInt(p,qend))
			 GError("Error parsing target end coordinate from GFF line:\n%s\n",l);
		 if (*p=='|' || *p==':') {
			 p++;
			 if (!parseUInt(p,qlen))
				 GError("Error parsing target length from GFF Qreg|: \n%s\n",l);
		 }
	 }//has Qreg attr
	 if (qlen==0 && (p=strifind(info,"Qlen="))!=NULL) {
		 p+=5;
		 if (!parseUInt(p,qlen))
			 GError("Error parsing target length from GFF Qlen:\n%s\n",l);
	 }
	  */
 } //GFF3
 else { // ----------------- GTF syntax ------------------
	 if (reader->transcripts_Only && !is_t_data) {
		 return; //alwasys skip unrecognized non-transcript features in GTF
	 }
	 if (is_gene) {
		 reader->gtf_gene=true;
		 ID = (gtf_tid!=NULL) ? gtf_tid : extractAttr("transcript_id", true, true); //Ensemble GTF might lack this
		 gene_id = (gtf_gid!=NULL) ? gtf_gid : extractAttr("gene_id", true, true);
		 if (ID==NULL) {
			 //no transcript_id -- this should not be valid GTF2 format, but Ensembl (Gencode?)
			 //has being known to add "gene" features with only gene_id in their GTF
			 if (gene_id!=NULL) { //likely a gene feature line (Ensembl!)
			 		 ID=Gstrdup(gene_id); //take over as ID (for defective GTF lacking transcript_id)
			 }
		 }
		 // else if (strcmp(gene_id, ID)==0) //GENCODE v20 gene feature ?
	 }
	 else if (is_transcript) {
		 ID = (gtf_tid!=NULL) ? gtf_tid : extractAttr("transcript_id", true, true);
		//gene_id=extractAttr("gene_id"); // for GTF this is the only attribute accepted as geneID
		 if (ID==NULL) {
			 	 //something is wrong here, cannot parse the GTF ID
				 GMessage("Warning: invalid GTF record, transcript_id not found:\n%s\n", l);
				 return;
		 }
		 gene_id = (gtf_gid!=NULL) ? gtf_gid : extractAttr("gene_id", true, true);
		if (gene_id!=NULL)
			Parent=Gstrdup(gene_id);
		reader->gtf_transcript=true;
		is_gtf_transcript=1;
	 } else { //must be an exon type
		 Parent = (gtf_tid!=NULL) ? gtf_tid : extractAttr("transcript_id", true, true);
		 gene_id = (gtf_gid!=NULL) ? gtf_gid : extractAttr("gene_id", true, true); // for GTF this is the only attribute accepted as geneID
		 //old pre-GTF2 formats like Jigsaw's (legacy support)
		 if (Parent==NULL && exontype==exgffExon) {
			 if (startsWith(track,"jigsaw")) {
				 is_cds=true;
				 strcpy(track,"jigsaw");
				 p=strchr(info,';');
				 if (p==NULL) { Parent=Gstrdup(info); info=NULL; }
				 else { Parent=Gstrdup(info,p-1);
				 info=p+1;
				 }
			 }
		 }
		 if (Parent==NULL) {
			 	 //something is wrong here couldn't parse the transcript ID for this feature
				 GMessage("Warning: invalid GTF record, transcript_id not found:\n%s\n", l);
				 return;
		 }
	 }
	 //more GTF attribute parsing
	 if (is_gene && gene_id==NULL && ID!=NULL)
    	 gene_id=Gstrdup(ID);
	 gene_name=getAttrValue("gene_name");
	 if (gene_name==NULL) {
		 gene_name=getAttrValue("gene_sym");
		 if (gene_name==NULL) {
			 gene_name=getAttrValue("gene");
			 if (gene_name==NULL)
				 gene_name=getAttrValue("genesymbol");
		 }
	 }
	 //*** IMPORTANT: prepare GTF for easy parseAttr by adding '=' character after the attribute name
	 //          for ALL attributes
	 p=info;
	 bool noed=true; //not edited after the last delim
	 bool nsp=false; //non-space found after last delim
	 while (*p!=0) {
		 if (*p==' ') {
			 if (nsp && noed) {
				 *p='=';
				 noed=false;
				 p++;
				 continue;
			 }
		 }
		 else nsp=true; //non-space
		 if (*p==';') { noed=true; nsp=false; }
		 p++;
	 }
	 //-- GTF prepare parents[] if Parent found
	 if (Parent!=NULL) { //GTF transcript_id found as a parent
		 _parents=Parent;
		 num_parents=1;
		 _parents_len=strlen(Parent)+1;
		 GMALLOC(parents, sizeof(char*));
		 parents[0]=_parents;
	 }
 } //GTF


 if (ID==NULL && parents==NULL) {
	 if (gffWarnings)
		 GMessage("Warning: discarding unrecognized feature (no ID or Parent):\n%s\n",dupline);
	 return; //skip
 }
 skipLine=false;
}

//FIXME - this should only be used AFTER finalize() was called, and must have cdss=NULL of course
void GffObj::setCDS(uint cd_start, uint cd_end, char phase) {
  if (cd_start<this->start) {
	  GMessage("Warning: setCDS() called for %s with an out of bounds CDS start %d!\n",
			  gffID, cd_start);
	  return;
  }
  if (cd_end>this->end) {
	  GMessage("Warning: setCDS() called for %s with an out of bounds CDS end %d!\n",
			  gffID, cd_end);
	  return;
  }
  this->CDstart=cd_start;
  this->CDend=cd_end;
  this->CDphase=phase;
  isTranscript(true);
  subftype_id=gff_fid_exon;
  if (monoFeature()) {
     if (exons.Count()==0) addExon(this->start, this->end, exgffExon);
            else exons[0]->exontype=exgffExon;
  }
}

void GffObj::setCDS(GffObj* t) {
	//copy the cdss as well
	uint cd_start=t->CDstart;
	uint cd_end=t->CDend;
	uint phase=t->CDphase;
	if (cd_start<this->start) {
		  GMessage("Warning: setCDS() called for %s with an out of bounds CDS start %d!\n",
				  gffID, cd_start);
		  return;
	}
	if (cd_end>this->end) {
		  GMessage("Warning: setCDS() called for %s with an out of bounds CDS end %d!\n",
				  gffID, cd_end);
		  return;
	}
	this->CDstart=cd_start;
	this->CDend=cd_end;
	this->CDphase=phase;
	isTranscript(true);
	subftype_id=gff_fid_exon;
	if (monoFeature()) {
	   if (exons.Count()==0) addExon(this->start, this->end, exgffExon);
	          else exons[0]->exontype=exgffExon;
	}
	if (t->cdss!=NULL) {
       if (this->cdss!=NULL) delete cdss;
       cdss=new GList<GffExon>(true, true, false);
       for (int i=0;i<t->cdss->Count();i++) {
    	   cdss->Add(new GffExon(*(t->cdss->Get(i))));
       }
	}
}

int GffObj::readExon(GffReader& reader, GffLine& gl) {
  // -- this should only be called before ::finalize()!
  //should make sure to get the right subftype_id!
  if (!isTranscript() && gl.exontype>exgffNone) {
     //subfeature recognized as exon-like, so this should be considered a transcript!
     isTranscript(true);
  }
  if (isTranscript()) {
     if (subftype_id<0) {//exon_ftype_id=gff_fid_exon;
          if (gl.exontype>exgffNone) subftype_id=gff_fid_exon;
                         else subftype_id=names->feats.addName(gl.ftype);
     }
     //any recognized exon-like segment gets the generic "exon" type (also applies to CDS)
     if (gl.exontype==exgffNone && !gl.is_transcript) {
          //extraneous mRNA feature, discard
          if (reader.gff_warns)
            GMessage("Warning: discarding unrecognized transcript subfeature '%s' of %s\n",
                gl.ftype, gffID);
          return -1;
          }
  }
  else { //non-mRNA parent feature, check this subf type
    int subf_id=names->feats.addName(gl.ftype);
    if (subftype_id<0 || exons.Count()==0) //never assigned a subfeature type before (e.g. first exon being added)
       subftype_id=subf_id;
    else {
       if (subftype_id!=subf_id) {
         if (subftype_id==ftype_id && exons.Count()==1 && exons[0]->start==start && exons[0]->end==end) {
            //the existing exon was just a dummy one created by default, discard it?
            exons.Clear();
            covlen=0;
            subftype_id=subf_id; //allow the new subfeature to completely takeover
            }
         else { //multiple subfeatures, prefer those exon-like
             if (reader.gff_warns)
               GMessage("Warning: multiple subfeatures (%s and %s) found for %s, discarding ",
                  names->feats.getName(subf_id), names->feats.getName(subftype_id),gffID);
            if (gl.exontype>exgffNone) { //new feature is an exon, discard previously parsed subfeatures
               if (reader.gff_warns) GMessage("%s.\n", names->feats.getName(subftype_id));
               subftype_id=subf_id;
               exons.Clear();
               covlen=0;
            }
              else { //discard new feature
               if (reader.gff_warns) GMessage("Warning: skipping subfeature %s.\n", names->feats.getName(subf_id));
               return -1; //skip this 2nd subfeature type for this parent!
            }
         }
       } //incoming subfeature is of different type
    } //new subfeature type
  } //non-mRNA parent
  int eidx=-1;
  GList<GffExon>* segs=NULL; //either cds or &exons
  if (gl.is_cds) {
     if (cdss==NULL)
       cdss=new GList<GffExon>(true, true, false);
     segs=cdss;
  } else {
     segs=&exons;
  }
  eidx=addExon(*segs, gl);
  if (eidx<0) {
	  GMessage("Warning: addExon() failed for GFF line:\n%s\n",gl.dupline);
	  return eidx; //this should never happen!
  }
  if (reader.keep_Attrs) {
     if (reader.noExonAttrs) {
           parseAttrs(attrs, gl.info, true);
     }
     else { //need all exon-level attributes
         parseAttrs((*segs)[eidx]->attrs, gl.info, true, gl.is_cds);
     }
  }
  return eidx;
}

int GffObj::addExon(GList<GffExon>& segs, GffLine& gl, int8_t exontype_override) {
	int ex_type=(exontype_override!=exgffNone) ? exontype_override : gl.exontype;
	GffScore exon_score(gl.score, gl.score_decimals);
	int eidx=addExon(gl.fstart, gl.fend, ex_type, gl.phase, exon_score, &segs);
	if (&segs==cdss && isGene() && gl.ID!=NULL && eidx>=0) {
     //special NCBI cases where CDS can be treated as discontiguous features, grouped by their ID
	 //-- used for genes with X_gene_segment features
	 //char* cds_id=Gstrdup(gl.ID);
	 //segs[eidx]->uptr=cds_id;
	 segs[eidx]->uptr=gl.ID;
	 gl.ID=NULL;
	}
	return eidx;
}

int GffObj::exonOverlapIdx(GList<GffExon>& segs, uint s, uint e, int* ovlen, int start_idx) {
	//return the exons' index for the overlapping OR ADJACENT exon
	//ovlen, if given, will return the overlap length
	//if (s>e) Gswap(s,e);
	for (int i=start_idx;i<segs.Count();i++) {
		if (segs[i]->start>e+1) break;
		if (s-1>segs[i]->end) continue;
		//-- overlap/adjacent if we are here:
		if (ovlen!=NULL) {
			int ovlend= (segs[i]->end>e) ? e : segs[i]->end;
			*ovlen= ovlend - ((s>segs[i]->start)? s : segs[i]->start)+1;
		}
		return i;
	} //for each exon
	*ovlen=0;
	return -1;
}

void GffObj::transferCDS(GffExon* cds) {
	//direct adding of a cds to the cdss pointer, without checking
	 if (cdss==NULL) cdss=new GList<GffExon>(true, true, false);
	 cdss->Add(cds); //now the caller must forget this exon!
	 if (CDstart==0 || CDstart>cds->start) CDstart=cds->start;
}

int GffObj::addExon(uint segstart, uint segend, int8_t exontype, char phase, GffScore exon_score, GList<GffExon>* segs) {
   if (segstart>segend) { Gswap(segstart, segend); }
   if (segs==NULL) segs=&exons;
	if (exontype!=exgffNone) { //check for overlaps between exon/CDS-type segments
		//addExonSegment(gl.fstart, gl.fend, gl.score, gl.phase, gl.is_cds, exontype_override);
		int ovlen=0;
		int oi=-1;
	    while ((oi=exonOverlapIdx(*segs, segstart, segend, &ovlen, oi+1))>=0) {
	        //note: ovlen==0 for adjacent segments
		    if ((*segs)[oi]->exontype>exgffNone &&
		    		(*segs)[oi]->start<=segstart && (*segs)[oi]->end>=segend) {
		    		//existing feature contains this segment, so we do NOT need to add it
		    	    //-- unless its the annoying NCBI exception: gene with multiple alternate
		    	    //        _gene_segment CDS features!
		    		if (!(this->isGene() && exontype==exgffCDS &&
		    				(*segs)[oi]->exontype==exgffCDS ))
		    			return oi;
		    }
		    if (ovlen==0 || !(exontype==exgffCDS && (*segs)[oi]->exontype==exgffCDS)) {
		    	//always merge adjacent features
		    	//but NEVER merge two overlapping CDS (CDS programmed ribosomal shift aware)
		    	int8_t segtype=((*segs)[oi]->exontype==exgffCDS || exontype==exgffCDS) ? exgffCDS : exgffExon;
		    	//if expanded upward, may overlap the segment(s) above
		    	expandSegment(*segs, oi, segstart, segend, segtype);
		    	return oi;
		    }
	    }
	} //exon overlap/adjacent check
   //new exon/CDS, not merged in a previous one
   GffExon* enew=new GffExon(segstart, segend, exontype, phase, exon_score.score, exon_score.precision);
   int eidx=segs->Add(enew);
   if (eidx<0) {
    //this would actually be possible if the object is a "Gene" and "exons" are in fact isoforms
     delete enew;
     hasErrors(true);
     return -1;
   }
   if (start>segs->First()->start) start=segs->First()->start;
   if (end<segs->Last()->end) end=segs->Last()->end;
   if (isFinalized() && segs==&exons) {
	   covlen+=(int)(exons[eidx]->end-exons[eidx]->start)+1;
   }
   return eidx;
}

void GffObj::expandSegment(GList<GffExon>& segs, int oi, uint segstart, uint segend, int8_t exontype) {
  //oi is the index of the *first* overlapping segment found that must be enlarged
  covlen-=segs[oi]->len();
  if (segstart<segs[oi]->start)
	  segs[oi]->start=segstart;
  //if (qs && qs<exons[oi]->qstart) exons[oi]->qstart=qs;
  if (segend>segs[oi]->end)
	  segs[oi]->end=segend;
  //if (qe && qe>exons[oi]->qend) exons[oi]->qend=qe;
  //warning: score cannot be properly adjusted! e.g. if it's a p-value it's just going to get worse
  //if (sc!=0) segs[oi]->score=sc;
  //covlen+=exons[oi]->len();
  //if (exons[oi]->exontype< exontype) -- always true
  segs[oi]->exontype = exontype;
  //if (exontype==exgffCDS) exons[oi]->phase=fr;
  //we must check if any more exons are also overlapping this
  int ni=oi+1; //next exon index after oi
  while (ni<segs.Count() && segs[ni]->start<=segend+1) { // next segment overlaps OR adjacent to newly enlarged segment
	  if (segs[ni]->exontype>0 &&
		(segs[ni]->start==segend+1 || segs[ni]->exontype!=exgffCDS || exontype!=exgffCDS)) {
         if (segs[ni]->start<segs[oi]->start) {
        	 segs[oi]->start=segs[ni]->start;
        	 if (strand=='+') segs[oi]->phase=segs[ni]->phase;
         }
         if (segs[ni]->end>segs[oi]->end) {
        	 segs[oi]->end=segs[ni]->end;
        	 if (strand=='-') segs[oi]->phase=segs[ni]->phase;
         }

         segs.Delete(ni);
      } else ++ni;
  } //until no more overlapping/adjacent segments found
  // -- make sure any other related boundaries are updated:
  if (isFinalized()) {
	  if (&segs==&exons) {
		start=exons.First()->start;
		end=exons.Last()->end;
		//recalculate covlen
		covlen=0;
		for (int i=0;i<exons.Count();++i) covlen+=exons[i]->len();
	  }
  }
  else {
	  if (start>segs.First()->start) start=segs.First()->start;
	  if (end<segs.Last()->end) end=segs.Last()->end;
  }
}

void GffObj::removeExon(int idx) {
  if (idx<0 || idx>=exons.Count()) return;
  int segstart=exons[idx]->start;
  int segend=exons[idx]->end;
  exons.Delete(idx);
  if (isFinalized()) {
    covlen -= (int)(segend-segstart)+1;
    start=exons.First()->start;
    end=exons.Last()->end;
    if (isCDSOnly()) { CDstart=start; CDend=end; }
  }
}

void GffObj::removeExon(GffExon* p) {
	for (int idx=0;idx<exons.Count();idx++) {
		if (exons[idx]==p) {
			int segstart=exons[idx]->start;
			int segend=exons[idx]->end;
			exons.Delete(idx);
			covlen -= (int)(segend-segstart)+1;

			if (exons.Count() > 0) {
				start=exons.First()->start;
				end=exons.Last()->end;
				if (isCDSOnly()) { CDstart=start; CDend=end; }
			}
			return;
		}
	}
}

GffObj::GffObj(GffReader& gfrd, BEDLine& bedline):GSeg(0,0),
		exons(true,true,false), cdss(NULL), gscore() {
	uptr=NULL;
	ulink=NULL;
	parent=NULL;
	udata=0;
	flags=0;
	CDstart=0;
	CDend=0;
	CDphase=0;
	attrs=NULL;
	gffID=NULL;
	track_id=-1;
	gseq_id=-1;
	//ftype_id=-1;
	//subftype_id=-1;
	strand='.';
	gffnames_ref(names);
	//qlen=0;qstart=0;qend=0;
	covlen=0;
	geneID=NULL;
	gene_name=NULL;
	ftype_id=gff_fid_transcript;
	subftype_id=gff_fid_exon;
	start=bedline.fstart;
	end=bedline.fend;
	gseq_id=names->gseqs.addName(bedline.gseqname);
	track_id=names->tracks.addName("BED");
	strand=bedline.strand;
	//setup flags from gffline
	isGene(false);
	isTranscript(true);
	gffID=Gstrdup(bedline.ID);
	for (int i=0;i<bedline.exons.Count();++i) {
		int eidx=this->addExon(bedline.exons[i].start, bedline.exons[i].end, exgffExon);
	    if (eidx<0 && gfrd.showWarnings())
	       GMessage("Warning: failed adding segment %d-%d for %s (discarded)!\n",
	    		   bedline.exons[i].start, bedline.exons[i].end, gffID);

	}
	if (bedline.cds_start>0) {
		CDstart=bedline.cds_start;
		CDend=bedline.cds_end;
		if (CDstart>0 && bedline.cds_phase)
			CDphase=bedline.cds_phase;
	}
	if (gfrd.keep_Attrs && bedline.info!=NULL) this->parseAttrs(attrs, bedline.info);
}

GffObj::GffObj(GffReader &gfrd, GffLine& gffline):
     GSeg(0,0), exons(true,true,false), cdss(NULL), children(1,false), gscore() {
  uptr=NULL;
  ulink=NULL;
  parent=NULL;
  udata=0;
  flags=0;
  CDstart=0;
  CDend=0;
  CDphase=0;
  geneID=NULL;
  gene_name=NULL;
  attrs=NULL;
  gffID=NULL;
  track_id=-1;
  gseq_id=-1;
  //ftype_id=-1;
  subftype_id=-1;
  strand='.';
  gffnames_ref(names);
  //qlen=0;qstart=0;qend=0;
  covlen=0;
  ftype_id=gffline.ftype_id;
  start=gffline.fstart;
  end=gffline.fend;
  gseq_id=names->gseqs.addName(gffline.gseqname);
  track_id=names->tracks.addName(gffline.track);
  strand=gffline.strand;
  /*
  qcov=0;
  qlen=gffline.qlen;
  qstart=gffline.qstart;
  qend=gffline.qend;
  */
  //setup flags from gffline
  isCDSOnly(gffline.is_cds); //for now
  isGene(gffline.is_gene);
  isTranscript(gffline.is_transcript || gffline.exontype!=exgffNone);
  //fromGff3(gffline.is_gff3);
  isGeneSegment(gffline.is_gene_segment);
  if (gffline.parents!=NULL && !gffline.is_transcript) {
    //GTF style -- create a GffObj directly by subfeature
    //(also possible orphan GFF3 exon line, or an exon given before its parent (chado))
    if (gffline.exontype!=exgffNone) { //recognized exon-like feature
       ftype_id=gff_fid_transcript; //so this is some sort of transcript
       subftype_id=gff_fid_exon; //subfeatures MUST be exons
       //typical GTF2 without "transcript" line
       // or malformed GFF3 with orphan or premature exon feature (found before the transcript line)
       gffID=Gstrdup(gffline.parents[0]);
       this->createdByExon(true);
       if (gfrd.is_gff3 && gfrd.showWarnings())
     	   GMessage("Warning: exon feature found before transcript ID %s\n",gffID);
       //this is the first exon/segment of the transcript
       readExon(gfrd, gffline);
    }
    else {//unrecognized (non-exon) subfeature
       //make this GffObj of the same feature type
       ftype_id=names->feats.addName(gffline.ftype);
       if (gffline.ID!=NULL) { //unrecognized non-exon feature ? use the ID instead
            this->hasGffID(true);
            gffID=Gstrdup(gffline.ID);
            if (gfrd.keep_Attrs) this->parseAttrs(attrs, gffline.info);
       }
       else { //no ID, just Parent
           GMessage("Warning: unrecognized parented feature without ID found before its parent:\n%s\n", gffline.dupline);
           gffID=Gstrdup(gffline.parents[0]);
           this->createdByExon(true);
           readExon(gfrd, gffline);
       }
    } //unrecognized (non-exon) feature
  } //non-transcript parented subfeature given directly
  else {
    //non-parented feature OR a recognizable transcript
    //create a parent feature in its own right
    gscore.score=gffline.score;
    gscore.precision=gffline.score_decimals;
    if (gffline.ID==NULL || gffline.ID[0]==0)
      GError("Error: no valid ID found for GFF record\n");
    this->hasGffID(true);
    gffID=Gstrdup(gffline.ID); //there must be an ID here
    //if (gffline.is_transcript) ftype_id=gff_fid_mRNA;
      //else
    if (gffline.is_transcript) {
      subftype_id=gff_fid_exon;
      if (ftype_id<0)
        ftype_id=names->feats.addName(gffline.ftype);
      if (gfrd.is_gff3) {
		  if (gffline.exons.Count()>0) {
			  //for compact GFF-like transcript line format (TLF), exons were already found as attributes
				for (int i=0;i<gffline.exons.Count();++i) {
					int eidx=this->addExon(gffline.exons[i].start, gffline.exons[i].end, exgffExon, '.', gscore);
				    if (eidx<0 && gfrd.showWarnings())
				       GMessage("Warning: failed adding exon %d-%d for %s (discarded)!\n",
				    		   gffline.exons[i].start, gffline.exons[i].end, gffID);

				}
		  }
		  if (gffline.cds_start>0) {
				CDstart=gffline.cds_start;
				CDend=gffline.cds_end;
		  }
		  if (gffline.phase!=0) CDphase=gffline.phase;
		  if (gffline.cdss.Count()>0) {
			    //for compact GFF-like transcript line format (TLF), CDS might be already found as attributes
			    if (cdss==NULL) cdss=new GList<GffExon>(true, true, false);
				for (int i=0;i<gffline.cdss.Count();++i) {
					int eidx=this->addExon(gffline.cdss[i].start, gffline.cdss[i].end, exgffCDS, 0, GFFSCORE_NONE, cdss);
				    if (eidx<0 && gfrd.showWarnings())
				       GMessage("Warning: failed adding CDS segment %d-%d for %s (discarded)!\n",
				    		   gffline.cdss[i].start, gffline.cdss[i].end, gffID);

				}
		  }
      }
    } //is_transcript
    if (gfrd.keep_Attrs) this->parseAttrs(attrs, gffline.info);
    if (gfrd.is_gff3 && gffline.parents==NULL && gffline.exontype!=exgffNone) {
       //special case with bacterial genes just given as a CDS/exon, without parent!
       this->createdByExon(true);
       if (ftype_id<0) ftype_id=gff_fid_mRNA;
       readExon(gfrd, gffline);
    }
    if (ftype_id<0)
        ftype_id=names->feats.addName(gffline.ftype);
  }//no parent OR recognizable transcript

  if (gffline.gene_name!=NULL) {
     gene_name=Gstrdup(gffline.gene_name);
     }
  if (gffline.gene_id) { //only for gene features or GTF2 gene_id attribute
     geneID=Gstrdup(gffline.gene_id);
  }
  /*//we cannot assume parents[0] is a gene! for NCBI miRNA, parent can be a primary_transcript feature!
  else if (gffline.is_transcript && gffline.parents!=NULL) {
	 geneID=Gstrdup(gffline.parents[0]);
  }
  */
}

BEDLine* GffReader::nextBEDLine() {
 if (bedline!=NULL) return bedline; //caller should free gffline after processing
 while (bedline==NULL) {
	int llen=0;
	buflen=GFF_LINELEN-1;
	char* l=fgetline(linebuf, buflen, fh, &fpos, &llen);
	if (l==NULL) return NULL;
	int ns=0; //first nonspace position
	while (l[ns]!=0 && isspace(l[ns])) ns++;
	if (l[ns]=='#' || llen<7) continue;
	bedline=new BEDLine(this, l);
	if (bedline->skip) {
	  delete bedline;
	  bedline=NULL;
	  continue;
	}
 }
 return bedline;
}


GffLine* GffReader::nextGffLine() {
 if (gffline!=NULL) return gffline; //caller should free gffline after processing
 while (gffline==NULL) {
    int llen=0;
    buflen=GFF_LINELEN-1;
    char* l=fgetline(linebuf, buflen, fh, &fpos, &llen);
    if (l==NULL) {
         return NULL; //end of file
         }
#ifdef CUFFLINKS
     _crc_result.process_bytes( linebuf, llen );
#endif
    int ns=0; //first nonspace position
    bool commentLine=false;
    while (l[ns]!=0 && isspace(l[ns])) ns++;
    if (l[ns]=='#') {
    	commentLine=true;
    	if (llen<10) {
    		if (commentParser!=NULL) (*commentParser)(l, &gflst);
    		continue;
    	}
    }
    gffline=new GffLine(this, l);
    if (gffline->skipLine) {
       if (commentLine && commentParser!=NULL) (*commentParser)(gffline->dupline, &gflst);
       delete gffline;
       gffline=NULL;
       continue;
    }
    if (gffline->ID==NULL && gffline->parents==NULL)  { //it must have an ID
        //this might not be needed, already checked in the GffLine constructor
        if (gff_warns)
            GMessage("Warning: malformed GFF line, no parent or record Id (kipping\n");
        delete gffline;
        gffline=NULL;
        //continue;
        }
    }
return gffline;
}


char* GffReader::gfoBuildId(const char* id, const char* ctg) {
//caller must free the returned pointer
 char* buf=NULL;
 int idlen=strlen(id);
 GMALLOC(buf, idlen+strlen(ctg)+2);
 strcpy(buf, id);
 buf[idlen]='~';
 strcpy(buf+idlen+1, ctg);
 return buf;
}

GffObj* GffReader::gfoAdd(GffObj* gfo) {
 GPVec<GffObj>* glst=phash.Find(gfo->gffID);
 if (glst==NULL)
	 glst=new GPVec<GffObj>(false);
 int i=glst->Add(gfo);
 phash.Add(gfo->gffID, glst);
 return glst->Get(i);
}

GffObj* GffReader::gfoAdd(GPVec<GffObj>& glst, GffObj* gfo) {
 int i=glst.Add(gfo);
 return glst[i];
}

GffObj* GffReader::gfoReplace(GPVec<GffObj>& glst, GffObj* gfo, GffObj* toreplace) {
 for (int i=0;i<glst.Count();++i) {
	 if (glst[i]==toreplace) {
		 //glst.Put(i,gfo);
		 glst[i]=gfo;
		 break;
	 }
 }
 return gfo;
}

bool GffReader::pFind(const char* id, GPVec<GffObj>*& glst) {
	glst = phash.Find(id);
	return (glst!=NULL);
}

GffObj* GffReader::gfoFind(const char* id, GPVec<GffObj>*& glst,
		const char* ctg, char strand, uint start, uint end) {
	GPVec<GffObj>* gl=NULL;
	if (glst) {
		gl=glst;
	} else {
		gl = phash.Find(id);
	}
	GffObj* gh=NULL;
	if (gl) {
		for (int i=0;i<gl->Count();i++) {
			GffObj& gfo = *(gl->Get(i));
			if (ctg!=NULL && strcmp(ctg, gfo.getGSeqName())!=0)
				continue;
			if (strand && gfo.strand!='.' && strand != gfo.strand)
				continue;
			if (start>0) {
				if (abs((int)start-(int)gfo.start)> (int)GFF_MAX_LOCUS)
					continue;
				if (end>0 && (gfo.start>end || gfo.end<start))
					continue;
			}
			//must be the same transcript, according to given comparison criteria
			gh=&gfo;
			break;
		}
	}
	if (!glst) glst=gl;
	return gh;
}

GffObj* GffReader::updateParent(GffObj* newgfo, GffObj* parent) {
  //assert(parent);
  //assert(newgfo);
  parent->children.Add(newgfo);
  if (newgfo->parent==NULL) newgfo->parent=parent;
  newgfo->setLevel(parent->getLevel()+1);
  //if (parent->isGene()) {
  if (parent->gene_name!=NULL && newgfo->gene_name==NULL)
      newgfo->gene_name=Gstrdup(parent->gene_name);
  if (parent->geneID!=NULL && newgfo->geneID==NULL)
      newgfo->geneID=Gstrdup(parent->geneID);
  //}

  return newgfo;
}

GffObj* GffReader::newGffRec(GffLine* gffline, GffObj* parent, GffExon* pexon, GPVec<GffObj>* glst, bool replace_parent) {
  GffObj* newgfo=new GffObj(*this, *gffline);
  GffObj* r=NULL;
  gflst.Add(newgfo);
  //tag non-transcripts to be discarded later
  if (this->transcripts_Only && this->is_gff3 && gffline->ID!=NULL &&
		  gffline->exontype==exgffNone && !gffline->is_gene && !gffline->is_transcript) {
	  //unrecognized non-exon entity, should be discarded
	  newgfo->isDiscarded(true);
	  this->discarded_ids.Add(gffline->ID, new int(1));
  }
  if (replace_parent && glst) {
	r=gfoReplace(*glst, newgfo, parent);
	updateParent(r, parent);
  } else { //regular case of new GffObj creation
	  r=(glst) ? gfoAdd(*glst, newgfo) : gfoAdd(newgfo);
	  if (parent!=NULL) {
		updateParent(r, parent);
		if (pexon!=NULL) parent->removeExon(pexon);
	  }
  }
  return r;
}

GffObj* GffReader::newGffRec(BEDLine* bedline, GPVec<GffObj>* glst) {
  GffObj* newgfo=new GffObj(*this, *bedline);
  GffObj* r=NULL;
  gflst.Add(newgfo);
  r=(glst) ? gfoAdd(*glst, newgfo) : gfoAdd(newgfo);
  return r;
}

GffObj* GffReader::updateGffRec(GffObj* prevgfo, GffLine* gffline) {
 if (prevgfo==NULL) return NULL;
 //prevgfo->gffobj->createdByExon(false);
 if (gffline->ftype_id>=0)
	 prevgfo->ftype_id=gffline->ftype_id;
 else
	 prevgfo->ftype_id=prevgfo->names->feats.addName(gffline->ftype);
 prevgfo->start=gffline->fstart;
 prevgfo->end=gffline->fend;
 prevgfo->isGene(gffline->is_gene);
 prevgfo->isTranscript(gffline->is_transcript || gffline->exontype!=exgffNone);
 prevgfo->hasGffID(gffline->ID!=NULL);
 if (keep_Attrs) {
   if (prevgfo->attrs!=NULL) prevgfo->attrs->Clear();
   prevgfo->parseAttrs(prevgfo->attrs, gffline->info);
   }
 return prevgfo;
}


bool GffReader::readExonFeature(GffObj* prevgfo, GffLine* gffline, GHash<CNonExon>* pex) {
	//this should only be called before prevgfo->finalize()!
	bool r=true;
	if (gffline->strand!=prevgfo->strand) {
		if (prevgfo->strand=='.') {
			prevgfo->strand=gffline->strand;
		}
		else {
			GMessage("Error at %s (%c): exon %d-%d (%c) found on different strand; discarded.\n",
					prevgfo->gffID, prevgfo->strand, gffline->fstart, gffline->fend, gffline->strand,
					prevgfo->getGSeqName());
			return true;
		}
	}
	int gdist=(gffline->fstart>prevgfo->end) ? gffline->fstart-prevgfo->end :
			((gffline->fend<prevgfo->start)? prevgfo->start-gffline->fend :
					0 );
	if (gdist>(int)GFF_MAX_LOCUS) { //too far apart, most likely this is a duplicate ID
		GMessage("Error: duplicate GFF ID '%s' (or exons too far apart)!\n",prevgfo->gffID);
		//validation_errors = true;
		r=false;
		if (!gff_warns) exit(1);
	}
	int eidx=prevgfo->readExon(*this, *gffline);
	if (pex!=NULL && eidx>=0) {
		//if (eidx==0 && gffline->exontype>0) prevgfo->isTranscript(true);
		if (gffline->ID!=NULL && gffline->exontype==exgffNone)
		   subfPoolAdd(*pex, prevgfo);
	}
	return r;
}

CNonExon* GffReader::subfPoolCheck(GffLine* gffline, GHash<CNonExon>& pex, char*& subp_name) {
  CNonExon* subp=NULL;
  subp_name=NULL;
  for (int i=0;i<gffline->num_parents;i++) {
    if (transcripts_Only && discarded_ids.Find(gffline->parents[i])!=NULL)
        continue;
    subp_name=gfoBuildId(gffline->parents[i], gffline->gseqname); //e.g. mRNA name
    subp=pex.Find(subp_name);
    if (subp!=NULL)
       return subp;
    GFREE(subp_name);
    }
  return NULL;
}

void GffReader::subfPoolAdd(GHash<CNonExon>& pex, GffObj* newgfo) {
//this might become a parent feature later
if (newgfo->exons.Count()>0) {
   char* xbuf=gfoBuildId(gffline->ID, gffline->gseqname);
   pex.Add(xbuf, new CNonExon(newgfo, newgfo->exons[0], *gffline));
   GFREE(xbuf);
   }
}

GffObj* GffReader::promoteFeature(CNonExon* subp, char*& subp_name, GHash<CNonExon>& pex) {
  GffObj* prevp=subp->parent; //grandparent of gffline (e.g. gene)
  //if (prevp!=gflst[subp->idx])
  //  GError("Error promoting subfeature %s, gflst index mismatch?!\n", subp->gffline->ID);
  subp->gffline->discardParent();
  GffObj* gfoh=newGffRec(subp->gffline, prevp, subp->exon);
  pex.Remove(subp_name); //no longer a potential parent, moved it to phash already
  prevp->promotedChildren(true);
  return gfoh; //returns the holder of newly promoted feature
}

//In the rare cases where the GFF/GTF stream is properly formatted
// i.e. when all sub-features are grouped with (and preceded by) their parent!
GffObj* GffReader::readNext() { //user must free the returned GffObj*
 GffObj* gfo=NULL;
 //GSeg tseg(0,0); //transcript boundaries
 char* lastID=NULL;
 if (is_BED) {
	 if (nextBEDLine()) {
		 gfo=new GffObj(*this, *bedline);
		 //tseg.start=gfo->start;
		 //tseg.end=gfo->end;
		 delete bedline;
		 bedline=NULL;
	 } //else return NULL;
 } else { //GFF parsing
    while (nextGffLine()!=NULL) {
    	char* tid=gffline->ID;
    	if (gffline->is_exon) tid=gffline->parents[0];
    	else //not an exon
    		if (!(gffline->is_transcript || gffline->is_gene))
    			tid=NULL; //WARNING: only parsing transcript && gene records here
    	//if (tid==NULL || gffline->num_parents>1) {
    	if (tid==NULL) { //not a suitable transcript ID found, skip this line
    		delete gffline;
    		gffline=NULL;
    		continue;
    	}
    	bool sameID=(lastID!=NULL && strcmp(lastID, tid)==0);
    	if (sameID) {
    		if (gfo==NULL) GError("Error: same transcript ID but GffObj not initialized?!(%s)\n", tid);
    		//TODO: if gffline->is_transcript: trans-splicing!
    		if (!gffline->is_exon) {
    			GMessage("Warning: skipping unexpected non-exon record with previously seen ID:\n%s\n", gffline->dupline);
    			delete gffline;
    			gffline=NULL;
    			continue;
    		}
    		readExonFeature(gfo, gffline); //also takes care of adding CDS segments
    	} else { //new transcript
    		if (gfo==NULL) {
    			//start gathering this transcript's data now
    			gfo=new GffObj(*this, *gffline);
    			//GFREE(lastID);
    			lastID=Gstrdup(tid);
    			/*if (gffline->is_transcript) {
    				tseg.start=gffline->fstart;
    				tseg.end=gffline->fend;
    			}*/
    		}
    		else {
    			//this gffline is for the next transcript!
    			//return what we've got so far
    			//return gfo;
    			break;
    		}
    	} //transcript ID change
    	//gffline processed, move on
		delete gffline;
		gffline=NULL;
    } //while nextgffline()
 } //GFF records
 GFREE(lastID);
 //gfo populated with all its sub-features (or eof reached)
 if (gfo!=NULL) {
	gfo->finalize(this);
 }
 return gfo;
}

//Usually we have to parse the whole file because exons and other subfeatures can be scattered, unordered in the input
// (thanks to the annoyingly loose GFF3 specs)
//Trans-splicing and fusions shall only be accepted in proper GFF3 format, i.e. multiple transcript features
//with the same ID but NOT overlapping/continuous
//  *** BUT (exception): proximal xRNA features with the same ID, on the same strand, will be merged
//  and the segments will be treated like exons (e.g. TRNAR15 (rna1940) in RefSeq)
void GffReader::readAll() {
	bool validation_errors = false;
	if (is_BED) {
		while (nextBEDLine()) {
			GPVec<GffObj>* prevgflst=NULL;
			GffObj* prevseen=gfoFind(bedline->ID, prevgflst, bedline->gseqname, bedline->strand, bedline->fstart);
			if (prevseen) {
			//duplicate ID -- but this could also be a discontinuous feature according to GFF3 specs
			  //e.g. a trans-spliced transcript - but segments should not overlap
				if (prevseen->overlap(bedline->fstart, bedline->fend)) {
					//overlapping feature with same ID is going too far
					GMessage("Error: overlapping duplicate BED feature (ID=%s)\n", bedline->ID);
					//validation_errors = true;
					if (gff_warns) { //validation intent: just skip the feature, allow the user to see other errors
						delete bedline;
						bedline=NULL;
						continue;
					}
					else exit(1);
				}
				//create a separate entry (true discontinuous feature?)
				prevseen=newGffRec(bedline, prevgflst);
				if (gff_warns) {
					GMessage("Warning: duplicate BED feature ID %s (%d-%d) (discontinuous feature?)\n",
							bedline->ID, bedline->fstart, bedline->fend);
				}
			}
			else {
				newGffRec(bedline, prevgflst);
			}
			delete bedline;
			bedline=NULL;
		}
	}
	else { //regular GFF/GTF or perhaps TLF?
		//loc_debug=false;
		GHash<CNonExon> pex; //keep track of any parented (i.e. exon-like) features that have an ID
		//and thus could become promoted to parent features
		while (nextGffLine()!=NULL) {
			GffObj* prevseen=NULL;
			GPVec<GffObj>* prevgflst=NULL;
			if (gffline->ID && gffline->exontype==exgffNone) {
				//parent-like feature ID (mRNA, gene, etc.) not recognized as an exon feature
				//check if this ID was previously seen on the same chromosome/strand within GFF_MAX_LOCUS distance
				prevseen=gfoFind(gffline->ID, prevgflst, gffline->gseqname, gffline->strand, gffline->fstart);
				if (prevseen) {
					//same ID seen in the same locus/region
					if (prevseen->createdByExon()) {
						if (gff_warns && (prevseen->start<gffline->fstart ||
								prevseen->end>gffline->fend))
							GMessage("Warning: invalid coordinates for %s parent feature (ID=%s)\n", gffline->ftype, gffline->ID);
						//an exon of this ID was given before
						//this line has the main attributes for this ID
						updateGffRec(prevseen, gffline);
					}
					else { //possibly a duplicate ID -- but this could also be a discontinuous feature according to GFF3 specs
					    //e.g. a trans-spliced transcript - though segments should not overlap!
						bool gtf_gene_dupID=(prevseen->isGene() && gffline->is_gtf_transcript);
						if (prevseen->overlap(gffline->fstart, gffline->fend) && !gtf_gene_dupID) {
							//in some GTFs a gene ID may actually be the same with the parented transcript ID (thanks)
							//overlapping feature with same ID is going too far
							GMessage("Error: discarding overlapping duplicate %s feature (%d-%d) with ID=%s\n", gffline->ftype,
									gffline->fstart, gffline->fend, gffline->ID);
							//validation_errors = true;
							if (gff_warns) { //validation intent: just skip the feature, allow the user to see other errors
								delete gffline;
								gffline=NULL;
								continue;
							}
							//else exit(1);
						}
						if (gtf_gene_dupID) {
							//special GTF case where parent gene_id matches transcript_id (sigh)
							prevseen=newGffRec(gffline, prevseen, NULL, prevgflst, true);
						}
						else {
							//create a separate entry (true discontinuous feature)
							prevseen=newGffRec(gffline, prevseen->parent, NULL, prevgflst);
							if (gff_warns) {
								GMessage("Warning: duplicate feature ID %s (%d-%d) (discontinuous feature?)\n",
										gffline->ID, gffline->fstart, gffline->fend);
							}
						}
					} //duplicate ID in the same locus
				} //ID seen previously in the same locus
			} //parent-like ID feature (non-exon)
			if (gffline->parents==NULL) {
				//top level feature (transcript, gene), no parents (or parents can be ignored)
				if (!prevseen) newGffRec(gffline, NULL, NULL, prevgflst);
			}
			else { //--- it's a child feature (exon/CDS or even a mRNA with a gene as parent)
				//updates all the declared parents with this child
				bool found_parent=false;
				if (gffline->is_gtf_transcript && prevseen && prevseen->parent) {
					found_parent=true; //parent already found in special GTF case
				}
				else {
					GffObj* newgfo=prevseen;
					GPVec<GffObj>* newgflst=NULL;
					GVec<int> kparents; //kept parents (non-discarded)
					GVec< GPVec<GffObj>* > kgflst(false);
					GPVec<GffObj>* gflst0=NULL;
					for (int i=0;i<gffline->num_parents;i++) {
						newgflst=NULL;
						//if (transcriptsOnly && (
						if (discarded_ids.Find(gffline->parents[i])!=NULL) continue;
						if (!pFind(gffline->parents[i], newgflst))
							continue; //skipping discarded parent feature
						kparents.Add(i);
						if (i==0) gflst0=newgflst;
						kgflst.Add(newgflst);
					}
					if (gffline->num_parents>0 && kparents.Count()==0) {
						kparents.cAdd(0);
						kgflst.Add(gflst0);
					}
					for (int k=0;k<kparents.Count();k++) {
						int i=kparents[k];
						newgflst=kgflst[k];
						GffObj* parentgfo=NULL;
						if (gffline->is_transcript || gffline->exontype==exgffNone) {//likely a transcript
							//parentgfo=gfoFind(gffline->parents[i], newgflst, gffline->gseqname,
							//		gffline->strand, gffline->fstart, gffline->fend);
							if (newgflst!=NULL && newgflst->Count()>0)
								parentgfo = newgflst->Get(0);
						}
						else {
							//for exon-like entities we only need a parent to be in locus distance,
							//on the same strand
							parentgfo=gfoFind(gffline->parents[i], newgflst, gffline->gseqname,
									gffline->strand, gffline->fstart);
						}
						if (parentgfo!=NULL) { //parent GffObj parsed earlier
							found_parent=true;
							if ((parentgfo->isGene() || parentgfo->isTranscript()) && (gffline->is_transcript ||
									 gffline->exontype==exgffNone)) {
								//not an exon, but could be a transcript parented by a gene
								// *or* by another transcript (! miRNA -> primary_transcript)
								if (newgfo) {
									updateParent(newgfo, parentgfo);
								}
								else {
									newgfo=newGffRec(gffline, parentgfo);
								}
							}
							else { //potential exon subfeature?
								bool addingExon=false;
								if (transcripts_Only) {
									if (gffline->exontype>0) addingExon=true;
								}
								else { //always discard silly "intron" features
									if (! (gffline->exontype==exgffIntron && (parentgfo->isTranscript() || parentgfo->exons.Count()>0)))
									  addingExon=true;
								}
								if (addingExon)
									if (!readExonFeature(parentgfo, gffline, &pex))
									   validation_errors=true;

							}
						} //overlapping parent feature found
					} //for each parsed parent Id
					if (!found_parent) { //new GTF-like record starting directly here as a subfeature
						//or it could be some chado GFF3 barf with exons coming BEFORE their parent :(
						//or it could also be a stray transcript without a parent gene defined previously
						//check if this feature isn't parented by a previously stored "child" subfeature
						char* subp_name=NULL;
						CNonExon* subp=NULL;
						if (!gffline->is_transcript) { //don't bother with this check for obvious transcripts
							if (pex.Count()>0) subp=subfPoolCheck(gffline, pex, subp_name);
							if (subp!=NULL) { //found a subfeature that is the parent of this (!)
								//promote that subfeature to a full GffObj
								GffObj* gfoh=promoteFeature(subp, subp_name, pex);
								//add current gffline as an exon of the newly promoted subfeature
								if (!readExonFeature(gfoh, gffline, &pex))
									validation_errors=true;
							}
						}
						if (subp==NULL) { //no parent subfeature seen before
							//loc_debug=true;
							GffObj* ngfo=prevseen;
							if (ngfo==NULL) {
								//if it's an exon type, create directly the parent with this exon
								//but if it's recognized as a transcript, the object itself is created
								ngfo=newGffRec(gffline, NULL, NULL, newgflst);
							}
							if (!ngfo->isTranscript() &&
									gffline->ID!=NULL && gffline->exontype==0)
								subfPoolAdd(pex, ngfo);
							//even those with errors will be added here!
						}
						GFREE(subp_name);
					} //no previous parent found
				}
			} //parented feature
			//--
			delete gffline;
			gffline=NULL;
		}//while gff lines
	}
	if (gflst.Count()>0) {
		gflst.finalize(this); //force sorting by locus if so constructed
	}
	// all gff records are now loaded in GList gflst
	// so we can free the hash
	phash.Clear();
	//tids.Clear();
	if (validation_errors) {
		exit(1);
	}
}

void GfList::finalize(GffReader* gfr) { //if set, enforce sort by locus
  GList<GffObj> discarded(false,true,false);
  for (int i=0;i<Count();i++) {
    //finalize the parsing of each GffObj
    fList[i]->finalize(gfr);
    if (fList[i]->isDiscarded()) {
       discarded.Add(fList[i]);
       //inform parent that thiis child is removed
       if (fList[i]->parent!=NULL) {
    	   GPVec<GffObj>& pchildren=fList[i]->parent->children;
    	   for (int c=0;c<pchildren.Count();c++) {
    		   if (pchildren[c]==fList[i]) {
    			   pchildren.Delete(c);
    			   break;
    		   }
    	   }
       }
       if (fList[i]->children.Count()>0) { //inform children that the parent was removed
      	 for (int c=0;c<fList[i]->children.Count();c++) {
      		 fList[i]->children[c]->parent=NULL;
      		 if (gfr->keep_Attrs)
      			 //inherit the attributes of discarded parent (e.g. pseudo=true; )
      			 fList[i]->children[c]->copyAttrs(fList[i]);
      	 }
       }
       this->Forget(i);
    }
  }
  if (discarded.Count()>0) {
          this->Pack();
  }
  if (gfr->sortByLoc) {
    this->setSorted(false);
    if (gfr->refAlphaSort)
      this->setSorted((GCompareProc*)gfo_cmpByLoc);
    else
      this->setSorted((GCompareProc*)gfo_cmpRefByID);
  }
}

bool GffObj::reduceExonAttrs(GList<GffExon>& segs) {
	bool attrs_discarded=false;
	for (int a=0;a<segs[0]->attrs->Count();a++) {
		int attr_id=segs[0]->attrs->Get(a)->attr_id;
		char* attr_name=names->attrs.getName(attr_id);
		char* attr_val =segs[0]->attrs->Get(a)->attr_val;
		bool sameExonAttr=true;
		bool discardAll=(GstrEq("exon_id", attr_name) || GstrEq("exon_number", attr_name));
		if (!discardAll)
			for (int i=1;i<segs.Count();i++) {
				char* ov=segs[i]->getAttr(attr_id);
				if (ov==NULL || (strcmp(ov,attr_val)!=0)) {
					sameExonAttr=false;
					break;
				}
			}
		if (sameExonAttr) {
			//delete this attribute from exon level
			attrs_discarded=true;
			if (!discardAll) {
				//add the attribute to transcript level
				//rename it if it exists and is different for the transcript!
				char* t_val=NULL;
				bool same_aval=false;
				if (this->attrs!=NULL &&
						(t_val=this->attrs->getAttr(attr_id))!=NULL) {
					//same attribute name already exists for the transcript!
					//write it using CDS_ or exon_ prefix
					same_aval=(strcmp(attr_val, t_val)==0);
					if (!same_aval) {
						//add renamed attribute
						const char* prefix = (&segs==cdss) ? "CDS_" : "exon_";
						char* new_attr_name=NULL;
						GMALLOC(new_attr_name, strlen(prefix)+strlen(attr_name)+1);
						new_attr_name[0]=0;
						strcat(new_attr_name, prefix);
						strcat(new_attr_name, attr_name);
						this->attrs->add_or_update(names, new_attr_name, attr_val);
						GFREE(new_attr_name);
					}
				}
				else { //no such attribute exists for the transcript, copy it from the exon
						this->addAttr(attr_name, attr_val);
				}
			}
			for (int i=1;i<segs.Count();i++) {
				removeExonAttr(*(segs[i]), attr_id);
			}
			segs[0]->attrs->freeItem(a);
		} //sameExonAttr
	}
	if (attrs_discarded) segs[0]->attrs->Pack();
	return attrs_discarded;
}
//return the segs index of segment containing coord:
int GffObj::whichExon(uint coord, GList<GffExon>* segs) {
	 //segs MUST be sorted by GSeg order (start coordinate)
	if (segs==NULL) segs=&exons;
	if (segs->Count()==0) return -1;
	if (coord<segs->First()->start || coord>segs->Last()->end)
		return -1;
	if (segs->Count()<6) {
		//simple scan
		for (int i=0;i<segs->Count();i++)
			if ((*segs)[i]->overlap(coord))
				return i;
		return -1;
	}
	else { //use quick search
		int i=0;
		int l=0; //lower boundary
		int h=segs->Count()-1; //higher boundary
		while (l<=h) {
			i = (l+h) >> 1; //range midpoint
			if (coord > segs->Get(i)->end)
				l=i+1;
			else { //coord <= segs->Get(i)->end
				if (coord >= segs->Get(i)->start) {
					return i;
				}
				//here: coord < segs->Get(i)->start
				h = i-1;
			}
		}
	}
	return -1;
}

bool GffObj::processGeneSegments(GffReader* gfr) {
	/* procedure:
	 1)store the info about any X_gene_segment entries in a GVec<int>
	     (just storing their index in gene->children[] list)
     2)for each CDS, group them by ID in GHash<GeneCDSChain> (and a GPVec<GeneCDSChain> for storage)
     3)for each GeneCDSChain, collect _gene_segments having a containment-relationship and rank them by lowest noncov
     4)for each GeneCDSChain, pick best _gene_segment match (if any) and transfer CDSs to it
	*/
    GVec<int> geneSegs; //X_gene_segment features (children transcripts of this gene)
    GHash<GeneCDSChain> cdsChainById(false); // hash of CDS chains: CDS feature grouped by ID
    GPVec<GeneCDSChain> cdsChains; // CDS chains storage
	if (cdss==NULL || cdss->Count()==0 || children.Count()==0)
		return false; //we shouldn't be here
	//check if we have any _gene_segment children for this gene
    for (int i=0;i<children.Count();i++)
    	if (children[i]->flag_GENE_SEGMENT) {
    		if (children[i]->hasCDS() || children[i]->cdss!=NULL) {
    			GMessage("Warning: will not transfer CDS from %s to gene_segment %s which already has its own\n",
    					gffID, children[i]->gffID);
    			continue;
    		}
    		geneSegs.Add(i);
    	}
    if (geneSegs.Count()==0) {
       if (gfr->gff_warns)
    	   GMessage("Warning: gene %s has CDS and transcripts but no suitable _gene_segment features\n",gffID);
       return false; //nothing to do
    }
    //group CDSs into CDS chains by their ID:
    for (int i=0;i<cdss->Count();i++) {
    	char* id=(char*)(cdss->Get(i)->uptr);
    	if (id==NULL) continue; //should never happen
    	GeneCDSChain *gcc=cdsChainById.Find(id);
    	if (gcc!=NULL)
             gcc->addCDS(i, cdss->Get(i)->start, cdss->Get(i)->end);
    	else { //new CDS chain:
    	     gcc=new GeneCDSChain(i, cdss->Get(i)->start, cdss->Get(i)->end);
    	     cdsChains.Add(gcc);
    	     cdsChainById.shkAdd(id, gcc);
    	}
    }
    for (int i=0;i<cdss->Count();i++) {
 	   GFREE(cdss->Get(i)->uptr); //no CDS ID no longer needed
    }

    //collect _gene_segment containers for each CDS chain
    int cds_moved=0;
    for (int i=0;i<cdsChains.Count();i++) {
    	GeneCDSChain &gc=*(cdsChains[i]);
        for (int si=0;si<geneSegs.Count();si++) {
        	GffObj& t=*(children[geneSegs[si]]);
        	//if (t.hasCDS() || t.cdss!=NULL) continue; //already transferred
        	int novl=-1;
        	if (gc.containedBy(t, novl))
        		gc.addMatch(geneSegs[si], novl, si);
        }
        //-- assess the collected matches for this CDS chain:
        if (gc.mxs.Count()==0) {
    		GMessage("Warning: could not find the corresponding gene segment for CDS chain %d-%d of gene %s\n",
    		    				gc.start, gc.end, gffID);
        	continue;
        }
        GffObj* t=children[gc.mxs.First().child_idx];
		for (int c=0;c<gc.cdsList.Count();c++) {
			t->transferCDS(cdss->Get(gc.cdsList[c].idx));
			cdss->Forget(gc.cdsList[c].idx);
			cds_moved++;
		}
		// also remove it from the list of gene_segments to be mapped
		geneSegs.Delete(gc.mxs.First().gsegidx); //assigned, should no longer be checked against other CDS chains
		if (t->isFinalized()) t->finalize(gfr);

    }
    if (cds_moved>0) cdss->Pack();
    if (cdss->Count()==0) {
    	delete cdss;
    	cdss=NULL;
    	if (exons.Count()==0) isTranscript(false);
    }
	return true;
}

GffObj* GffObj::finalize(GffReader* gfr) {
	if (this->createdByExon() && this->end-this->start<10 && this->exons.Count()<=1) {
		//? misleading exon-like feature parented by an exon or CDS mistakenly
		//  interpreted as a standalone transcript
		// example: GENCODE gff3 feature "stop_codon_redefined_as_selenocysteine" which is
		// parented by a CDS !
		if (cdss==NULL || cdss->Count()<=1) {
			if (gfr->showWarnings()) {
			 GMessage("Warning: discarding suspicious '%s' record (ID=%s)\n",this->getFeatureName(),gffID);
			}
		   isDiscarded(true);
		}
	}
	if (!isDiscarded()) {
		bool noExons=(exons.Count()==0 && (cdss==NULL || cdss->Count()==0));
		if (noExons) {
			if (isTranscript() || (isGene() && children.Count()==0 && gfr->gene2exon)) {
				//add exon feature to an exonless transcript/gene
				addExon(this->start, this->end, exgffExon);
				//effectively this becomes a transcript (even childless genes if gene2exon)
				isTranscript(true);
			}
		}
		else { //it has exons or CDSs
			if (cdss!=NULL && isGene() && children.Count()>0) {
				//check for X_gene_segment processing
				processGeneSegments(gfr);//distribute the cdss to children _gene_segments
			} // _gene_segment processing
		}
	}
	if (cdss!=NULL && isGene()) //in case we stored IDs for gene_segment features
		for (int i=0;i<cdss->Count();i++) {
			GFREE(cdss->Get(i)->uptr);
		}
	if (gfr->transcripts_Only && !isTranscript() &&
			!(gfr->keep_Genes && isGene())) {
		//discard non-transcripts, unless it's a gene and keepGenes was specified
		isDiscarded(true);
	}
	isFinalized(true);

	if (isDiscarded()) {
		//just in case we have cds with uptr in use (X_gene_segment), free them
		uptr=NULL;
		udata=0;
		return this;
	}
	if (isTranscript()) {
		isCDSOnly(cdss!=NULL && exons.Count()==0 && cdss->Count()>0);
		subftype_id=isCDSOnly() ? gff_fid_CDS : gff_fid_exon;
	}
	if (cdss!=NULL && cdss->Count()>0) {
		CDstart=cdss->First()->start;
		CDend=cdss->Last()->end;
		CDphase=(strand=='-')? cdss->Last()->phase : cdss->First()->phase;
		bool updatePhase=(CDphase=='.' || CDphase==0);
		if (!updatePhase)
			for (int i=0;i<cdss->Count();++i)
				if ((*cdss)[i]->phase<'0') {
					updatePhase=true;
					break;
				}
		if (updatePhase) updateCDSPhase(*cdss);
		//there are GFFs out there which only provide UTR and CDS records instead of full exons
		//so make sure we add all CDS segments to exons, if they are not already there
		for (int i=0;i<cdss->Count();++i) {
			int eidx=addExon((*cdss)[i]->start, (*cdss)[i]->end, exgffExon, 0, (*cdss)[i]->score);
			if (eidx<0) GError("Error: could not reconcile CDS %d-%d with exons of transcript %s\n",
					(*cdss)[i]->start, (*cdss)[i]->end, gffID);
		}
	}
	else if (CDstart==0) {//no CDS, no phase
		CDphase=0;
		CDend=0;
	}
	//-- attribute reduction for some records which
	//   repeat the exact same attr=value for every exon
	bool reduceAttributes=(gfr->keep_Attrs && !gfr->noExonAttrs &&
			!gfr->keep_AllExonAttrs && exons.Count()>0 && exons[0]->attrs!=NULL);
	if (reduceAttributes) {
		//for each attribute of the 1st exon, if it has the
		//same value for all other exons, move it to transcript level
		//bool reduced=reduceExonAttrs(exons);
		reduceExonAttrs(exons);
		//if (gfr->showWarnings() && reduced)
		//	GMessage("Info: duplicate exon attributes reduced for %s\n", gffID);
		//do the same for CDS segments, if any
		if (cdss!=NULL && cdss->Count()>0 && (*cdss)[0]->attrs!=NULL) {
			//reduced=
			reduceExonAttrs(*cdss);
			//if (gfr->showWarnings() && reduced)
			//	GMessage("Info: duplicate CDS attributes reduced for %s\n", gffID);
		}
	}
	//merge close exons if requested
	if (exons.Count()>0 && isTranscript()) {
		if (gfr->merge_CloseExons) {
			for (int i=0;i<exons.Count()-1;i++) {
				int ni=i+1;
				uint mend=exons[i]->end;
				while (ni<exons.Count()) {
					int dist=(int)(exons[ni]->start-mend-1); //<0 = overlap, 0 = adjacent, >0 = bases apart
					if (dist>GFF_MIN_INTRON) break; //no merging with next segment
					if (gfr!=NULL && gfr->gff_warns && dist!=0 && (exons[ni]->exontype!=exgffUTR && exons[i]->exontype!=exgffUTR)) {
						GMessage("Warning: merging adjacent/overlapping segments (distance=%d) of %s on %s (%d-%d, %d-%d)\n",
								dist, gffID, getGSeqName(), exons[i]->start, exons[i]->end,exons[ni]->start, exons[ni]->end);
					}
					mend=exons[ni]->end;
					exons[i]->end=mend;
					if (exons[ni]->attrs!=NULL && (exons[i]->attrs==NULL ||
							exons[i]->attrs->Count()<exons[ni]->attrs->Count())) {
						//use the other exon attributes, if it has more
						delete(exons[i]->attrs);
						exons[i]->attrs=exons[ni]->attrs;
						exons[ni]->attrs=NULL;
					}
					exons.Delete(ni);
				} //check for merge with next exon
			} //for each exon
		} //merge close exons
		if (isCDSOnly() && exons.Count()!=cdss->Count())
			isCDSOnly(false);
	}
	//-- check features vs their exons' span
	if (isTranscript()) {
	   if (exons.Count()>0) {
		 if (gfr->gff_warns && (this->start!=exons.First()->start ||
				 this->end!=exons.Last()->end) )
			 GMessage("Warning: adjusted transcript %s boundaries according to terminal exons.\n",
					 gffID);
	     this->start=exons.First()->start;
	     this->end=exons.Last()->end;
	   }
	}
	else { //non-transcripts just have to be at least as wide as their sub-features
	  if (exons.Count()>0) {
		  bool adj=false;
		  if (this->start>exons.First()->start) {
			  this->start=exons.First()->start;
			  adj=true;
		  }
		  if (this->end<exons.Last()->end) {
			  this->end=exons.First()->end;
			  adj=true;
		  }
		  if (gfr->gff_warns && adj)
			  GMessage("Warning: adjusted %s %s boundaries according to terminal sub-features.\n",
			  					 this->getFeatureName(), gffID);
	  }
	}
	//-- update covlen
	covlen=0;
	for (int i=0;i<exons.Count();++i) covlen+=exons[i]->len();
	//-- check if CDS segments are different from exons and thus worth keeping separately in cdss
	if (cdss!=NULL && cdss->Count()>0) {
		bool cds_exComp=true; //CDSs are exon-compatible (no need to keep them separately)
		if (cdss->Count()==1) {
			//check that the CDS segment is within a single exon
			int start_eidx=-1;
			int end_eidx=-1;
			for (int i=0;i<exons.Count();i++) {
				//GMessage("[DBG:] checking if CDS %d-%d is within exon %d-%d\n", CDstart, CDend, exons[i]->start,
				//		exons[i]->end);
				if (CDstart>=exons[i]->start && CDstart<=exons[i]->end) {
					start_eidx=i;
				}
				if (CDend>=exons[i]->start || CDend<=exons[i]->end ) {
					end_eidx=i;
				}
				if (start_eidx>=0 && end_eidx>=0) break;
			}
			cds_exComp=(start_eidx==end_eidx && start_eidx>=0);
			if (!cds_exComp) GMessage("Warning: transcript %s has incorrect CDS segment definition (%d-%d)!\n",
					gffID, CDstart, CDend);
			cds_exComp=true; //just to free cdss, even though it's wrong
		} else {
			if (cdss->Count()>exons.Count()) {
				cds_exComp=false;
			} else { //2 or more CDS segments
				//CDSs should be intron compatible with exons, and CDS ends should be within exons
				int imax=exons.Count()-1;
				int jmax=cdss->Count()-1;
				int i=0;
				int j=0;
				//find which exon has CDstart
				for (i=0;i<=imax;++i)
					if (CDstart>=exons[i]->start
							&& CDstart<=exons[i]->end) break;
				if (i>imax) cds_exComp=false;
				else { //check the introns now
					while (i<imax && j<jmax) {
						if (exons[i]->end!=(*cdss)[j]->end ||
								exons[i+1]->start!=(*cdss)[j+1]->start) {
							cds_exComp=false;
							break;
						}
						++i;
						++j;
					}
					//now j must be the last segment of cdss and CDend must be within exon[i]
					if (cds_exComp)
						if (j!=jmax || CDend>exons[i]->end || CDend<exons[i]->start)
							cds_exComp=false;
				}
			}
		} //multiple CDS segments
		if (cds_exComp) {
			if (isCDSOnly() && cdss->Count()==exons.Count())
				for (int i=0;i<cdss->Count();i++)
					exons[i]->phase=cdss->Get(i)->phase;
			if (gfr->keep_Attrs && !gfr->noExonAttrs) {
				int eidx=whichExon((*cdss)[0]->start, &exons);
				if (eidx<0)
					GError("Error finding CDS coordinate inside exons (?) for %s\n",
						    gffID);
				for (int i=0;i<cdss->Count();i++) {
					if (isCDSOnly()) //eidx should be the same with i
						exons[eidx]->phase=cdss->Get(i)->phase;
					if ((*cdss)[i]->attrs!=NULL && (*cdss)[i]->attrs->Count()>0) {
						if (exons[eidx]->attrs==NULL)
							exons[eidx]->attrs=new GffAttrs();
						exons[eidx]->attrs->copyAttrs((*cdss)[i]->attrs, true);
						if (exons[eidx]->attrs->Count()==0) {
							delete exons[eidx]->attrs;
							exons[eidx]->attrs=NULL;
						}
					}
					++eidx;
				}
			}
			delete cdss;
			cdss=NULL;
			//this->isXCDS(false);
		} else this->isXCDS(true);
	}//cdss check

	//--- collect stats for the reference genomic sequence
	if (gfr->gseqtable.Count()<=gseq_id) {
		gfr->gseqtable.setCount(gseq_id+1);
	}
	GSeqStat* gsd=gfr->gseqtable[gseq_id];
	if (gsd==NULL) {
		gsd=new GSeqStat(gseq_id,names->gseqs.getName(gseq_id));
		//gfr->gseqtable.Put(gseq_id, gsd);
		gfr->gseqtable[gseq_id]=gsd;
		gfr->gseqStats.Add(gsd);
	}
	gsd->fcount++;
	if (start<gsd->mincoord) gsd->mincoord=start;
	if (end>gsd->maxcoord) gsd->maxcoord=end;
	if (this->len()>gsd->maxfeat_len) {
		gsd->maxfeat_len=this->len();
		gsd->maxfeat=this;
	}
	uptr=NULL;
	udata=0;
	return this;
}

void GffObj::printExonList(FILE* fout) {
	//print comma delimited list of exon intervals
	for (int i=0;i<exons.Count();++i) {
		if (i>0) fprintf(fout, ",");
		fprintf(fout, "%d-%d",exons[i]->start, exons[i]->end);
	}
}

void GffObj::printCDSList(FILE* fout) {
	//print comma delimited list of CDS intervals
	if (!hasCDS()) return;
	GVec<GffExon> cds;
	this->getCDSegs(cds); //also uses/prepares the CDS phase for each CDS segment
	for (int i=0;i<cds.Count();i++) {
		if (i>0) fprintf(fout, ",");
		fprintf(fout, "%d-%d", cds[i].start, cds[i].end);
	}
}

void BED_addAttribute(FILE* fout, int& acc, const char* format,... ) {
	++acc;
	if (acc==1) fprintf(fout, "\t");
	       else fprintf(fout, ";");
    va_list arguments;
    va_start(arguments,format);
    vfprintf(fout,format,arguments);
    va_end(arguments);
}

void GffObj::printBED(FILE* fout, bool cvtChars, char* dbuf, int dbuf_len) {
//print a BED-12 line + GFF3 attributes in 13th field
 int cd_start=CDstart>0? CDstart-1 : start-1;
 int cd_end=CDend>0 ? CDend : end;
 char cdphase=(CDphase>0) ? CDphase : '0';
 fprintf(fout, "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t%c,0,0", getGSeqName(), start-1, end, getID(),
		 100, strand, cd_start, cd_end, cdphase);
 if (exons.Count()>0) {
	 int i;
	 fprintf(fout, "\t%d\t", exons.Count());
	 for (i=0;i<exons.Count();++i)
		 fprintf(fout,"%d,",exons[i]->len());
	 fprintf(fout, "\t");
	 for (i=0;i<exons.Count();++i)
		 fprintf(fout,"%d,",exons[i]->start-start);
 } else { //no-exon feature(!), shouldn't happen
	 fprintf(fout, "\t1\t%d,\t0,", len());
 }
 //now add the GFF3 attributes for in the 13th field
 int numattrs=0;
 if (CDstart>0) BED_addAttribute(fout, numattrs,"CDS=%d:%d",CDstart-1, CDend);
 if (CDphase>0) BED_addAttribute(fout, numattrs,"CDSphase=%c", CDphase);
 if (geneID!=NULL)
	 BED_addAttribute(fout, numattrs, "geneID=%s",geneID);
 if (gene_name!=NULL)
    fprintf(fout, ";gene_name=%s",gene_name);
 if (attrs!=NULL) {
    for (int i=0;i<attrs->Count();i++) {
      const char* attrname=names->attrs.getName(attrs->Get(i)->attr_id);
      const char* attrval=attrs->Get(i)->attr_val;
      if (attrval==NULL || attrval[0]=='\0') {
    	  BED_addAttribute(fout, numattrs,"%s",attrname);
    	  continue;
      }
      if (cvtChars) {
    	  decodeHexChars(dbuf, attrval, dbuf_len-1);
    	  BED_addAttribute(fout, numattrs, "%s=%s", attrname, dbuf);
      }
      else
    	  BED_addAttribute(fout, numattrs,"%s=%s", attrname, attrs->Get(i)->attr_val);
    }
 }
 fprintf(fout, "\n");
}

void GffObj::parseAttrs(GffAttrs*& atrlist, char* info, bool isExon, bool CDSsrc) {
  if (names==NULL)
     GError(ERR_NULL_GFNAMES, "parseAttrs()");
  if (atrlist==NULL) {
      atrlist=new GffAttrs();
  }
  bool exon2transcript=(isExon && atrlist==this->attrs);
  char* endinfo=info+strlen(info);
  char* start=info;
  char* pch=start;
  while (start<endinfo) {
    while (*start==' ' && start<endinfo) start++;
    pch=strchr(start, ';');
    if (pch==NULL) pch=endinfo;
       else {  *pch='\0'; pch++; }
    char* ech=strchr(start,'=');
    if (ech!=NULL) { // attr=value format found
       *ech='\0';
       if (exon2transcript) { //we do NOT want these exon attributes at transcript level
          if (startsiWith(start, "exon_") || strcmp(start, "exon")==0)  {
        	  start=pch; continue;
          }
       }
       ech++;
       while (*ech==' ' && ech<endinfo) ech++;//skip extra spaces after the '='
       /*
       if (isExon && startsiWith(start, "protein")) {
         //--protein info should never be left at exon level
         // ( attribute policing much? :p )
         this->addAttr(start, ech);
         start=pch;
         continue;
       }
       */
       if (exon2transcript)
            atrlist->add_if_new(this->names, start, ech); //never override transcript attribute with exon's
       else atrlist->add_or_update(this->names, start, ech, CDSsrc); //overwrite previous attr with the same name
    }
    start=pch;
  } //while info characters
  if (atrlist->Count()==0) { delete atrlist; atrlist=NULL; }
}

void GffObj::addAttr(const char* attrname, const char* attrvalue) {
  if (this->attrs==NULL)
      this->attrs=new GffAttrs();
  //this->attrs->Add(new GffAttr(names->attrs.addName(attrname),attrvalue));
  this->attrs->add_or_update(names, attrname, attrvalue);
}

void GffObj::copyAttrs(GffObj* from) { //typically from is the parent gene, and this is a transcript
	if (from==NULL || from->attrs==NULL || from->attrs->Count()==0) return;
	if (this->attrs==NULL) {
		this->attrs=new GffAttrs();
	}
	//special RefSeq case
	int desc_attr_id=names->attrs.getId("description"); //from gene
	int prod_attr_id=names->attrs.getId("product"); //from transcript (this)
	char* prod = (prod_attr_id>=0) ? this->attrs->getAttr(prod_attr_id) : NULL;

	for (int i=0;i<from->attrs->Count();++i) {
		//this->attrs->add_no_update(names, from->attrs->Get(i)->attr_id, from->attrs->Get(i)->attr_val);
		int aid=from->attrs->Get(i)->attr_id;
		//special case for GenBank refseq genes vs transcripts:
		if (prod && aid==desc_attr_id && strcmp(from->attrs->getAttr(desc_attr_id), prod)==0)
			continue; //skip description if product already there and the same
		bool haveit=false;
		for (int ai=0;ai<this->attrs->Count();++ai) {
			//do we have it already?
			if (aid==this->attrs->Get(ai)->attr_id) {
				haveit=true;
				break; //skip this, don't replace
			}
		}
		if (!haveit)
			this->attrs->Add(new GffAttr(aid, from->attrs->Get(i)->attr_val));
	}
}

void GffObj::setFeatureName(const char* feature) {
 //change the feature name/type for a transcript
 int fid=names->feats.addName(feature);
 if (monoFeature() && exons.Count()>0)
    this->subftype_id=fid;
 this->ftype_id=fid;
}

void GffObj::setRefName(const char* newname) {
 //change the feature name/type for a transcript
 int rid=names->gseqs.addName(newname);
 this->gseq_id=rid;
}



int GffObj::removeAttr(const char* attrname, const char* attrval) {
  if (this->attrs==NULL || attrname==NULL || attrname[0]==0) return 0;
  int aid=this->names->attrs.getId(attrname);
  if (aid<0) return 0;
  int delcount=0;  //could be more than one ?
  for (int i=0;i<this->attrs->Count();i++) {
     if (aid==this->attrs->Get(i)->attr_id) {
       if (attrval==NULL ||
          strcmp(attrval, this->attrs->Get(i)->attr_val)==0) {
             delcount++;
             this->attrs->freeItem(i);
             }
       }
     }
  if (delcount>0) this->attrs->Pack();
  return delcount;
}

int GffObj::removeAttr(int aid, const char* attrval) {
  if (this->attrs==NULL || aid<0) return 0;
  int delcount=0;  //could be more than one ?
  for (int i=0;i<this->attrs->Count();i++) {
     if (aid==this->attrs->Get(i)->attr_id) {
       if (attrval==NULL ||
          strcmp(attrval, this->attrs->Get(i)->attr_val)==0) {
             delcount++;
             this->attrs->freeItem(i);
             }
       }
     }
  if (delcount>0) this->attrs->Pack();
  return delcount;
}


int GffObj::removeExonAttr(GffExon& exon, const char* attrname, const char* attrval) {
  if (exon.attrs==NULL || attrname==NULL || attrname[0]==0) return 0;
  int aid=this->names->attrs.getId(attrname);
  if (aid<0) return 0;
  int delcount=0;  //could be more than one
  for (int i=0;i<exon.attrs->Count();i++) {
     if (aid==exon.attrs->Get(i)->attr_id) {
       if (attrval==NULL ||
          strcmp(attrval, exon.attrs->Get(i)->attr_val)==0) {
             delcount++;
             exon.attrs->freeItem(i);
             }
       }
     }
  if (delcount>0) exon.attrs->Pack();
  return delcount;
}

int GffObj::removeExonAttr(GffExon& exon, int aid, const char* attrval) {
  if (exon.attrs==NULL || aid<0) return 0;
  int delcount=0;  //could be more than one
  for (int i=0;i<exon.attrs->Count();i++) {
     if (aid==exon.attrs->Get(i)->attr_id) {
       if (attrval==NULL ||
          strcmp(attrval, exon.attrs->Get(i)->attr_val)==0) {
             delcount++;
             exon.attrs->freeItem(i);
             }
       }
     }
  if (delcount>0) exon.attrs->Pack();
  return delcount;
}


char* GffObj::getUnspliced(GFaSeqGet* faseq, int* rlen, GMapSegments* seglst) {

    if (faseq==NULL) { GMessage("Warning: getUnspliced(NULL,.. ) called!\n");
        return NULL;
    }
    //restore normal coordinates:
    if (exons.Count()==0) return NULL;
    int fspan=end-start+1;
    const char* gsubseq=faseq->subseq(start, fspan);
    if (gsubseq==NULL) {
        GError("Error getting subseq for %s (%d..%d)!\n", gffID, start, end);
    }
    char* unspliced=NULL;

    int seqstart=exons.First()->start;
    int seqend=exons.Last()->end;

    int unsplicedlen = 0;
    if (seglst)
    	seglst->Clear(strand);
    unsplicedlen += seqend - seqstart + 1;

    GMALLOC(unspliced, unsplicedlen+1); //allocate more here
    //uint seqstart, seqend;
    int s = 0; //resulting nucleotide counter
    if (strand=='-')
    {
        if (seglst!=NULL)
            seglst->add(s+1,s+1+seqend-seqstart, seqstart, seqend);
        for (int i=seqend;i>=seqstart;i--)   {
            unspliced[s] = ntComplement(gsubseq[i-start]);
            s++;
        }//for each nt
    } // - strand
    else {
      // + strand
        if (seglst!=NULL)
            seglst->add(s+1,s+1+seqend-seqstart, seqstart, seqend);
        for (int i=seqstart;i<=seqend;i++) {
            unspliced[s]=gsubseq[i-start];
            s++;
        }//for each nt
    } // + strand
    //assert(s <= unsplicedlen);
    unspliced[s]=0;
    if (rlen!=NULL) *rlen=s;
    return unspliced;
}

 void GffObj::addPadding(int padLeft, int padRight) {
	 this->start-=padLeft;
	 this->end+=padRight;
	 if (exons.Count()>0) {
		 exons[0]->start-=padLeft;
		 exons.Last()->end+=padRight;
	 }
	 covlen+=padLeft+padRight;
 }

 void GffObj::removePadding(int padLeft, int padRight) {
	 this->start+=padLeft;
	 this->end-=padRight;
	 if (exons.Count()>0) {
		 exons[0]->start+=padLeft;
		 exons.Last()->end-=padRight;
	 }
	 covlen-=padLeft+padRight;
 }

char* GffObj::getSpliced(GFaSeqGet* faseq, bool CDSonly, int* rlen, uint* cds_start, uint* cds_end,
          GMapSegments* seglst, bool cds_open) {
	//cds_open only makes sense when CDSonly is true by overriding CDS 3'end such that the end of
	//the sequence beyond the 3' CDS end is also returned (the 3' UTR is appended to the CDS)
  if (CDSonly && CDstart==0) {
	  GMessage("Warning: getSpliced(CDSOnly) requested for transcript with no CDS (%s)!\n", gffID); //should never happen
	  return NULL;
  }
  if (faseq==NULL) {
	  GMessage("Warning: getSpliced() called with uninitialized GFaSeqGet object!\n"); //should never happen
      return NULL;
  }
  GList<GffExon>* xsegs=&exons;
  if (CDSonly && this->cdss!=NULL)
	  xsegs=this->cdss;
  if (xsegs->Count()==0) return NULL;
  int fspan=end-start+1;
  const char* gsubseq=faseq->subseq(start, fspan);
  if (gsubseq==NULL) {
        GError("Error getting subseq for %s (%d..%d)!\n", gffID, start, end);
  }
  if (fspan<(int)(end-start+1)) {
	  //special case: stop coordinate was extended past the gseq length, must adjust
     int endadj=end-start+1-fspan;
     uint prevend=end;
     end-=endadj;
     if (CDend>end) CDend=end;
     if (xsegs->Last()->end>end) {
         xsegs->Last()->end=end; //this could be trouble if exon start is also > end
         if (xsegs->Last()->start>xsegs->Last()->end) {
            GError("GffObj::getSpliced() error: improper genomic coordinate %d on %s for %s\n",
                  prevend,getGSeqName(), getID());
         }
         covlen-=endadj;
     }
  }
  char* spliced=NULL;
  GMALLOC(spliced, covlen+1); //IMPORTANT: covlen must be correct here!
  uint g_start=0, g_end=0;
  int cdsadj=0;
  if (CDphase=='1' || CDphase=='2') {
      cdsadj=CDphase-'0';
  }
  uint CDS_start=CDstart;
  uint CDS_stop=CDend;
  if (cdsadj>0) {
     if (strand=='-') CDS_stop-=cdsadj;
     else CDS_start+=cdsadj;
  }
  if (CDSonly) {
    g_start=CDS_start;
    g_end=CDS_stop;
    if (g_end-g_start<3)
    	 GMessage("Warning: CDS %d-%d too short for %s, check your data.\n",
    			 g_start, g_end, gffID);
  } else { //all exon content, not just CDS
    g_start=xsegs->First()->start;
    g_end=xsegs->Last()->end;
    cds_open=false; //override mistaken user request
  }
  if (seglst!=NULL) seglst->Clear(strand);
  int s=0; //resulting nucleotide counter
  if (strand=='-') {
    if (cds_open) {// appending 3'UTR
    	g_start=xsegs->First()->start;
    	//CDS_start=g_start;
    }
    for (int x=xsegs->Count()-1;x>=0;x--) {
       uint sgstart=xsegs->Get(x)->start;
       uint sgend=xsegs->Get(x)->end;
       if (g_end<sgstart || g_start>sgend) continue;
       if (g_start>=sgstart && g_start<=sgend)
          sgstart=g_start; //3' end within this segment
       if (g_end>=sgstart && g_end<=sgend)
          sgend=g_end; //5' end within this segment
       if (seglst!=NULL)
          seglst->add(s+1,s+1+sgend-sgstart,sgend,sgstart);
       for (uint i=sgend;i>=sgstart;i--) {
         spliced[s] = ntComplement(gsubseq[i-start]);
         s++;
       }//for each nt
       //--update local CDS start-end coordinates
       if (cds_start!=NULL && CDS_stop>=sgstart && CDS_stop<=sgend) {
         //CDS start in this segment
         *cds_start=s-(CDS_stop-sgstart);
       }
       if (cds_end!=NULL && CDS_start>=sgstart && CDS_start<=sgend) {
         //CDS stop in this segment
         *cds_end=s-(CDS_start-sgstart);
       }
      } //for each exon
    } // - strand
   else { // + strand
    if (cds_open) { // appending 3'UTR
      	g_end=xsegs->Last()->end;
      	//CDS_stop=g_end;
    }
    for (int x=0;x<xsegs->Count();x++) {
      uint sgstart=xsegs->Get(x)->start;
      uint sgend=xsegs->Get(x)->end;
      if (g_end<sgstart || g_start>sgend) continue;
      if (g_start>=sgstart && g_start<=sgend)
            sgstart=g_start; //seqstart within this segment
      if (g_end>=sgstart && g_end<=sgend)
            sgend=g_end; //seqend within this segment
      if (seglst!=NULL)
          seglst->add(s+1,s+1+sgend-sgstart, sgstart, sgend);
      for (uint i=sgstart;i<=sgend;i++) {
          spliced[s]=gsubseq[i-start];
          s++;
      }//for each nt
      //--update local CDS start-end coordinates
      if (cds_start!=NULL && CDS_start>=sgstart && CDS_start<=sgend) {
        //CDS start in this segment
        *cds_start=s-(sgend-CDS_start);
      }
      if (cds_end!=NULL && CDS_stop>=sgstart && CDS_stop<=sgend) {
        //CDS stop in this segment
        *cds_end=s-(sgend-CDS_stop);
      }
    } //for each exon
  } // + strand
  spliced[s]=0;
  if (rlen!=NULL) *rlen=s;
  return spliced;
}

void GffObj::printSummary(FILE* fout) {
 if (fout==NULL) fout=stdout;
 fprintf(fout, "%s\t%c\t%d\t%d\t", gffID,
          strand, start, end);
 gscore.print(fout);
 fprintf(fout, "\n");
}
//TODO we should also have an escapeChars function for some situations
//when we want to write a GFF3 strictly compliant to the dang specification
void GffObj::decodeHexChars(char* dbuf, const char* s, int maxlen) {
	int dlen=0;
	dbuf[0]=0;
	if (s==NULL) return;
	for (const char* p=s;(*p)!=0 && dlen<maxlen;++p) {
		if (p[0]=='%' && isxdigit(p[1]) && isxdigit(p[2])) {
			int a=*(++p);
			if (a>'Z') a^=0x20; //toupper()
			if (a>'9') a=10+(a-'A');
			      else a-='0';
			int b=*(++p);
			if (b>'Z') b^=0x20;
			if (b>'9') b=10+(b-'A');
			      else b-='0';
			char c=(char)((a<<4)+b);
			if (c=='%') {
				dbuf[dlen]='p';
				++dlen;
				dbuf[dlen]='r';
				++dlen;
				c='c';
			}
			else if (c==';') c='.';
			else if (c<='\t') c=' ';
			if (c>=' ') {
				dbuf[dlen]=c;
				++dlen;
				continue;
			}
		}
		dbuf[dlen]=*p;
		++dlen;
	}
	dbuf[dlen]=0;
}

void GffObj::printGTab(FILE* fout, char** extraAttrs) {
	fprintf(fout, "%s\t%c\t%d\t%d\t%s\t", this->getGSeqName(), this->strand,
			this->start, this->end, this->getID());
	if (exons.Count()) printExonList(fout);
	else fprintf(fout, ".");
	if (extraAttrs!=NULL) {
		//print a list of "attr=value;" pairs here as the last column
		//for requested attributes
		bool t1=true;
		for (int i=0;extraAttrs[i]!=NULL;++i) {
			const char* v=this->getAttr(extraAttrs[i]);
			if (v==NULL) continue;
			if (t1) { fprintf(fout, "\t"); t1=false; }
			fprintf(fout, "%s=%s;", extraAttrs[i], v);
		}
	}
	fprintf(fout,"\n");
}

void GffObj::printGxfExon(FILE* fout, const char* tlabel, const char* gseqname, bool iscds,
                             GffExon* exon, bool gff3, bool cvtChars,
							 char* dbuf, int dbuf_len) {
  //strcpy(dbuf,".");
  //if (exon->score>0) sprintf(dbuf,"%.2f", exon->score);
  exon->score.sprint(dbuf);
  if (exon->phase==0 || !iscds) exon->phase='.';
  const char* ftype=iscds ? "CDS" : getSubfName();
  const char* attrname=NULL;
  const char* attrval=NULL;
  if (gff3) {
    fprintf(fout,
      "%s\t%s\t%s\t%d\t%d\t%s\t%c\t%c\tParent=%s",
      gseqname, tlabel, ftype, exon->start, exon->end, dbuf, strand,
	  exon->phase, gffID);
    if (exon->attrs!=NULL) {
      for (int i=0;i<exon->attrs->Count();i++) {
        if (exon->attrs->Get(i)->cds!=iscds) continue;
        attrname=names->attrs.getName(exon->attrs->Get(i)->attr_id);
        if (cvtChars) {
          decodeHexChars(dbuf, exon->attrs->Get(i)->attr_val, dbuf_len-1);
          fprintf(fout,";%s=%s", attrname, dbuf);
        } else {
          fprintf(fout,";%s=%s", attrname, exon->attrs->Get(i)->attr_val);
        }
      }
    }
    fprintf(fout, "\n");
    } //GFF3
  else {//GTF
    fprintf(fout, "%s\t%s\t%s\t%d\t%d\t%s\t%c\t%c\ttranscript_id \"%s\";",
           gseqname, tlabel, ftype, exon->start, exon->end, dbuf, strand, exon->phase, gffID);
    if (geneID)
      fprintf(fout," gene_id \"%s\";",geneID);
    if (gene_name!=NULL) {
       fprintf(fout," gene_name \"%s\";",gene_name);
    }
    if (exon->attrs!=NULL) {
       bool trId=false;
       bool gId=false;
       for (int i=0;i<exon->attrs->Count();i++) {
            if (exon->attrs->Get(i)->attr_val==NULL) continue;
            if (exon->attrs->Get(i)->cds!=iscds) continue;
            attrname=names->attrs.getName(exon->attrs->Get(i)->attr_id);
            if (strcmp(attrname, "transcriptID")==0) {
            	if (trId) continue;
            	trId=true;
            }
            if (strcmp(attrname, "transcript_id")==0 && !trId) {
            	attrname="transcriptID";
            	trId=true;
            }
            if (strcmp(attrname, "geneID")==0) {
            	if (gId) continue;
            	gId=true;
            }
            if (strcmp(attrname, "gene_id")==0 && !gId) {
            	attrname="geneID";
            	gId=true;
            }
            if (Gstricmp(attrname, "gene_name")==0 && gene_name!=NULL) {
            	continue;
            }
            fprintf(fout, " %s ",attrname);
            if (cvtChars) {
              decodeHexChars(dbuf, exon->attrs->Get(i)->attr_val, dbuf_len-1);
              attrval=dbuf;
            } else {
              attrval=exon->attrs->Get(i)->attr_val;
            }

            if (attrval[0]=='"') fprintf(fout, "%s;",attrval);
                           else fprintf(fout, "\"%s\";",attrval);
        }
    }
    //for GTF, also append the GffObj attributes to each exon line
    // - do not do this when the transcript line is also printed!
    /*
    if (attrs!=NULL) {
         for (int i=0;i<attrs->Count();i++) {
            if (attrs->Get(i)->attr_val==NULL) continue;
            attrname=names->attrs.getName(attrs->Get(i)->attr_id);
            fprintf(fout, " %s ",attrname);
            if (cvtChars) {
              decodeHexChars(dbuf, attrs->Get(i)->attr_val, dbuf_len-1);
              attrval=dbuf;
            } else {
              attrval=attrs->Get(i)->attr_val;
            }
            if (attrval[0]=='"') fprintf(fout, "%s;",attrval);
                           else fprintf(fout, "\"%s\";",attrval);
         }
    }
    */
    fprintf(fout, "\n");
 }//GTF
}

void GffObj::printGxf(FILE* fout, GffPrintMode gffp,
                   const char* tlabel, const char* gfparent, bool cvtChars) {
 const int DBUF_LEN=1024; //there should not be attribute values longer than 1K!
 char dbuf[DBUF_LEN];
 if (tlabel==NULL) {
    tlabel=track_id>=0 ? names->tracks.Get(track_id)->name :
         (char*)"gffobj" ;
    }
 if (gffp==pgffBED) {
	 printBED(fout, cvtChars, dbuf, DBUF_LEN);
	 return;
 }
 const char* gseqname=names->gseqs.Get(gseq_id)->name;
 bool gff3 = (gffp>=pgffAny && gffp<=pgffTLF);
 bool showCDS = (gffp==pgtfAny || gffp==pgtfCDS || gffp==pgffCDS || gffp==pgffAny || gffp==pgffBoth);
 bool showExon = (gffp<=pgtfExon || gffp==pgffAny || gffp==pgffExon || gffp==pgffBoth);
 //if (gscore>0.0) sprintf(dbuf,"%.2f", gscore);
 //       else strcpy(dbuf,".");
 gscore.sprint(dbuf);
 if (gffp<=pgtfCDS && gffp>=pgtfAny) { //GTF output
	   fprintf(fout,
	     "%s\t%s\ttranscript\t%d\t%d\t%s\t%c\t.\ttranscript_id \"%s\"",
	     gseqname, tlabel, start, end, dbuf, strand, gffID);
	   char* gid=NULL;
	   if (geneID!=NULL) {
	      gid=geneID;
	   }
	   else {
		   gid=getAttr("gene_id");
		   if (gid==NULL)
			   gid=gffID; //last resort, write gid the same with gffID
	   }
	   if (gid!=NULL) fprintf(fout, "; gene_id \"%s\"",gid);
	   if (gene_name!=NULL && getAttr("gene_name")==NULL && getAttr("GENE_NAME")==NULL)
	      fprintf(fout, "; gene_name \"%s\"",gene_name);
	   if (attrs!=NULL) {
		    bool trId=false;
		    //bool gId=false;
		    for (int i=0;i<attrs->Count();i++) {
		      const char* attrname=names->attrs.getName(attrs->Get(i)->attr_id);
		      const char* attrval=attrs->Get(i)->attr_val;
		      if (attrval==NULL || attrval[0]=='\0') continue;
		      if (strcmp(attrname, "transcriptID")==0) {
	            	if (trId) continue;
	            	trId=true;
		      }
		      if (strcmp(attrname, "transcript_id")==0 && !trId) {
	            	attrname="transcriptID";
	            	trId=true;
		      }
		      if (Gstrcmp(attrname, "geneID")==0 && gid!=NULL &&
	            		strcmp(attrval, gid)==0) continue;
		      if (strcmp(attrname, "gene_id")==0) continue;
		      if (cvtChars) {
		    	  decodeHexChars(dbuf, attrval, DBUF_LEN-1);
		    	  fprintf(fout,"; %s \"%s\"", attrname, dbuf);
		      }
		      else
		    	 fprintf(fout,"; %s \"%s\"", attrname, attrs->Get(i)->attr_val);
		    }
	   }
	   fprintf(fout,";\n");
 }
 else if (gff3) {
   //print GFF3 transcript line:
   uint pstart, pend;
   if (gffp==pgffCDS) {
      pstart=CDstart;
      pend=CDend;
   }
   else { pstart=start;pend=end; }
   //const char* ftype=isTranscript() ? "mRNA" : getFeatureName();
   const char* ftype=getFeatureName();
   fprintf(fout,
     "%s\t%s\t%s\t%d\t%d\t%s\t%c\t.\tID=%s",
     gseqname, tlabel, ftype, pstart, pend, dbuf, strand, gffID);
   bool parentPrint=false;
   if (gfparent!=NULL && gffp!=pgffTLF) {
      //parent override - also prevents printing gene_name and gene_id
      fprintf(fout, ";Parent=%s",gfparent);
      parentPrint=true;
   }
   else if (parent!=NULL && !parent->isDiscarded() && gffp!=pgffTLF) {
           fprintf(fout, ";Parent=%s",parent->getID());
           if (parent->isGene()) parentPrint=true;
   }
   if (gffp==pgffTLF) {
	   fprintf(fout, ";exonCount=%d",exons.Count());
	   if (exons.Count()>0)
		   fprintf(fout, ";exons=%d-%d", exons[0]->start, exons[0]->end);
	   for (int i=1;i<exons.Count();++i) {
		   fprintf(fout, ",%d-%d",exons[i]->start, exons[i]->end);
	   }
   }
   if (CDstart>0 && (gffp==pgffTLF || !showCDS)) {
	   if (cdss==NULL) fprintf(fout,";CDS=%d:%d",CDstart,CDend);
	   else {
		   fprintf(fout, ";CDS=");
		   for (int i=0;i<cdss->Count();++i) {
			   if (i>0) fprintf(fout, ",");
			   fprintf(fout, "%d-%d", (*cdss)[i]->start, (*cdss)[i]->end);
		   }
	   }
   }
   if (CDphase>0 && (gffp==pgffTLF || !showCDS)) fprintf(fout,";CDSphase=%c", CDphase);
   char* g_id=NULL;
   if (geneID!=NULL && !parentPrint && getAttr("geneID")==NULL &&
		   ((g_id=getAttr("gene_id"))==NULL || strcmp(g_id, geneID)!=0))
      fprintf(fout, ";geneID=%s",geneID);
   if (gene_name!=NULL && !parentPrint && getAttr("gene_name")==NULL && getAttr("GENE_NAME")==NULL)
      fprintf(fout, ";gene_name=%s",gene_name);
   if (attrs!=NULL) {
	    for (int i=0;i<attrs->Count();i++) {
	      const char* attrname=names->attrs.getName(attrs->Get(i)->attr_id);
	      const char* attrval=attrs->Get(i)->attr_val;
	      if (attrval==NULL || attrval[0]=='\0') continue;
	    	  //fprintf(fout,";%s",attrname);
	      if (cvtChars) {
	    	  decodeHexChars(dbuf, attrval, DBUF_LEN-1);
	    	  fprintf(fout,";%s=%s", attrname, dbuf);
	      }
	      else
	    	 fprintf(fout,";%s=%s", attrname, attrs->Get(i)->attr_val);
	    }
   }
   fprintf(fout,"\n");
 }// gff3 transcript line
 if (gffp==pgffTLF) return;
 bool is_cds_only = (gffp==pgffBoth) ? false : isCDSOnly();
 if (showExon) {
    //print exons
    for (int i=0;i<exons.Count();i++) {
      printGxfExon(fout, tlabel, gseqname, is_cds_only, exons[i], gff3, cvtChars, dbuf, DBUF_LEN);
    }
 }//printing exons
 if (showCDS && !is_cds_only && CDstart>0) {
	GVec<GffExon> cds;
	getCDSegs(cds); //also uses/prepares the CDS phase for each CDS segment
	for (int i=0;i<cds.Count();i++) {
		printGxfExon(fout, tlabel, gseqname, true, &(cds[i]), gff3, cvtChars, dbuf, DBUF_LEN);
	}
  } //showCDS
}

void GffObj::updateCDSPhase(GList<GffExon>& segs) {
  int cdsacc=0;
  if (CDphase=='1' || CDphase=='2') {
      cdsacc+= 3-(CDphase-'0');
  }
  else CDphase='0';
  if (strand=='-') { //reverse strand
     for (int i=segs.Count()-1;i>=0;i--) {
         segs[i]->phase='0'+ (3-cdsacc%3)%3;
         cdsacc+=segs[i]->end-segs[i]->start+1;
     }
  }
    else { //forward strand
     for (int i=0;i<segs.Count();i++) {
         segs[i]->phase='0'+ (3-cdsacc%3)%3;
         cdsacc+=segs[i]->end-segs[i]->start+1;
     }
  }
}

void GffObj::getCDSegs(GVec<GffExon>& cds) {
  //like updateCDSPhase() above, also updates phase for each segment
  GffExon cdseg(true);
  cds.Clear();
  if (cdss!=NULL) {
	//copy directly from cdss list
	for (int i=0;i<cdss->Count();i++) {
		cdseg=(*cdss->Get(i));
		cdseg.sharedAttrs=true;
		cds.Add(cdseg);
	}
    return;
  }
  int cdsacc=0;
  if (CDphase=='1' || CDphase=='2') {
      cdsacc+= 3-(CDphase-'0');
  }
  if (strand=='-') {
     for (int x=exons.Count()-1;x>=0;x--) {
        uint sgstart=exons[x]->start;
        uint sgend=exons[x]->end;
        if (CDend<sgstart || CDstart>sgend) continue;
        if (CDstart>=sgstart && CDstart<=sgend)
              sgstart=CDstart; //cdstart within this segment
        if (CDend>=sgstart && CDend<=sgend)
              sgend=CDend; //cdend within this segment
        cdseg.start=sgstart;
        cdseg.end=sgend;
        //cdseg.phase='0'+(cdsacc>0 ? (3-cdsacc%3)%3 : 0);
        cdseg.phase='0'+ (3-cdsacc%3)%3;
        cdsacc+=sgend-sgstart+1;
        cdseg.attrs=exons[x]->attrs;
        cdseg.sharedAttrs=true;
        cds.Add(cdseg);
       } //for each exon
     cds.Reverse();
     } // - strand
    else { // + strand
     for (int x=0;x<exons.Count();x++) {
       uint sgstart=exons[x]->start;
       uint sgend=exons[x]->end;
       if (CDend<sgstart || CDstart>sgend) continue;
       if (CDstart>=sgstart && CDstart<=sgend)
             sgstart=CDstart; //seqstart within this segment
       if (CDend>=sgstart && CDend<=sgend)
             sgend=CDend; //seqend within this segment
       cdseg.start=sgstart;
       cdseg.end=sgend;
       //cdseg.phase='0'+(cdsacc>0 ? (3-cdsacc%3)%3 : 0);
       cdseg.phase='0' + (3-cdsacc%3)%3 ;
       cdsacc+=sgend-sgstart+1;
       cdseg.attrs=exons[x]->attrs;
       cdseg.sharedAttrs=true;
       cds.Add(cdseg);
       } //for each exon
   } // + strand

}

//-- transcript match/overlap classification functions


char transcriptMatch(GffObj& a, GffObj& b, int& ovlen) {
	//return '=' if exact exon match, '~' if intron-chain match (or 80% overlap for single-exon)
	// or 0 otherwise
	int imax=a.exons.Count()-1;
	int jmax=b.exons.Count()-1;
	ovlen=0;
	if (imax!=jmax) return false; //different number of exons, cannot match
	if (imax==0) //single-exon mRNAs
	    return (singleExonTMatch(a,b,ovlen));
	if ( a.exons[imax]->start<b.exons[0]->end ||
		b.exons[jmax]->start<a.exons[0]->end )
		return 0; //intron chains do not overlap at all
	//check intron overlaps
	ovlen=a.exons[0]->end-(GMAX(a.start,b.start))+1;
	ovlen+=(GMIN(a.end,b.end))-a.exons.Last()->start;
	for (int i=1;i<=imax;i++) {
		if (i<imax) ovlen+=a.exons[i]->len();
		if ((a.exons[i-1]->end!=b.exons[i-1]->end) ||
			(a.exons[i]->start!=b.exons[i]->start)) {
			return 0; //intron mismatch
		}
	}
	//--- full intron chain match:
	if (a.exons[0]->start==b.exons[0]->start &&
		a.exons.Last()->end==b.exons.Last()->end)
		   return '=';
	return '~';
}


char singleExonTMatch(GffObj& m, GffObj& r, int& ovlen) {
 //return '=' if exact match, '~' if the overlap is >=80% of the longer sequence length
 // return 0 if there is no overlap
 GSeg mseg(m.start, m.end);
 ovlen=mseg.overlapLen(r.start,r.end);
 if (ovlen<=0) return 0;
 // fuzzy matching for single-exon transcripts:
 // matching = overlap is at least 80% of the length of the longer transcript
 // *OR* in case of reverse containment (reference contained in m)
 //   it's also considered "matching" if the overlap is at least 80% of
 //   the reference len AND at least 70% of the query len
 if (m.start==r.start && m.end==r.end) return '=';
 if (m.covlen>r.covlen) {
   if ( (ovlen >= m.covlen*0.8) ||
		   (ovlen >= r.covlen*0.8 && ovlen >= m.covlen* 0.7 ) )
		   //allow also some fuzzy reverse containment
           return '~';
 } else {
   if (ovlen >= r.covlen*0.8) return '~';
 }
 return 0;
}

//formerly in gffcompare
char getOvlCode(GffObj& m, GffObj& r, int& ovlen, bool strictMatch) {
	ovlen=0; //total actual exonic overlap
	if (!m.overlap(r.start,r.end)) return 0;
	int jmax=r.exons.Count()-1;
	//int iovlen=0; //total m.exons overlap with ref introns
	char rcode=0;
	if (m.exons.Count()==1) { //single-exon transfrag
		GSeg mseg(m.start, m.end);
		if (jmax==0) { //also single-exon ref
			//ovlen=mseg.overlapLen(r.start,r.end);
			char eqcode=0;
			if ((eqcode=singleExonTMatch(m, r, ovlen))>0) {
				if (strictMatch) return eqcode;
				            else return '=';
			}
			if (m.covlen<r.covlen)
			   { if (ovlen >= m.covlen*0.8) return 'c'; } // fuzzy containment
			else
				if (ovlen >= r.covlen*0.8 ) return 'k';   // fuzzy reverse containment
			return 'o'; //just plain overlapping
		}
		//-- single-exon qry overlaping multi-exon ref
		//check full pre-mRNA case (all introns retained): code 'm'
		if (m.start<=r.exons[0]->end && m.end>=r.exons[jmax]->start)
			return 'm';

		for (int j=0;j<=jmax;j++) {
			//check if it's ~contained by an exon
			int exovlen=mseg.overlapLen(r.exons[j]);
			if (exovlen>0) {
				ovlen+=exovlen;
				if (m.start>r.exons[j]->start-4 && m.end<r.exons[j]->end+4) {
					return 'c'; //close enough to be considered contained in this exon
				}
			}
			if (j==jmax) break; //last exon here, no intron to check
			//check if it fully covers an intron (retained intron)
			if (m.start<r.exons[j]->end && m.end>r.exons[j+1]->start)
				return 'n';
			//check if it's fully contained by an intron
			if (m.end<r.exons[j+1]->start && m.start>r.exons[j]->end)
				return 'i';
			// check if it's a potential pre-mRNA transcript
			// (if overlaps this intron at least 10 bases)
			uint introvl=mseg.overlapLen(r.exons[j]->end+1, r.exons[j+1]->start-1);
			//iovlen+=introvl;
			if (introvl>=10 && mseg.len()>introvl+10) { rcode='e'; }
		} //for each ref exon
		if (rcode>0) return rcode;
		return 'o'; //plain overlap, uncategorized
	} //single-exon transfrag
	//-- multi-exon transfrag --
	int imax=m.exons.Count()-1;// imax>0 here
	if (jmax==0) { //single-exon reference overlap
		//any exon overlap?
		GSeg rseg(r.start, r.end);
		for (int i=0;i<=imax;i++) {
			//check if it's ~contained by an exon
			int exovlen=rseg.overlapLen(m.exons[i]);
			if (exovlen>0) {
				ovlen+=exovlen;
				if (r.start>m.exons[i]->start-4 && r.end<m.exons[i]->end+4) {
					return 'k'; //reference contained in this assembled exon
				}
			}
			if (i==imax) break;
			if (r.end<m.exons[i+1]->start && r.start>m.exons[i]->end)
				return 'y'; //ref contained in this transfrag intron
		}
		return 'o';
	}
	// * check if transfrag contained by a ref intron
	for (int j=0;j<jmax;j++) {
		if (m.end<r.exons[j+1]->start && m.start>r.exons[j]->end)
			return 'i';
	}
	if (m.exons[imax]->start<r.exons[0]->end) {
		//qry intron chain ends before ref intron chain starts
		//check if last qry exon plugs the 1st ref intron
		if (m.exons[imax]->start<=r.exons[0]->end &&
			m.exons[imax]->end>=r.exons[1]->start) return 'n';
		return 'o'; //only terminal exons overlap
	}
	else if (r.exons[jmax]->start<m.exons[0]->end) {
		//qry intron chain starts after ref intron chain ends
		//check if first qry exon plugs the last ref intron
		if (m.exons[0]->start<=r.exons[jmax-1]->end &&
			m.exons[0]->end>=r.exons[jmax]->start) return 'n';
		return 'o'; //only terminal exons overlap
	}
	//check intron chain overlap (match, containment, intron retention etc.)
	int i=1; //index of exon to the right of current qry intron
	int j=1; //index of exon to the right of current ref intron
	bool intron_conflict=false; //overlapping introns have at least a mismatching splice site
	//from here on we check all qry introns against ref introns
	bool junct_match=false; //true if at least a junction match is found
	bool ichain_match=false; //if there is intron (sub-)chain match, to be updated by any mismatch
	bool intron_ovl=false; //if any intron overlap is found
	bool intron_retention=false; //if any ref intron is covered by a qry exon
	//intron chain (partial) match exon-index boundaries:
	int imfirst=0; //index of exon after first intron match in query (valid>0)
	int jmfirst=0; //index of exon after first intron match in reference (valid>0)
	int imlast=0;  //index of exon after last intron match in query
	int jmlast=0;  //index of  exon after last intron match in reference
	//--keep track of the last overlapping introns in both qry and ref:
	//int q_last_iovl=0;
	//int r_last_iovl=0;

	//check for intron matches
	while (i<=imax && j<=jmax) {
		uint mstart=m.exons[i-1]->end; //qry intron start-end
		uint mend=m.exons[i]->start;
		uint rstart=r.exons[j-1]->end; //ref intron start-end
		uint rend=r.exons[j]->start;
		if (rend<mstart) { //qry intron starts after ref intron ends
			if (!intron_conflict && r.exons[j]->overlap(mstart+1, mend-1))
				intron_conflict=true; //next ref exon overlaps this qry intron
			if (!intron_retention && rstart>=m.exons[i-1]->start && rend<=m.exons[i-1]->end)
				intron_retention=true; //this ref intron is covered by previous qry exons[i-1]
			if (intron_ovl) ichain_match=false;
			j++;
			continue;
		} //no intron overlap, skipping ref intron
		if (rstart>mend) { //qry intron ends before ref intron starts
			//if qry intron overlaps the exon on the left, we have an intron conflict
			if (!intron_conflict && r.exons[j-1]->overlap(mstart+1, mend-1))
				intron_conflict=true;
			if (!intron_retention && rstart>=m.exons[i]->start && rend<=m.exons[i]->end)
				intron_retention=true;
			if (intron_ovl) ichain_match=false;
			i++;
			continue;
		} //no intron overlap, skipping qry intron
		intron_ovl=true;
		//q_last_iovl=i; //keep track of the last overlapping introns in both qry and ref
		//r_last_iovl=j;
		//overlapping introns, test junction matching
		bool smatch=(mstart==rstart);
		bool ematch=(mend==rend);
		if (smatch || ematch) junct_match=true;
		if (smatch && ematch) {
			//perfect match for this intron
			if (jmfirst==0) {
				ichain_match=true;
				jmfirst=j;
				imfirst=i;
			}
			if (ichain_match) {
  		       imlast=i;
			   jmlast=j;
			}
			i++; j++;
			continue;
		}
		//intron overlapping but not fully matching
		intron_conflict=true;
		ichain_match=false;
		if (mend>rend) j++; else i++;
	} //while checking intron overlaps
	/*** additional checking needed for intron retention when there is no ichain_match or overlap ?
    if (!intron_retention && r_last_iovl<jmax) {
	   //-- check the remaining ref introns not checked yet for retention
       int i=q_last_iovl;
       for (int j=r_last_iovl+1;j<=jmax && i<=imax;++j) {
   		uint rstart=r.exons[j-1]->end; //ref intron start-end
   		uint rend=r.exons[j]->start;
   		if (rend<m.exons[i]->start) {
   			i++;
   			continue;
   		}
   		if (rstart>m.exons[i]->end)
   			continue;
   		//overlap between ref intron and m.exons[i]
   		if (rstart>=m.exons[i]->start && rend<=m.exons[i]->end) {
   			intron_retention=true;
   			break;
   		}
       }
    }
    ***/
	// --- when qry intron chain is contained within ref intron chain
	//     qry terminal exons may poke (overhang) into ref's other introns
	int l_iovh=0;   // overhang of q left boundary beyond the end of ref intron on the left
	int r_iovh=0;   // same type of overhang through the ref intron on the right
	int qry_intron_poking=0;
	// --- when ref intron chain is contained within qry intron chain,
	//     terminal exons of ref may poke (overhang) into qry other introns
	int l_jovh=0;   // overhang of q left boundary beyond the end of ref intron to the left
	int r_jovh=0;   // same type of overhang through the ref intron on the right
	int ref_intron_poking=0;
	if (ichain_match) { //intron (sub-)chain compatible so far (but there could still be conflicts)
		if (imfirst==1 && imlast==imax) { // qry full intron chain match
			if (jmfirst==1 && jmlast==jmax) {//identical intron chains
				if (strictMatch) return (r.exons[0]->start==m.exons[0]->start &&
						              r.exons.Last()->end && m.exons.Last()->end) ? '=' : '~';
				else return '=';
			}
			// -- a partial intron chain match
			if (jmfirst>1) {
				//find if m.start falls within any ref intron before jmfirst
				for (int j=jmfirst-1;j>0;--j)
					if (m.start<r.exons[j]->start) {
						if (m.start>r.exons[j-1]->end) { //m.start within this ref intron
							l_iovh = r.exons[j]->start - m.start;
							break;
						}
						else { intron_retention=true; ichain_match=false; }
					}
			}
			if (jmlast<jmax) {
				for (int j=jmlast;j<jmax;++j)
					if (m.end > r.exons[j]->end) {
						if (m.end < r.exons[j+1]->start) { //m.end within this ref intron
							r_iovh = m.end - r.exons[j]->end;
						    break;
						}
						else { intron_retention=true; ichain_match=false; }
					}
			}
			if (ichain_match && l_iovh<4 && r_iovh<4) return 'c';
			qry_intron_poking=GMAX(l_iovh, r_iovh);
		} else if ((jmfirst==1 && jmlast==jmax)) {//ref intron chain match
			//check if the reference j-chain is contained in qry i-chain
			//check for ref ends poking into qry introns
			if (imfirst>1)  {
				for (int i=imfirst-1;i>0;--i)
					if (m.exons[i]->start>r.start) {
						if (r.start>m.exons[i-1]->end) {
							l_jovh = m.exons[i]->start - r.start;
							break;
						}
						else { ichain_match = false; }
					}
			}
			if (imlast<imax) {
				for (int i=imlast;i<imax;++i)
					if (r.end > m.exons[i]->end) {
						if (r.end < m.exons[i+1]->start)
							 { r_jovh = r.end - m.exons[i]->end; break; }
						else { ichain_match = false; }
					}
			}
			if (ichain_match && l_jovh<4 && r_jovh<4) return 'k'; //reverse containment
			ref_intron_poking=GMAX(l_jovh, r_jovh);
		}
	}
	//'=', 'c' and 'k' were checked and assigned, check for 'm' and 'n' before falling back to 'j'
	if (intron_retention) {
			//ref is boundary contained with qry intron chain ? that's not required for 'm'
		    //GMessage("r_jovh=%d, r_iovh=%d, l_jovh=%d, l_iovh=%d\n", r_jovh, r_iovh, l_jovh, l_iovh);
		    //GMessage("m.start=%d, r.exons[0]->end=%d, m.end=%d, r.exons[jmax]->start=%d\n",
		    //		m.start, r.exons[0]->end, m.end, r.exons[jmax]->start);
		    //if (ref_intron_poking>0 && )
		//we just need to have no intron poking going on
		if (!intron_conflict && ref_intron_poking<4 && qry_intron_poking<4) return 'm';
		else return 'n';
	}
	if (junct_match) return 'j';
	//we could have 'o' or 'y' here
	//any real exon overlaps?
	ovlen=m.exonOverlapLen(r);
	if (ovlen>4) return 'o';
	return 'y'; //all reference exons are within transfrag introns!
}
