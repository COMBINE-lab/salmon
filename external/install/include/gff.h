#ifndef GFF_H
#define GFF_H

//#define CUFFLINKS 1

#include "GBase.h"
#include "gdna.h"
#include "codons.h"
#include "GFaSeqGet.h"
#include "GList.hh"
#include "GHash.hh"

#ifdef CUFFLINKS
#include <boost/crc.hpp>  // for boost::crc_32_type
#endif

//reserved Gffnames::feats entries -- basic feature types
extern int gff_fid_mRNA; // "mRNA" feature name
extern int gff_fid_transcript; // *RNA, *transcript feature name
extern int gff_fid_exon;

extern const uint GFF_MAX_LOCUS;
extern const uint GFF_MAX_EXON;
extern const uint GFF_MAX_INTRON;
extern const int CLASSCODE_OVL_RANK;
//extern const uint gfo_flag_LEVEL_MSK; //hierarchical level: 0 = no parent
//extern const byte gfo_flagShift_LEVEL;

//extern bool gff_show_warnings;

#define GFF_LINELEN 4096
#define ERR_NULL_GFNAMES "Error: GffObj::%s requires a non-null GffNames* names!\n"


enum GffExonType {
  exgffIntron=-1, // useless "intron" feature
  exgffNone=0,  //not recognizable or unitialized exonic segment
  exgffStartCodon, //from "start_codon" feature (within CDS)
  exgffStopCodon, //from "stop_codon" feature (may be outside CDS, but should)
  exgffCDS,  //from "CDS" feature
  exgffUTR,  //from "UTR" feature
  exgffCDSUTR, //from a merge of UTR and CDS feature
  exgffExon, //from "exon" feature
};

extern const char* exonTypes[];

const char* strExonType(char xtype);
class GfList;

typedef void GFFCommentParser(const char* cmline, GfList* gflst); //comment parser callback
//Useful for parsing/maintaining ref seq info from comment lines like this:
//##sequence-region chr1 1 24895642

class GffReader;
class GffObj;

//---transcript overlapping - utility functions:
int classcode_rank(char c); //returns priority value for class codes

char getOvlCode(GffObj& m, GffObj& r, int& ovlen, bool strictMatch=false); //returns: class code

char transcriptMatch(GffObj& a, GffObj& b, int& ovlen); //generic transcript match test
// -- return '=', '~'  or 0
char singleExonTMatch(GffObj& m, GffObj& r, int& ovlen); //single-exon transcript match test

//---
// -- tracking exon/CDS segments from local mRNA to genome coordinates
class GMapSeg:public GSeg {
public:
	uint gstart; //genome start location
	uint gend;   //genome end location
	//gend<gstart when mapped on reverse strand !
	GMapSeg(uint s=0, uint e=0, uint gs=0, uint ge=0):GSeg(s,e),
			 gstart(gs), gend(ge) {};
    int g_within(uint gc) {
    	//return 0 if not within gstart-gend intervals
    	//or offset from gstart otherwise (always positive)
    	if (gstart>gend) { //reverse strand mapping
    		if  (gc<gend || gc>gstart) return 0;
			return (gstart-gc);
    	}
    	else {
    		if (gc<gstart || gc>gend) return 0;
    		return (gc-gstart);
    	}
    }
};

struct GffScore {
	float score;
	int8_t precision;
	GffScore(float sc=0, int8_t prec=-1):score(sc),precision(prec) { }
	void print(FILE* outf) {
		if (precision<0) fprintf(outf, ".");
		else fprintf(outf, "%.*f", precision, score);
	}
	void sprint(char* outs) {
		if (precision<0) sprintf(outs, ".");
		else sprintf(outs, "%.*f", precision, score);
	}
	bool operator<(GffScore& v) {
		return this->score<v.score;
	}
	bool operator<=(GffScore& v) {
		return this->score<=v.score;
	}
	bool operator>(GffScore& v) {
		return this->score>v.score;
	}
	bool operator>=(GffScore& v) {
		return this->score>=v.score;
	}
	bool operator==(GffScore& v) {
		return this->score==v.score;
	}
};

extern const GffScore GFFSCORE_NONE;

class GMapSegments:public GVec<GMapSeg> {
  public:
	int dir; //-1 or +1 (reverse/forward for genome coordinates)
	GSeg lreg; // always 1,max local coord
	GSeg greg; // genomic min,max coords
	GMapSegments(char strand='+'):lreg(0,0),greg(0,0) {
		dir=(strand=='-') ? -1 : 1;
	}
	void Clear(char strand='+') {
		lreg.start=0;lreg.end=0;
		greg.start=0;greg.end=0;
		dir = (strand=='-') ? -1 : 1;;
		GVec<GMapSeg>::Clear();
	}
    int add(uint s, uint e, uint gs, uint ge) {
    	if (dir<0) {
    		if (gs<ge) {
    			Gswap(gs, ge);
    		}
    		if (gs>greg.end) greg.end=gs;
    		if (ge<greg.start || greg.start==0) greg.start=ge;
    	} else {
    		if (ge>greg.end) greg.end=ge;
    		if (gs<greg.start || greg.start==0) greg.start=gs;
    	}
    	GMapSeg gm(s, e, gs, ge);
		if (gm.end>lreg.end) lreg.end=gm.end;
		if (gm.start<lreg.start || lreg.start==0) lreg.start=gm.start;
    	return GVec<GMapSeg>::Add(gm);
    }
    uint gmap(uint lc) { //takes a local coordinate and returns its mapping to genomic coordinates
    	//returns 0 if mapping cannot be performed!
    	if (lc==0 || fCount==0 || lc<lreg.start || lc>lreg.end) return 0;
    	//find local segment containing this coord
    	int i=0;
    	while (i<fCount) {
    		if (lc>=fArray[i].start && lc<=fArray[i].end)
    			return (fArray[i].gstart+dir*(lc-fArray[i].start));
    		++i;
        }
    	return 0;
    }
    uint lmap(uint gc) { //takes a genome coordinate and returns its mapping to local coordinates
    	if (gc==0 || fCount==0 || gc<greg.start || gc>greg.end) return 0;
    	//find genomic segment containing this coord
    	int i=0;
    	while (i<fCount) {
    		int ofs=fArray[i].g_within(gc);
    		if (ofs!=0)
    			return (fArray[i].start+ofs);
    		++i;
        }
    	return 0;
    }
};

//reading a whole transcript from a BED-12 line
class BEDLine {
 public:
    bool skip;
    char* dupline; //duplicate of original line
    char* line; //this will have tabs replaced by \0
    int llen;
    char* gseqname;
    uint fstart;
    uint fend;
    char strand;
    char* ID; //transcript ID from BED-12 (4th column)
    char* info; //13th column - these could be GFF3 attributes
    uint cds_start;
    uint cds_end;
    char cds_phase;
    GVec<GSeg> exons;
    BEDLine(GffReader* r=NULL, const char* l=NULL);
    ~BEDLine() {
    	GFREE(dupline);
    	GFREE(line);
    }
};

class GffLine {
 protected:
    char* _parents; //stores a copy of the Parent attribute value,
       //with commas replaced by \0
    int _parents_len;
    bool parseSegmentList(GVec<GSeg>& segs, char* str);
 public:
    char* dupline; //duplicate of original line
    char* line; //this will have tabs replaced by \0
    int llen;
    char* gseqname;
    char* track;
    char* ftype; //feature name: mRNA/gene/exon/CDS
    int ftype_id;
    char* info; //the last, attributes' field, unparsed
    uint fstart;
    uint fend;
    /*
    uint qstart; //overlap coords on query, if available
    uint qend;
    uint qlen; //query len, if given
    */
    float score;
    int8_t score_decimals;
    char strand;
    union {
    	unsigned int flags;
    	struct {
    		bool is_exonlike:2; //CDS,codon, UTR, exon
    	};
    	struct {
    	    bool is_cds:1; //"cds" or "start/stop_codon" features
    	    bool is_exon:1; //"exon" and "utr" features
    	    bool is_transcript:1; //if current feature is *RNA or *transcript
    	    bool is_gene:1; //current feature is *gene
    	    //bool is_gff3:1; //line appears to be in GFF3 format (0=GTF)
    	    bool is_gtf_transcript:1; //GTF transcript line with Parents parsed from gene_id
    	    bool skipLine:1;
    	    bool gffWarnings:1;
    	    bool is_gene_segment:1; //for NCBI's D/J/V/C_gene_segment
    	};
    };
    int8_t exontype; // gffExonType
    char phase;  // '.' , '0', '1' or '2', can be also given as CDSphase attribute in TLF
    uint cds_start; //if TLF: CDS=start:end attribute
    uint cds_end;
    GVec<GSeg> exons; //if TLF: exons= attribute
    GVec<GSeg> cdss; //if TLF: CDS=segment_list attribute
    char* gene_name; //value of gene_name attribute (GTF) if present or Name attribute of a gene feature (GFF3)
    char* gene_id; //GTF only: value of "gene_id" attribute if present
    char** parents; //for GTF only parents[0] is used
    int num_parents;
    char* ID;     // if a ID=.. attribute was parsed, or a GTF with 'transcript' line (transcript_id)
    GffLine(GffReader* reader, const char* l); //parse the line accordingly
    void discardParent() {
    	GFREE(_parents);
    	_parents_len=0;
    	num_parents=0;
    	GFREE(parents);
    	parents=NULL;
    }
    static char* extractGFFAttr(char*& infostr, const char* oline, const char* pre, bool caseStrict=false,
    		bool enforce_GTF2=false, int* rlen=NULL, bool deleteAttr=true);
    char* extractAttr(const char* pre, bool caseStrict=false, bool enforce_GTF2=false, int* rlen=NULL){
    	return extractGFFAttr(info, dupline, pre, caseStrict, enforce_GTF2, rlen, true);
    }
    char* getAttrValue(const char* pre, bool caseStrict=false, bool enforce_GTF2=false, int* rlen=NULL) {
    	return extractGFFAttr(info, dupline, pre, caseStrict, enforce_GTF2, rlen, false);
    }
    GffLine(GffLine& l): _parents(NULL), _parents_len(l._parents_len),
    		dupline(NULL), line(NULL), llen(l.llen), gseqname(NULL), track(NULL),
    		ftype(NULL), ftype_id(l.ftype_id), info(NULL), fstart(l.fstart), fend(l.fend),
			//qstart(l.fstart), qend(l.fend), qlen(l.qlen),
			score(l.score), score_decimals(l.score_decimals), strand(l.strand), flags(l.flags), exontype(l.exontype),
			phase(l.phase), cds_start(l.cds_start), cds_end(l.cds_end), exons(l.exons), cdss(l.cdss),
			gene_name(NULL), gene_id(NULL), parents(NULL), num_parents(l.num_parents), ID(NULL) {
    	//if (l==NULL || l->line==NULL)
    	//	GError("Error: invalid GffLine(l)\n");
    	//memcpy((void*)this, (void*)l, sizeof(GffLine));
    	GMALLOC(line, llen+1);
    	memcpy(line, l.line, llen+1);
    	GMALLOC(dupline, llen+1);
    	memcpy(dupline, l.dupline, llen+1);
    	//--offsets within line[]
    	gseqname=line+(l.gseqname-l.line);
    	track=line+(l.track-l.line);
    	ftype=line+(l.ftype-l.line);
    	info=line+(l.info-l.line);
    	if (num_parents>0 && parents) {
    		GMALLOC(parents, num_parents*sizeof(char*));
    		//_parents_len=l->_parents_len; copied above
    		_parents=NULL; //re-init, forget pointer copy
    		GMALLOC(_parents, _parents_len);
    		memcpy(_parents, l._parents, _parents_len);
    		for (int i=0;i<num_parents;i++) {
    			parents[i]=_parents+(l.parents[i] - l._parents);
    		}
    	}
    	//-- allocated string copies:
    	ID=Gstrdup(l.ID);
    	if (l.gene_name!=NULL)
    		gene_name=Gstrdup(l.gene_name);
    	if (l.gene_id!=NULL)
    		gene_id=Gstrdup(l.gene_id);
    }
    GffLine(): _parents(NULL), _parents_len(0),
    		dupline(NULL), line(NULL), llen(0), gseqname(NULL), track(NULL),
    		ftype(NULL), ftype_id(-1), info(NULL), fstart(0), fend(0), //qstart(0), qend(0), qlen(0),
    		score(0), score_decimals(-1), strand(0), flags(0), exontype(0), phase(0), cds_start(0), cds_end(0),
			exons(), cdss(),  gene_name(NULL), gene_id(NULL), parents(NULL), num_parents(0), ID(NULL) {
    }
    ~GffLine() {
    	GFREE(dupline);
    	GFREE(line);
    	GFREE(_parents);
    	GFREE(parents);
    	GFREE(ID);
    	GFREE(gene_name);
    	GFREE(gene_id);
    }
};

class GffAttr {
 public:
  union {
    int id_full;
    struct {
      bool      cds:1;
      int  attr_id:31;
    };
  };
  char* attr_val;
  GffAttr(int an_id, const char* av=NULL, bool is_cds=false):id_full(0), attr_val(NULL) {
	 attr_id=an_id;
     setValue(av, is_cds);
  }
  ~GffAttr() {
     GFREE(attr_val);
  }
  void setValue(const char* av, bool is_cds=false) {
     if (attr_val!=NULL) {
        GFREE(attr_val);
     }
     if (av==NULL || av[0]==0) return;
     //trim spaces
     const char* vstart=av;
     while (*vstart==' ') av++;
     const char* vend=vstart;
     bool keep_dq=false;
     while (vend[1]!=0) {
        if (*vend==' ' && vend[1]!=' ') keep_dq=true;
          else if (*vend==';') keep_dq=true;
        vend++;
     }
     //remove spaces at the end:
     while (*vend==' ' && vend!=vstart) vend--;
     //practical clean-up: if it doesn't have any internal spaces just strip those useless double quotes
     if (!keep_dq && *vstart=='"' && *vend=='"') {
               vend--;
               vstart++;
     }
     attr_val=Gstrdup(vstart, vend);
     cds=is_cds;
  }
  bool operator==(GffAttr& d){
      return (this==&d);
  }
  bool operator>(GffAttr& d){
     return (this>&d);
  }
  bool operator<(GffAttr& d){
    return (this<&d);
  }

 };

class GffNameList;
class GffNames;

class GffNameInfo {
  friend class GffNameList;
 public:
   int idx;
   char* name;
   GffNameInfo(const char* n=NULL):idx(-1),name(NULL) {
     if (n) name=Gstrdup(n);
     }

   ~GffNameInfo() {
      GFREE(name);
     }

   bool operator==(GffNameInfo& d){
       return (strcmp(this->name, d.name)==0);
       }
   bool operator<(GffNameInfo& d){
     return (strcmp(this->name, d.name)<0);
     }
};

class GffNameList:public GPVec<GffNameInfo> {
  friend class GffNameInfo;
  friend class GffNames;
protected:
  GHash<GffNameInfo> byName;//hash with shared keys
  int idlast; //fList index of last added/reused name
  int addStatic(const char* tname) {// fast add
     GffNameInfo* f=new GffNameInfo(tname);
     idlast=this->Add(f);
     f->idx=idlast;
     byName.shkAdd(f->name,f);
     return idlast;
     }
public:
 //GffNameList(int init_capacity=6):GList<GffNameInfo>(init_capacity, false,true,true), byName(false) {
  GffNameList(int init_capacity=6):GPVec<GffNameInfo>(init_capacity, true), byName(false) {
    idlast=-1;
    setCapacity(init_capacity);
    }
 char* lastNameUsed() { return idlast<0 ? NULL : Get(idlast)->name; }
 int lastNameId() { return idlast; }
 char* getName(int nid) { //retrieve name by its ID
   if (nid<0 || nid>=fCount)
         GError("GffNameList Error: invalid index (%d)\n",nid);
   return fList[nid]->name;
   }

 int addName(const char* tname) {//returns or create an id for the given name
   //check idlast first, chances are it's the same feature name checked
   /*if (idlast>=0 && strcmp(fList[idlast]->name,tname)==0)
       return idlast;*/
   GffNameInfo* f=byName.Find(tname);
   int fidx=-1;
   if (f!=NULL) fidx=f->idx;
     else {//add new entry
      f=new GffNameInfo(tname);
      fidx=this->Add(f);
      f->idx=fidx;
      byName.shkAdd(f->name,f);
      }
   idlast=fidx;
   return fidx;
   }

 int addNewName(const char* tname) {
    GffNameInfo* f=new GffNameInfo(tname);
    int fidx=this->Add(f);
    f->idx=fidx;
    byName.shkAdd(f->name,f);
    return fidx;
    }

 int getId(const char* tname) { //only returns a name id# if found
    GffNameInfo* f=byName.Find(tname);
    if (f==NULL) return -1;
    return f->idx;
    }
 int removeName() {
   GError("Error: removing names from GffNameList not allowed!\n");
   return -1;
   }
};

class GffNames {
 public:
   int numrefs;
   GffNameList tracks;
   GffNameList gseqs;
   GffNameList attrs;
   GffNameList feats; //feature names: 'mRNA', 'exon', 'CDS' etc.
   GffNames():tracks(),gseqs(),attrs(), feats() {
    numrefs=0;
    //the order below is critical!
    //has to match: gff_fid_mRNA, gff_fid_exon
    gff_fid_mRNA = feats.addStatic("mRNA");//index 0=gff_fid_mRNA
    gff_fid_transcript=feats.addStatic("transcript");//index 1=gff_fid_transcript
    gff_fid_exon=feats.addStatic("exon");//index 1=gff_fid_exon
    //feats.addStatic("CDS"); //index 2=gff_fid_CDS
    }
};

void gffnames_ref(GffNames* &n);
void gffnames_unref(GffNames* &n);

enum GffPrintMode {
  pgtfAny, //print record as read, if GTF
  pgtfExon, //print exon only features
  pgtfCDS,  //print CDS and exon features
  pgffAny, //print record as read (if isCDSonly() prints only CDS)
  pgffExon,
  pgffCDS,
  pgffBoth, //enforce exon printing if isCDSOnly()
  pgffTLF,  //exon and CDS data shown as additional GFF attributes
            //in the transcript line (Transcript Line Format)
            //every line has the whole transcript data
  pgffBED //print a BED line with all other GFF attributes in column 13
};


class GffAttrs:public GList<GffAttr> {
  public:
    GffAttrs():GList<GffAttr>(false,true,true) { }
    void add_if_new(GffNames* names, const char* attrname, const char* attrval) {
        //adding a new value without checking for cds status
        int nid=names->attrs.getId(attrname);
        if (nid>=0) { //attribute name found in the dictionary
           for (int i=0;i<Count();i++)
              if (nid==Get(i)->attr_id) { return; } //don't update existing
        }
        else { //adding attribute name to global attr name dictionary
           nid=names->attrs.addNewName(attrname);
        }
        this->Add(new GffAttr(nid, attrval));
    }

    void add_if_new(GffNames* names, const char* attrname, const char* attrval, bool is_cds) {
        int nid=names->attrs.getId(attrname);
        if (nid>=0) { //attribute name found in the dictionary
           for (int i=0;i<Count();i++)
              if (nid==Get(i)->attr_id && is_cds==Get(i)->cds) { return; } //don't update existing
        }
        else { //adding attribute name to global attr name dictionary
           nid=names->attrs.addNewName(attrname);
        }
        this->Add(new GffAttr(nid, attrval, is_cds));
    }

    void add_or_update(GffNames* names, const char* attrname, const char* val) {
    //adding a new value without checking for cds status
        int aid=names->attrs.getId(attrname);
        if (aid>=0) {
           //attribute found in the dictionary
           for (int i=0;i<Count();i++) {
              //do we have it?
              if (aid==Get(i)->attr_id) {
                  //update the existing value for this attribute
                  Get(i)->setValue(val);
                  return;
                  }
              }
        }
        else { //adding attribute name to global attr name dictionary
           aid=names->attrs.addNewName(attrname);
        }
        this->Add(new GffAttr(aid, val));
    }

    void add_or_update(GffNames* names, const char* attrname, const char* val, bool is_cds) {
      int aid=names->attrs.getId(attrname);
      if (aid>=0) {
         //attribute found in the dictionary
         for (int i=0;i<Count();i++) {
            //do we have it?
            if (aid==Get(i)->attr_id && Get(i)->cds==is_cds) {
                //update the existing value for this attribute
                Get(i)->setValue(val, is_cds);
                return;
                }
            }
      }
      else { //adding attribute name to global attr name dictionary
         aid=names->attrs.addNewName(attrname);
      }
      this->Add(new GffAttr(aid, val, is_cds));
    }

    int haveId(int attr_id, bool is_cds=false) {
        for (int i=0;i<Count();i++)
           if (attr_id==Get(i)->attr_id && Get(i)->cds==is_cds)
        	   return i;
        return -1;
    }

    int haveId(const char* attrname, GffNames* names, bool is_cds=false) {
    	int aid=names->attrs.getId(attrname);
    	if (aid>=0) {
            for (int i=0;i<Count();i++)
               if (aid==Get(i)->attr_id && Get(i)->cds==is_cds)
            	   return i;
    	}
    	return -1;
    }

    char* getAttr(GffNames* names, const char* attrname) {
      int aid=names->attrs.getId(attrname);
      if (aid>=0)
        for (int i=0;i<Count();i++)
          if (aid==Get(i)->attr_id) return Get(i)->attr_val;
      return NULL;
    }

    char* getAttr(GffNames* names, const char* attrname, bool is_cds) {
      int aid=names->attrs.getId(attrname);
      if (aid>=0)
        for (int i=0;i<Count();i++)
          if (aid==Get(i)->attr_id && Get(i)->cds==is_cds) return Get(i)->attr_val;
      return NULL;
    }

    char* getAttr(int aid) {
      if (aid>=0)
        for (int i=0;i<Count();i++)
          if (aid==Get(i)->attr_id) return Get(i)->attr_val;
      return NULL;
    }

    char* getAttr(int aid, bool is_cds) {
      if (aid>=0)
        for (int i=0;i<Count();i++)
          if (aid==Get(i)->attr_id && Get(i)->cds==is_cds)
            return Get(i)->attr_val;
      return NULL;
    }

    void copyAttrs(GffAttrs* attrs, bool is_cds=false) {
 	 //deep copy attributes from another GffAttrs list
     // (only the ones which do not exist yet)
    	if (attrs==NULL) return;
    	for (int i=0;i<attrs->Count();i++) {
    		int aid=attrs->Get(i)->attr_id;
    		if (haveId(aid, is_cds)<0)
    			Add(new GffAttr(aid, attrs->Get(i)->attr_val, is_cds));
    	}
    }
};

class GffExon : public GSeg {
 public:
  bool sharedAttrs; //do not free attrs on destruct!
  GffAttrs* attrs; //other attributes kept for this exon/CDS
  GffScore score; // gff score column
  int8_t exontype;
  char phase; //GFF phase column - for CDS segments only!
	             // '.' = undefined (UTR), '0','1','2' for CDS exons
  void* uptr; //for associating extended user data to this exon

  char* getAttr(GffNames* names, const char* atrname) {
    if (attrs==NULL || names==NULL || atrname==NULL) return NULL;
    return attrs->getAttr(names, atrname);
  }

  char* getAttr(int aid) {
    if (attrs==NULL) return NULL;
    return attrs->getAttr(aid);
  }
  GffExon(bool share_attributes):GSeg(0,0), sharedAttrs(share_attributes), attrs(NULL), score(),
		  exontype(0), phase('.'), uptr(NULL){
  }
  GffExon(uint s=0, uint e=0, int8_t et=0, char ph='.', float sc=0, int8_t sc_prec=0):sharedAttrs(false), attrs(NULL),
		  score(sc,sc_prec), exontype(et), phase(ph), uptr(NULL) {
		if (s<e) { start=s; end=e; }
		    else { start=e; end=s; }

  } //constructor


  GffExon(const GffExon& ex):GSeg(ex.start, ex.end) { //copy constructor
      (*this)=ex; //use the default (shallow!) copy operator
      if (ex.attrs!=NULL) { //make a deep copy here
        attrs=new GffAttrs();
        attrs->copyAttrs(ex.attrs);
      }
  }

  GffExon& operator=(const GffExon& o) = default; //prevent gcc 9 warnings:
                                           //yes, I want a shallow copy here

  ~GffExon() { //destructor
     if (attrs!=NULL && !sharedAttrs) delete attrs;
  }

};

//only for mapping to spliced coding sequence:
class GffCDSeg:public GSeg {
 public:
  char phase;
  int exonidx;
};

//one GFF mRNA object -- e.g. a mRNA with its exons and/or CDS segments
class GffObj:public GSeg {
 protected:
   char* gffID; // ID name for mRNA (parent) feature
   char* gene_name; //value of gene_name attribute (GTF) if present or Name attribute of the parent gene feature (GFF3)
   char* geneID; //value of gene_id attribute (GTF) if present, or the ID attribute of a parent gene feature (GFF3)
   union {
      unsigned int flags;
      struct {
    	  bool flag_HAS_ERRORS        :1;
    	  bool flag_CHILDREN_PROMOTED :1;
    	  bool flag_IS_GENE           :1;
    	  bool flag_IS_TRANSCRIPT     :1;
    	  bool flag_HAS_GFF_ID        :1; //found transcript/RNA feature line (GFF3 or GTF2 with transcript line)
    	  bool flag_BY_EXON           :1; //created by subfeature (exon/CDS) directly
    	  bool flag_CDS_ONLY          :1; //transcript defined by CDS features only (GffObj::isCDS())
    	  bool flag_CDS_NOSTART       :1; //partial CDS at 5' end (no start codon)
    	  bool flag_CDS_NOSTOP        :1; //partial CDS at 3' end (no stop codon)
    	  bool flag_CDS_X             :1; //transcript having CDS with ribosomal shift (i.e. after merging exons)
    	                                  //CDS segments stored in ::cdss are incompatible with the exon segments
    	  bool flag_GENE_SEGMENT      :1; //a transcript-like C/D/J/V_gene_segment (NCBI's annotation)
    	  bool flag_TRANS_SPLICED     :1;
    	  bool flag_DISCONTINUOUS     :1; //discontinuous feature (e.g. cDNA_match) segments linked by same ID
    	  bool flag_TARGET_ONLY       :1; //Target= feature (e.g. from RepeatMasker output), lacks ID
    	  bool flag_DISCARDED         :1; //it will be  discarded from the final GffReader list
    	  bool flag_LST_KEEP          :1; //controlled by isUsed(); if set, this GffObj will not be
    	                                  //deallocated when GffReader is destroyed
    	  bool flag_FINALIZED         :1; //if finalize() was already called for this GffObj
    	  unsigned int gff_level      :4; //hierarchical level (0..15)
      };
   };
   //-- friends:
   friend class GffReader;
   friend class GffExon;
public:
  static GffNames* names; // dictionary storage that holds the various attribute names etc.
  int track_id; // index of track name in names->tracks
  int gseq_id; // index of genomic sequence name in names->gseqs
  int ftype_id; // index of this record's feature name in names->feats, or the special gff_fid_mRNA value
  int subftype_id; //index of child subfeature name in names->feats (subfeatures stored in "exons")
                   //if ftype_id==gff_fid_mRNA then this value is ignored
  GList<GffExon> exons; //for non-mRNA entries, these can be any subfeature of type subftype_id
  GList<GffExon>* cdss; //only !NULL for cases of "programmed frameshift" when CDS boundaries do not match
                      //exons boundaries
  GPVec<GffObj> children;
  GffObj* parent;
  int udata; //user data, flags etc.
  void* uptr; //user pointer (to a parent object, cluster, locus etc.)
  GffObj* ulink; //link to another GffObj (user controlled field)
  //---mRNA specific fields:
  //bool isCDS; //just a CDS, no UTRs
  uint CDstart; //CDS lowest coordinate
  uint CDend;   //CDS highest coordinate
  char CDphase; //initial phase for CDS start ('.','0'..'2')
                //CDphase is at CDend if strand=='-'
  static void decodeHexChars(char* dbuf, const char* s, int maxlen=1023);
  bool hasErrors() { return flag_HAS_ERRORS; }
  void hasErrors(bool v) { flag_HAS_ERRORS=v; }
  bool hasGffID() { return flag_HAS_GFF_ID; }
  void hasGffID(bool v) {flag_HAS_GFF_ID=v; }
  bool createdByExon() { return flag_BY_EXON; }
  void createdByExon(bool v) {flag_BY_EXON=v; }
  bool isCDSOnly() { return flag_CDS_ONLY; }
  void isCDSOnly(bool v) {  flag_CDS_ONLY=v; }
  bool isXCDS() { return flag_CDS_X; }
  void isXCDS(bool v) {  flag_CDS_X=v; }
  bool isFinalized() {  return flag_FINALIZED; }
  void isFinalized(bool v) {  flag_FINALIZED=v; }

  bool isGene() { return flag_IS_GENE; }
  void isGene(bool v) {flag_IS_GENE=v; }
  bool isDiscarded() { return flag_DISCARDED; }
  void isDiscarded(bool v) { flag_DISCARDED=v; }
  bool isUsed() { return flag_LST_KEEP; }
  void isUsed(bool v) {flag_LST_KEEP=v; }
  bool isTranscript() { return flag_IS_TRANSCRIPT; }
  void isTranscript(bool v) {flag_IS_TRANSCRIPT=v; }
  bool isGeneSegment() { return flag_GENE_SEGMENT; }
  void isGeneSegment(bool v) {flag_GENE_SEGMENT=v; }
  bool promotedChildren() { return flag_CHILDREN_PROMOTED; }
  void promotedChildren(bool v) { flag_CHILDREN_PROMOTED=v; }
  void setLevel(byte v) { gff_level=v; }
  byte getLevel() { return gff_level; }
  byte incLevel() { gff_level++; return gff_level; }

  bool isValidTranscript() {
    //return (ftype_id==gff_fid_mRNA && exons.Count()>0);
    return (isTranscript() && exons.Count()>0);
  }

  //return the index of exon containing coordinate coord, or -1 if not
  int whichExon(uint coord, GList<GffExon>* segs=NULL);
  int readExon(GffReader& reader, GffLine& gl);

  int addExon(GList<GffExon>& segs, GffLine& gl, int8_t exontype_override=exgffNone); //add to cdss or exons

  int addExon(uint segstart, uint segend, int8_t exontype=exgffNone, char phase='.',
		      GffScore exon_score=GFFSCORE_NONE, GList<GffExon>* segs=NULL);

protected:
  bool reduceExonAttrs(GList<GffExon>& segs);
  //utility segment-merging function for addExon()
  void expandSegment(GList<GffExon>&segs, int oi, uint segstart, uint segend,
       int8_t exontype);
  bool processGeneSegments(GffReader* gfr); //for genes that have _gene_segment features (NCBI annotation)
  void transferCDS(GffExon* cds);
public:
  void removeExon(int idx);
  void removeExon(GffExon* p);
  char  strand; //true if features are on the reverse complement strand
  GffScore gscore;
  int covlen; //total coverage of reference genomic sequence (sum of maxcf segment lengths)
  GffAttrs* attrs; //other gff3 attributes found for the main mRNA feature
   //constructor by gff line parsing:
  GffObj(GffReader& gfrd, BEDLine& bedline);
  GffObj(GffReader& gfrd, GffLine& gffline);
   //if gfline->Parent!=NULL then this will also add the first sub-feature
   // otherwise, only the main feature is created
  void copyAttrs(GffObj* from);
  void clearAttrs() {
    if (attrs!=NULL) {
      bool sharedattrs=(exons.Count()>0 && exons[0]->attrs==attrs);
      delete attrs; attrs=NULL;
      if (sharedattrs) exons[0]->attrs=NULL;
      }
    }
  GffObj(char* anid=NULL):GSeg(0,0), exons(true,true,false), cdss(NULL), children(1,false), gscore() {
                                   //exons: sorted, free, non-unique
       gffID=NULL;
       uptr=NULL;
       ulink=NULL;
       flags=0;
       udata=0;
       parent=NULL;
       ftype_id=-1;
       subftype_id=-1;
       if (anid!=NULL) gffID=Gstrdup(anid);
       gffnames_ref(names);
       CDstart=0; // hasCDS <=> CDstart>0
       CDend=0;
       CDphase=0;
       gseq_id=-1;
       track_id=-1;
       strand='.';
       attrs=NULL;
       covlen=0;
       geneID=NULL;
       gene_name=NULL;
   }
   ~GffObj() {
       GFREE(gffID);
       GFREE(gene_name);
       GFREE(geneID);
       delete cdss;
       clearAttrs();
       gffnames_unref(names);
       }
   //--------------
   GffObj* finalize(GffReader* gfr);
               //complete parsing: must be called in order to merge adjacent/close proximity subfeatures
   void parseAttrs(GffAttrs*& atrlist, char* info, bool isExon=false, bool CDSsrc=false);
   const char* getSubfName() { //returns the generic feature type of the entries in exons array
     //int sid=exon_ftype_id;
     //if (sid==gff_fid_exon && isCDS) sid=gff_fid_CDS;
     return names->feats.getName(subftype_id);
     }
   void setCDS(uint cd_start, uint cd_end, char phase=0);
   void setCDS(GffObj* t); //set CDS from another transcript

   bool monoFeature() {
     return (exons.Count()==0 ||
          (exons.Count()==1 &&  //exon_ftype_id==ftype_id &&
              exons[0]->end==this->end && exons[0]->start==this->start));
     }

   bool hasCDS() { return (CDstart>0); }

   const char* getFeatureName() {
     return names->feats.getName(ftype_id);
     }
   void setFeatureName(const char* feature);

   void addAttr(const char* attrname, const char* attrvalue);
   int removeAttr(const char* attrname, const char* attrval=NULL);
   int removeAttr(int aid, const char* attrval=NULL);
   int removeExonAttr(GffExon& exon, const char* attrname, const char* attrval=NULL);
   int removeExonAttr(GffExon& exon, int aid, const char* attrval=NULL);

   const char* getAttrName(int i) {
     if (attrs==NULL) return NULL;
     return names->attrs.getName(attrs->Get(i)->attr_id);
   }

   char* getAttr(const char* attrname, bool checkFirstExon=false) {
     if (names==NULL || attrname==NULL) return NULL;
     char* r=NULL;
     if (attrs==NULL) {
         if (!checkFirstExon) return NULL;
     } else
         r=attrs->getAttr(names, attrname);
     if (r!=NULL) return r;
     if (checkFirstExon && exons.Count()>0) {
        r=exons.First()->getAttr(names, attrname);
     }
     return r;
   }

   char* getExonAttr(GffExon* exon, const char* attrname) {
      if (exon==NULL || attrname==NULL) return NULL;
      return exon->getAttr(names, attrname);
   }

   char* getExonAttr(int exonidx, const char* attrname) {
      if (exonidx<0 || exonidx>=exons.Count() || attrname==NULL) return NULL;
      return exons[exonidx]->getAttr(names, attrname);
      }

   char* getAttrValue(int i) {
     if (attrs==NULL) return NULL;
     return attrs->Get(i)->attr_val;
     }
   const char* getGSeqName() {
     return names->gseqs.getName(gseq_id);
     }

   const char* getRefName() {
     return names->gseqs.getName(gseq_id);
     }
   void setRefName(const char* newname);

   const char* getTrackName() {
     return names->tracks.getName(track_id);
   }
   bool exonOverlap(uint s, uint e) {//check if ANY exon overlaps given segment
      //ignores strand!
      if (s>e) Gswap(s,e);
      for (int i=0;i<exons.Count();i++) {
         if (exons[i]->overlap(s,e)) return true;
         }
      return false;
   }
   bool exonOverlap(GffObj& m) {//check if ANY exon overlaps given segment
     //if (gseq_id!=m.gseq_id) return false;
     // ignores strand and gseq_id, must check in advance
     for (int i=0;i<exons.Count();i++) {
         for (int j=0;j<m.exons.Count();j++) {
            if (exons[i]->start>m.exons[j]->end) continue;
            if (m.exons[j]->start>exons[i]->end) break;
            //-- overlap if we are here:
            return true;
         }
     }
     return false;
   }

   int exonOverlapIdx(GList<GffExon>& segs, uint s, uint e, int* ovlen=NULL, int start_idx=0);

   int exonOverlapLen(GffObj& m) {
      if (start>m.end || m.start>end) return 0;
      int i=0;
      int j=0;
      int ovlen=0;
      while (i<exons.Count() && j<m.exons.Count()) {
        uint istart=exons[i]->start;
        uint iend=exons[i]->end;
        uint jstart=m.exons[j]->start;
        uint jend=m.exons[j]->end;
        if (istart>jend) { j++; continue; }
        if (jstart>iend) { i++; continue; }
        //exon overlap
        uint ovstart=GMAX(istart,jstart);
        if (iend<jend) {
           ovlen+=iend-ovstart+1;
           i++;
           }
        else {
           ovlen+=jend-ovstart+1;
           j++;
           }
        }//while comparing exons
      return ovlen;
      }

    bool exonOverlap(GffObj* m) {
      return exonOverlap(*m);
      }

   //---------------------
   bool operator==(GffObj& d){
       return (gseq_id==d.gseq_id && start==d.start && end==d.end && strcmp(gffID, d.gffID)==0);
       }
   bool operator>(GffObj& d){
      if (gseq_id!=d.gseq_id) return (gseq_id>d.gseq_id);
      if (start==d.start) {
         if (getLevel()==d.getLevel()) {
             if (end==d.end) return (strcmp(gffID, d.gffID)>0);
                        else return (end>d.end);
             } else return (getLevel()>d.getLevel());
         } else return (start>d.start);
      }
   bool operator<(GffObj& d){
     if (gseq_id!=d.gseq_id) return (gseq_id<d.gseq_id);
     if (start==d.start) {
         if (getLevel()==d.getLevel()) {
            if (end==d.end) return strcmp(gffID, d.gffID)<0;
                     else return end<d.end;
            } else return (getLevel()<d.getLevel());
        } else return (start<d.start);
     }
   char* getID() { return gffID; }

   char* getGeneID() { return geneID; }
   char* getGeneName() { return gene_name; }
   void setGeneName(const char* gname) {
        GFREE(gene_name);
        if (gname) gene_name=Gstrdup(gname);
   }
   void setGeneID(const char* gene_id) {
        GFREE(geneID);
        if (gene_id) geneID=Gstrdup(gene_id);
   }
   int addSeg(GffLine* gfline);
   int addSeg(int fnid, GffLine* gfline);
   void getCDSegs(GVec<GffExon>& cds);

   void updateCDSPhase(GList<GffExon>& segs); //for CDS-only features, updates GffExon::phase
   void printGTab(FILE* fout, char** extraAttrs=NULL);
   void printGxfExon(FILE* fout, const char* tlabel, const char* gseqname,
          bool iscds, GffExon* exon, bool gff3, bool cvtChars, char* dbuf, int dbuf_len);
   void printGxf(FILE* fout, GffPrintMode gffp=pgffExon,
             const char* tlabel=NULL, const char* gfparent=NULL, bool cvtChars=false);
   void printGtf(FILE* fout, const char* tlabel=NULL, bool cvtChars=false) {
      printGxf(fout, pgtfAny, tlabel, NULL, cvtChars);
   }
   void printGff(FILE* fout, const char* tlabel=NULL,
                                const char* gfparent=NULL, bool cvtChars=false) {
      printGxf(fout, pgffAny, tlabel, gfparent, cvtChars);
   }
   void printTranscriptGff(FILE* fout, char* tlabel=NULL,
                            bool showCDS=false, const char* gfparent=NULL, bool cvtChars=false) {
      if (isValidTranscript())
         printGxf(fout, showCDS ? pgffBoth : pgffExon, tlabel, gfparent, cvtChars);
      }
   void printExonList(FILE* fout); //print comma delimited list of exon intervals
   void printCDSList(FILE* fout); //print comma delimited list of CDS intervals

   void printBED(FILE* fout, bool cvtChars, char* dbuf, int dbuf_len);
       //print a BED-12 line + GFF3 attributes in 13th field
   void printSummary(FILE* fout=NULL);

   char* getSpliced(GFaSeqGet* faseq, bool CDSonly=false, int* rlen=NULL,
           uint* cds_start=NULL, uint* cds_end=NULL, GMapSegments* seglst=NULL,
		   bool cds_open=false);
    char* getUnspliced(GFaSeqGet* faseq, int* rlen, GMapSegments* seglst=NULL);

    void addPadding(int padLeft, int padRight); //change exons to include this padding on the sides
    void removePadding(int padLeft, int padRight);

   //bool validCDS(GFaSeqGet* faseq); //has In-Frame Stop Codon ?
   bool empty() { return (start==0); }
};

typedef bool GffRecFunc(GffObj* gobj, void* usrptr1, void* usrptr2);
//user callback after parsing a mapping object:
// Returns: "done with it" status:
//   TRUE if gobj is no longer needed so it's FREEd upon return
//   FALSE if the user needs the gobj pointer and is responsible for
//             collecting and freeing all GffObj objects


//GSeqStat: collect basic stats about a common underlying genomic sequence
//          for multiple GffObj
class GSeqStat {
 public:
   int gseqid; //gseq id in the global static pool of gseqs
   char* gseqname; //just a pointer to the name of gseq
   int fcount;//number of features on this gseq
   uint mincoord;
   uint maxcoord;
   uint maxfeat_len; //maximum feature length on this genomic sequence
   GffObj* maxfeat;
   GSeqStat(int id=-1, char* name=NULL) {
     gseqid=id;
     gseqname=name;
     fcount=0;
     mincoord=MAXUINT;
     maxcoord=0;
     maxfeat_len=0;
     maxfeat=NULL;
     }
   bool operator>(GSeqStat& g) {
    return (gseqid>g.gseqid);
    }
   bool operator<(GSeqStat& g) {
    return (gseqid<g.gseqid);
    }
   bool operator==(GSeqStat& g) {
    return (gseqid==g.gseqid);
    }
};


int gfo_cmpByLoc(const pointer p1, const pointer p2);
int gfo_cmpRefByID(const pointer p1, const pointer p2);

class GfList: public GList<GffObj> {
 public:

   GfList(bool sorted):GList<GffObj>(sorted,false,false) { }
   GfList():GList<GffObj>(false,false,false) {
     //GffObjs in this list are NOT deleted when the list is cleared
     //-- for deallocation of these objects, call freeAll() or freeUnused() as needed
   }
   void finalize(GffReader* gfr);

   void freeAll() {
     for (int i=0;i<fCount;i++) {
       delete fList[i];
       fList[i]=NULL;
       }
     Clear();
   }
   void freeUnused() {
     for (int i=0;i<fCount;i++) {
       if (fList[i]->isUsed()) continue;
       /*//inform the children?
       for (int c=0;c<fList[i]->children.Count();c++) {
          fList[i]->children[c]->parent=NULL;
       }
       */
       delete fList[i];
       fList[i]=NULL;
       }
     Clear();
     }
};

class CNonExon { //utility class used in subfeature promotion
 public:
   //int idx;
   GffObj* parent;
   GffExon* exon;
   GffLine* gffline;
   //CNonExon(int i, GffObj* p, GffExon* e, GffLine* gl) {
   CNonExon(GffObj* p, GffExon* e, GffLine& gl) {
     parent=p;
     exon=e;
     //idx=i;
     gffline=new GffLine(gl);
     }
  ~CNonExon() {
     delete gffline;
     }
 };


class GffReader {
  friend class GffObj;
  friend class GffLine;
  friend class GfList;
  char* linebuf;
  off_t fpos;
  int buflen;
 protected:
  union {
	unsigned int flags;
    unsigned int gff_type: 6;
    struct {
       bool is_gff3: 1;  //GFF3 syntax was detected
       bool is_gtf:1; //GTF syntax was detected
       bool gtf_transcript:1; //has "transcript" features (2-level GTF)
       bool gtf_gene:1; //has "gene" features (3-level GTF ..Ensembl?)
       bool is_BED:1; //input is BED-12 format, possibly with attributes in 13th field
       bool is_TLF:1; //input is GFF3-like Transcript Line Format with exons= attribute
       //--other flags
       bool transcripts_Only:1; //default ; only keep recognized transcript features
       bool keep_Genes:1; //for transcriptsOnly, do not discard genes from gflst
       bool keep_Attrs:1;
       bool keep_AllExonAttrs:1; //when keep_Attrs, do not attempt to reduce exon attributes
       bool noExonAttrs:1;
       bool ignoreLocus:1; //discard locus features and attributes from input
       bool merge_CloseExons:1;
       bool gene2exon:1;
       bool sortByLoc:1; //if records should be sorted by location
       bool refAlphaSort:1; //if sortByLoc, reference sequences are
                       // sorted lexically instead of their id#
       bool gff_warns:1;
    };
  };
  //char* lastReadNext;
  FILE* fh;
  char* fname;  //optional fasta file with the underlying genomic sequence to be attached to this reader
  GFFCommentParser* commentParser;
  GffLine* gffline;
  BEDLine* bedline;
  //bool transcriptsOnly; //keep only transcripts w/ their exon/CDS features
  //bool gene2exon;  // for childless genes: add an exon as the entire gene span
  GHash<int> discarded_ids; //for transcriptsOnly mode, keep track
                            // of discarded parent IDs
  GHash< GPVec<GffObj> > phash; //transcript_id => GPVec<GffObj>(false)
  //GHash<int> tids; //just for transcript_id uniqueness
  char* gfoBuildId(const char* id, const char* ctg);
  //void gfoRemove(const char* id, const char* ctg);
  GffObj* gfoAdd(GffObj* gfo);
  GffObj* gfoAdd(GPVec<GffObj>& glst, GffObj* gfo);
  GffObj* gfoReplace(GPVec<GffObj>& glst, GffObj* gfo, GffObj* toreplace);
  // const char* id, const char* ctg, char strand, GVec<GfoHolder>** glst, uint start, uint end
  bool pFind(const char* id, GPVec<GffObj>*& glst);
  GffObj* gfoFind(const char* id, GPVec<GffObj>* & glst, const char* ctg=NULL,
	                                         char strand=0, uint start=0, uint end=0);
  CNonExon* subfPoolCheck(GffLine* gffline, GHash<CNonExon>& pex, char*& subp_name);
  void subfPoolAdd(GHash<CNonExon>& pex, GffObj* newgfo);
  GffObj* promoteFeature(CNonExon* subp, char*& subp_name, GHash<CNonExon>& pex);

#ifdef CUFFLINKS
     boost::crc_32_type  _crc_result;
#endif
 public:
  GPVec<GSeqStat> gseqtable; //table with all genomic sequences, but only current GXF gseq ID indices will have non-NULL
  //GffNames* names; //just a pointer to the global static Gff names repository
  GfList gflst; //keeps track of all GffObj records being read (when readAll() is used)
  GffObj* newGffRec(GffLine* gffline, GffObj* parent=NULL, GffExon* pexon=NULL,
		       GPVec<GffObj>* glst=NULL, bool replace_parent=false);
  GffObj* newGffRec(BEDLine* bedline, GPVec<GffObj>* glst=NULL);
  //GffObj* replaceGffRec(GffLine* gffline, bool keepAttr, bool noExonAttr, int replaceidx);
  GffObj* updateGffRec(GffObj* prevgfo, GffLine* gffline);
  GffObj* updateParent(GffObj* newgfh, GffObj* parent);
  bool readExonFeature(GffObj* prevgfo, GffLine* gffline, GHash<CNonExon>* pex=NULL);
  GPVec<GSeqStat> gseqStats; //populated after finalize() with only the ref seqs in this file
  GffReader(FILE* f=NULL, bool t_only=false, bool sort=false):linebuf(NULL), fpos(0),
		  buflen(0), flags(0), fh(f), fname(NULL), commentParser(NULL), gffline(NULL),
		  bedline(NULL), discarded_ids(true), phash(true), gseqtable(1,true),
		  gflst(), gseqStats(1, false) {
      GMALLOC(linebuf, GFF_LINELEN);
      buflen=GFF_LINELEN-1;
      gffnames_ref(GffObj::names);
      //gff_warns=gff_show_warnings;
      transcripts_Only=t_only;
      sortByLoc=sort;
      noExonAttrs=true;
      //lastReadNext=NULL;
  }
  /*
  void init(FILE *f, bool t_only=false, bool sortbyloc=false, bool g2exon=false) {
      fname=NULL;
      fh=f;
      if (fh!=NULL) rewind(fh);
      fpos=0;
      flags=0;
      transcriptsOnly=t_only;
      gflst.sortedByLoc(sortbyloc);
      gene2exon=g2exon;
  }
  */
  void gene2Exon(bool v) { gene2exon=v;}
  void enableSorting(bool sorting=true) { sortByLoc=sorting; }
  bool getSorting() { return sortByLoc; }
  void isBED(bool v=true) { is_BED=v; } //should be set before any parsing!
  void isTLF(bool v=true) { is_TLF=v; } //should be set before any parsing!
  void keepAttrs(bool keep_attrs=true, bool discardExonAttrs=true, bool preserve_exon_attrs=false) {
	  keep_Attrs=keep_attrs;
	  noExonAttrs=discardExonAttrs;
	  keep_AllExonAttrs=preserve_exon_attrs;
  }
  void transcriptsOnly(bool t_only) { transcripts_Only=t_only; }
  bool transcriptsOnly() { return transcripts_Only; }
  void setIgnoreLocus(bool nolocus) { ignoreLocus=nolocus; }
  void keepGenes(bool keep_genes) {
	  keep_Genes=keep_genes;
  }
  bool keepGenes() { return keep_Genes; }
  void mergeCloseExons(bool merge_close_exons=true) {
	  merge_CloseExons=merge_close_exons;
  }
  void showWarnings(bool v) {
     gff_warns=v;
     //gff_show_warnings=v;
  }
  bool showWarnings() {
     return gff_warns;
  }
  void setRefAlphaSorted(bool v=true) {
	refAlphaSort=v;
	if (v) sortByLoc=true;
  }
  void setCommentParser(GFFCommentParser* cmParser=NULL) {
	  commentParser=cmParser;
  }

  GffReader(const char* fn, bool t_only=false, bool sort=false):linebuf(NULL), fpos(0),
	  		  buflen(0), flags(0), fh(NULL), fname(NULL), commentParser(NULL),
			  gffline(NULL), bedline(NULL), discarded_ids(true),
			  phash(true), gseqtable(1,true), gflst(), gseqStats(1,false) {
      //gff_warns=gff_show_warnings;
      gffnames_ref(GffObj::names);
      noExonAttrs=true;
      transcripts_Only=t_only;
      sortByLoc=sort;
      fname=Gstrdup(fn);
      fh=fopen(fname, "rb");
      GMALLOC(linebuf, GFF_LINELEN);
      buflen=GFF_LINELEN-1;
      //lastReadNext=NULL;
      }

 ~GffReader() {
      delete gffline;
      gffline=NULL;
      fpos=0;
      if (fh && fh!=stdin) fclose(fh);
      gflst.freeUnused();
      gflst.Clear();
      discarded_ids.Clear();
      phash.Clear();
      GFREE(fname);
      GFREE(linebuf);
      //GFREE(lastReadNext);
      gffnames_unref(GffObj::names);
      }


  GffLine* nextGffLine();
  BEDLine* nextBEDLine();

  // load all subfeatures, re-group them:
  void readAll();
  void readAll(bool keepAttr, bool mergeCloseExons=false, bool noExonAttr=true) {
	  this->keep_Attrs=keepAttr;
	  this->merge_CloseExons=mergeCloseExons;
	  this->noExonAttrs=noExonAttr;
	  readAll();
  }


  //only for well-formed files: BED or GxF where exons are strictly grouped by their transcript_id/Parent
  GffObj* readNext(); //user must free the returned GffObj* !

#ifdef CUFFLINKS
    boost::crc_32_type current_crc_result() const { return _crc_result; }
#endif

}; // end of GffReader

// ----------------------------------------------------------
// -- auxiliary classes for GffObj::processGeneSegments() --
class GSegMatch { //keep track of "matching" overlaps of a GeneCDSChain with multiple GeneSegment containers
  public:
   int child_idx; //index of matching _gene_segment GffObj in gene->children[] list
   int noncov; //number of "non-covered" bases in the GeneSegment
   int gsegidx; //index of _gene_segment in GVec<int> geneSegs
     // (i.e. UTRs + implied introns if exons are missing)
   bool operator<(GSegMatch& o) { return (noncov<o.noncov); }
   bool operator==(GSegMatch& o) { return (noncov==o.noncov); }
   GSegMatch(int cidx=-1, int ncov=-1, int gsidx=-1):child_idx(cidx),
    	   noncov(ncov), gsegidx(gsidx) { }
};

class GeneCDS: public GSeg {
  public:
	int idx; //index of this CDS entry in this gene->cdss[] list
	GeneCDS(int i=-1, uint cstart=0, uint cend=0):GSeg(cstart, cend), idx(i) {
	}
};

class GeneCDSChain: public GSeg { //keep track of CDS chains of the gene and their boundaries
  public:
	GVec<GeneCDS> cdsList; //all CDSs in this chain
	GArray<GSegMatch> mxs; //list of "matching" container X_gene_segment transcripts;
	GeneCDSChain():cdsList(),mxs() { }
	GeneCDSChain(int idx, uint cstart, uint cend):GSeg(cstart, cend),
	    	cdsList(),mxs(true) {
	    addCDS(idx, cstart, cend);

	}
	void addCDS(int idx, uint cstart, uint cend) {
	    GeneCDS cds(idx, cstart, cend);
	    cdsList.Add(cds);
	    expandInclude(cstart, cend);
	}
	void addMatch(int childidx, int ncov, int gsegidx) {
	    GSegMatch segmatch(childidx, ncov, gsegidx);
	    mxs.Add(segmatch);
	}
	bool singleExonCDSMatch(uint tstart, uint tend, int& ncov) {
	    if (start>=tstart && end<=tend) {
	    	ncov=start-tstart + tend-end;
	    	//add all CDS-"introns"
	    	if (cdsList.Count()>1)
	    		//shouldn't really consider this a valid "match"
	    		for (int i=1;i<cdsList.Count();i++)
	    			ncov+=cdsList[i].start-cdsList[i-1].end-1;
	    	return true;
	    }
	    return false;
	}
	bool singleCDStoExon(GffObj&t, int& ncov) {
	    //cdsList[0] must be contained in a single exon of t
	    int nc=0;
	    bool match=false;
	    for (int i=0;i<t.exons.Count();i++) {
	    	if (t.exons[i]->overlap(cdsList[0])) {
	    	   if (cdsList[0].start>=t.exons[i]->start &&
	    			cdsList[0].end<=t.exons[i]->end) {
	    		  match=true;
	    		  nc+=cdsList[0].start-t.exons[i]->start+t.exons[i]->end+cdsList[0].end;
	    	   } //contained in this exon
	    	   else return false; //overlap, but not contained
	    	   continue;
	    	}
	    	nc+=t.exons[i]->len();
	    }
	    if (!match) return false;
	    ncov=nc;
	    return true;
	}

	bool multiCDStoExon(GffObj &t, int& ncov) {
	    //multi-CDS vs multi-exon t
	    int nc=0;
	    int e=0, c=0;
	    int emax=t.exons.Count()-1;
	    int cmax=cdsList.Count()-1;
	    int mintrons=0; //matched introns
	    while (e<emax && c<cmax) {
	    	if (mintrons>0 &&
	    			(cdsList[c].end!=t.exons[e]->end ||
	    					cdsList[c+1].start!=t.exons[e+1]->start))
	    		return false;
	    	GSeg cintron(cdsList[c].end+1, cdsList[c+1].start-1);
	    	GSeg eintron(t.exons[e]->end+1, t.exons[e+1]->start-1);
	    	if (cintron.start>eintron.end) {
	    		nc+=t.exons[e]->len();
	    		e++;
	    		continue;
	    	}
	    	if (eintron.start<=cintron.end) {
	    		//intron overlap
	    		if (cintron.start==eintron.start &&
	    				cintron.end==eintron.end) {
	    			//intron match
	    			if (mintrons==0) {
	    				if (cdsList[c].start<t.exons[e]->start) return false;
	    				nc+=cdsList[c].start-t.exons[e]->start;
	    			}
	    			mintrons++;
	    			c++;e++;
	    			continue;
	    		}
	    		else return false;
	    	}
	    	c++; //should never get here, CDS shouldn't be have to catch up with e
	    }
	    if (mintrons<cdsList.Count()-1) return false;
        //c should be cmax, e should be the last exon with CDS
	    nc+=t.exons[e]->end-cdsList[c].end;
	    for(int i=e+1;i<t.exons.Count();i++)
	    	nc+=t.exons[i]->len();
	    ncov=nc;
	    return true;
	}

	bool containedBy(GffObj& t, int& ncov) {

	    // (Warning: t may have no defined exons!)
	    //if yes: ncov will be set to the number of non-CDS-covered bases in t
	    if (t.exons.Count()<2) {
           if (t.exons.Count()==0)
               //no exons defined, just check boundaries
               return singleExonCDSMatch(t.start, t.end, ncov);
           else //single-exon
               return singleExonCDSMatch(t.exons[0]->start, t.exons[0]->end, ncov);
	    } //single or no exon
	    else { //multi-exon transcript
	    	if (start<t.exons.First()->start || end>t.exons.Last()->end)
	    		return false; //no containment possible;
	    	if (cdsList.Count()==1)
	    		return singleCDStoExon(t, ncov);
            //check intron compatibility!
	    }
	    return true;
	}
};


#endif
