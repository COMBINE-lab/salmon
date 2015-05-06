#ifndef GFF_UTILS_H
#define GFF_UTILS_H
#include "gff.h"
#include "GStr.h"
#include "GFastaIndex.h"
#include "GFaSeqGet.h"

typedef bool GFValidateFunc(GffObj* gf, GList<GffObj>* gfadd);

class GeneInfo { //for Ensembl GTF conversion
 public:
   int flag;
   GffObj* gf;
   GList<GStr> gene_names;
   GList<GStr> transcripts; //list of transcript IDs
   GeneInfo():gene_names(true, true, true), transcripts(true,true,true) {
     gf=NULL;
     flag=0;
     }
   GeneInfo(GffObj* gfrec, bool ensembl_convert=false):gene_names(true, true, true), 
                    transcripts(true,true,true) {
     flag=0;
     if (gfrec->getGeneName())
        gene_names.Add(new GStr(gfrec->getGeneName()));
     transcripts.Add(new GStr(gfrec->getID()));
     create_gf(gfrec, ensembl_convert);
     }
     
   void create_gf(GffObj* gfrec, bool ensembl_convert) {
     gf=new GffObj(gfrec->getGeneID());
     gf->gseq_id=gfrec->gseq_id;
     gf->track_id=gfrec->track_id;
     gf->start=gfrec->start;
     gf->end=gfrec->end;
     gf->strand=gfrec->strand;
     gf->setFeatureName("gene");
     gf->isGene(true);
     gf->isUsed(true);
     gf->uptr=this;
     gfrec->incLevel();
     gfrec->parent=gf;
     gf->children.Add(gfrec);
     if (ensembl_convert) {
       //gf->addAttr("type", gf->getTrackName());
       const char* biotype=gfrec->getAttr("type");
       if (biotype) gf->addAttr("type", biotype);
       }
     //gf->children.Add(gfrec);
     }
   //~GeneInfo() {
   //  }
   void update(GffObj* gfrec) {
     if (transcripts.AddedIfNew(new GStr(gfrec->getID()))<0)
       return;
     gene_names.AddedIfNew(new GStr(gfrec->getGeneName()));
     if (gf==NULL) {
        GError("GeneInfo::update() called on uninitialized gf!\n");
        //create_gf(gfrec);
        //return;
        }
     gfrec->parent=gf;
     gf->children.Add(gfrec);
     gfrec->incLevel();
     if (gf->start>gfrec->start) 
           gf->start=gfrec->start;
     if (gf->end<gfrec->end) 
           gf->end=gfrec->end;
     }
    void finalize() {
     //prepare attributes for printing
     //must be called right before printing
     if (gf==NULL || transcripts.Count()==0) return;
     if (gene_names.Count()>0) {
       gf->addAttr("Name", gene_names[0]->chars());
       /*
       GStr s(gene_names[0]->chars());
       for (int i=1;i<gene_names.Count();i++) {
          s.append(",");
          s.append(gene_names[i]->chars());
          }
       gf->addAttr("genes", s.chars());
       */
       } //has gene names
       GStr t(transcripts[0]->chars());
       for (int i=1;i<transcripts.Count();i++) {
          t.append(",");
          t.append(transcripts[i]->chars());
          }
       gf->addAttr("transcripts", t.chars());
     }
};

//genomic fasta sequence handling
class GFastaDb {
 public:
  char* fastaPath;
  GFastaIndex* faIdx; //could be a cdb .cidx file
  int last_fetchid;
  GFaSeqGet* faseq;
  //GCdbYank* gcdb;
  char* getFastaFile(int gseq_id) {
     if (fastaPath==NULL) return NULL;
     GStr s(fastaPath);
     s.trimR('/');
     s.appendfmt("/%s",GffObj::names->gseqs.getName(gseq_id));
     GStr sbase(s);
     if (!fileExists(s.chars())) s.append(".fa");
     if (!fileExists(s.chars())) s.append("sta");
     if (fileExists(s.chars())) return Gstrdup(s.chars());
         else {
             GMessage("Warning: cannot find genomic sequence file %s{.fa,.fasta}\n",sbase.chars());
             return NULL;
             }
     }

   GFastaDb(const char* fpath=NULL) {
     //gcdb=NULL;
     fastaPath=NULL;
     faseq=NULL;
     faIdx=NULL;
     init(fpath);
     }

   void init(const char* fpath) {
     if (fpath==NULL || fpath[0]==0) return;
     last_fetchid=-1;
     if (!fileExists(fpath))
       GError("Error: file/directory %s does not exist!\n",fpath);
     fastaPath=Gstrdup(fpath);
     GStr gseqpath(fpath);
     if (fileExists(fastaPath)>1) { //exists and it's not a directory
            GStr fainame(fastaPath);
            if (fainame.rindex(".fai")==fainame.length()-4) {
               //.fai index file given directly
               fastaPath[fainame.length()-4]=0;
               if (!fileExists(fastaPath))
                  GError("Error: cannot find fasta file for index %s !\n", fastaPath);
               }
              else fainame.append(".fai");
            //GMessage("creating GFastaIndex with fastaPath=%s, fainame=%s\n", fastaPath, fainame.chars());
            faIdx=new GFastaIndex(fastaPath,fainame.chars());
            GStr fainamecwd(fainame);
            int ip=-1;
            if ((ip=fainamecwd.rindex(CHPATHSEP))>=0)
               fainamecwd.cut(0,ip+1);
            if (!faIdx->hasIndex()) { //could not load index
               //try current directory
                  if (fainame!=fainamecwd) {
                    if (fileExists(fainamecwd.chars())>1) {
                       faIdx->loadIndex(fainamecwd.chars());
                       }
                    }
                  } //tried to load index
            if (!faIdx->hasIndex()) {
                 GMessage("No fasta index found for %s. Rebuilding, please wait..\n",fastaPath);
                 faIdx->buildIndex();
                 if (faIdx->getCount()==0) GError("Error: no fasta records found!\n");
                 GMessage("Fasta index rebuilt.\n");
                 FILE* fcreate=fopen(fainame.chars(), "w");
                 if (fcreate==NULL) {
                   GMessage("Warning: cannot create fasta index %s! (permissions?)\n", fainame.chars());
                   if (fainame!=fainamecwd) fcreate=fopen(fainamecwd.chars(), "w");
                   if (fcreate==NULL)
                      GError("Error: cannot create fasta index %s!\n", fainamecwd.chars());
                   }
                 if (faIdx->storeIndex(fcreate)<faIdx->getCount())
                     GMessage("Warning: error writing the index file!\n");
                 } //index created and attempted to store it
            } //multi-fasta
     }
   GFaSeqGet* fetch(int gseq_id, bool checkFasta=false) {
     if (fastaPath==NULL) return NULL;
     if (gseq_id==last_fetchid && faseq!=NULL) return faseq;
     delete faseq;
     faseq=NULL;
     last_fetchid=-1;
     char* gseqname=GffObj::names->gseqs.getName(gseq_id);
     if (faIdx!=NULL) { //fastaPath was the multi-fasta file name
        GFastaRec* farec=faIdx->getRecord(gseqname);
        if (farec!=NULL) {
             faseq=new GFaSeqGet(fastaPath,farec->seqlen, farec->fpos,
                               farec->line_len, farec->line_blen);
             faseq->loadall(); //just cache the whole sequence, it's faster
             last_fetchid=gseq_id;
             }
        else {
          GMessage("Warning: couldn't find fasta record for '%s'!\n",gseqname);
          return NULL;
          }
        }
     else {
         char* sfile=getFastaFile(gseq_id);
         if (sfile!=NULL) {
            faseq=new GFaSeqGet(sfile,checkFasta);
            faseq->loadall();
            last_fetchid=gseq_id;
            GFREE(sfile);
            }
         } //one fasta file per contig
       return faseq;
     }

   ~GFastaDb() {
     GFREE(fastaPath);
     //delete gcdb;
     delete faIdx;
     delete faseq;
     }
};

class GffLocus;

class GTData { //transcript associated data
 public:
    GffObj* rna;
    GffLocus* locus;
    GffObj* replaced_by;
    GeneInfo* geneinfo;
    int flag;
    GTData(GffObj* t=NULL) {
        rna=t;
        flag=0;
        locus=NULL;
        replaced_by=NULL;
        geneinfo=NULL;
        if (rna!=NULL) {
            geneinfo=(GeneInfo*)rna->uptr; //take over geneinfo, if there
            rna->uptr=this;
            }
        }
   bool operator<(GTData& b) { return (rna < b.rna); }
   bool operator==(GTData& b) { return (rna==b.rna); }
};

class CGeneSym {
 public:
  GStr name;
  int freq;
  CGeneSym(const char* n=NULL, int f=0):name(n) {
    freq=f;
    }
  bool operator<(CGeneSym& b) {
     return (freq==b.freq)? ( (name.length()==b.name.length()) ? (name<b.name) :
         (name.length()<b.name.length()) ) : ( freq>b.freq );
     }
  bool operator==(CGeneSym& b) { return name==b.name; }
};

const char* getGeneDescr(const char* gsym);

void printLocus(GffLocus* loc, const char* pre=NULL);

class GffLocus:public GSeg {
public:
    int gseq_id; //id of underlying genomic sequence
    int locus_num;
    bool is_mrna;
    char strand;
    GffObj* t_maxcov;  //transcript with maximum coverage (for main "ref" transcript)
    GList<GffObj> rnas; //list of transcripts (isoforms) for this locus
    GArray<GSeg> mexons; //list of merged exons in this region
    GList<CGeneSym> gene_names;
    GList<CGeneSym> gene_ids;
    int v; //user flag/data
   /*
   bool operator==(GffLocus& d){
       return (gseq_id==d.gseq_id && strand==d.strand && start==d.start && end==d.end);
       }
   bool operator<(GffLocus& d){
     if (gseq_id!=d.gseq_id) return (gseq_id<d.gseq_id);
     if (start==d.start) {
        if (end==d.end) return strand<d.strand;
                     else return end<d.end;
        } else return (start<d.start);
     }
    */
    const char* getGeneName() {
         if (gene_names.Count()==0) return NULL;
         return gene_names.First()->name.chars();
         }
    const char* get_tmax_id() {
         return t_maxcov->getID();
         }
    const char* get_descr() {
       if (gene_names.Count()>0) {
          for (int i=0;i<gene_names.Count();i++) {
            const char* gn=getGeneDescr(gene_names.First()->name.chars());
            if (gn!=NULL) return gn;
            }
          }
       char* s=t_maxcov->getAttr("product");
       if (s!=NULL) return s;
       s=t_maxcov->getAttr("descr");
       if (s!=NULL) return s;
       s=t_maxcov->getAttr("description");
       if (s!=NULL) return s;
       s=t_maxcov->getAttr("info");
       if (s!=NULL) return s;
       return NULL;
       }

    GffLocus(GffObj* t=NULL):rnas(true,false,false),mexons(true,true),
           gene_names(true,true,false), gene_ids(true,true,false) {
        //this will NOT free rnas!
        t_maxcov=NULL;
        gseq_id=-1;
        v=0;
        locus_num=0;
        start=0;
        end=0;
        strand=0;
        is_mrna=false;
        if (t!=NULL) {
           start=t->exons.First()->start;
           end=t->exons.Last()->end;;
           gseq_id=t->gseq_id;
           GSeg seg;
           for (int i=0;i<t->exons.Count();i++) {
                seg.start=t->exons[i]->start;
                seg.end=t->exons[i]->end;
                mexons.Add(seg);
                }
           rnas.Add(t);
           ((GTData*)(t->uptr))->locus=this;
           t_maxcov=t;
           strand=t->strand;
           if (t->ftype_id==gff_fid_mRNA) {
              is_mrna=true;
              }
           }
    }

   void addMerge(GffLocus& locus, GffObj* lnkrna) {
     //add all the elements of the other locus (merging)
     //-- merge mexons
     GArray<int> ovlexons(true,true); //list of locus.mexons indexes overlapping existing mexons
     int i=0; //index of first mexons with a merge
     int j=0; //index current mrna exon
     while (i<mexons.Count() && j<locus.mexons.Count()) {
            uint istart=mexons[i].start;
            uint iend=mexons[i].end;
            uint jstart=locus.mexons[j].start;
            uint jend=locus.mexons[j].end;
            if (iend<jstart) { i++; continue; }
            if (jend<istart) { j++; continue; }
            ovlexons.Add(j);
            //extend mexons[i] as needed
            if (jstart<istart) mexons[i].start=jstart;
            if (jend>iend) { //mexons[i] end extend
                mexons[i].end=jend;
                //now this could overlap the next mexon(s), so we have to merge them all
                while (i<mexons.Count()-1 && mexons[i].end>mexons[i+1].start) {
                    uint nextend=mexons[i+1].end;
                    mexons.Delete(i+1);
                    if (nextend>mexons[i].end) {
                        mexons[i].end=nextend;
                        break; //no need to check next mexons
                    }
                } //while next mexons merge
            } // mexons[i] end extend
            j++; //check the next locus.mexon
        }
        //-- add the rest of the non-overlapping mexons:
        GSeg seg;
        for (int i=0;i<locus.mexons.Count();i++) {
            seg.start=locus.mexons[i].start;
            seg.end=locus.mexons[i].end;
            if (!ovlexons.Exists(i)) mexons.Add(seg);
        }
     // -- add locus.rnas
     for (int i=0;i<locus.rnas.Count();i++) {
          ((GTData*)(locus.rnas[i]->uptr))->locus=this;
          if (locus.rnas[i]!=lnkrna) rnas.Add(locus.rnas[i]);
          }
        // -- adjust start/end as needed
     if (start>locus.start) start=locus.start;
     if (end<locus.end) end=locus.end;
     if (locus.is_mrna) is_mrna=true;
     if (t_maxcov->covlen<locus.t_maxcov->covlen)
            t_maxcov=locus.t_maxcov;
     }

    bool exonOverlap(GffLocus& loc) {
        //check if any mexons overlap!
        if (strand!=loc.strand || loc.start>end || start>loc.end) return false;
        int i=0;
        int j=0;
        while (i<mexons.Count() && j<loc.mexons.Count()) {
            uint istart=mexons[i].start;
            uint iend=mexons[i].end;
            uint jstart=loc.mexons[j].start;
            uint jend=loc.mexons[j].end;
            if (iend<jstart) { i++; continue; }
            if (jend<istart) { j++; continue; }
            //exon overlap found if we're here:
            return true;
        }
        return false;
    }

    bool add_RNA(GffObj* t) {
        //if (rnas.Count()==0) return true; //? should never be called on an empty locus
        if (t->gseq_id!=gseq_id || t->strand!=strand || t->start>end || start>t->end)
              return false; //rna must be on the same genomic seq
        //check for exon overlap with existing mexons
        //also update mexons accordingly if t is to be added
        bool hasovl=false;
        int i=0; //index of first mexons with a merge
        int j=0; //index current t exon
        GArray<int> ovlexons(true,true); //list of mrna exon indexes overlapping mexons
        while (i<mexons.Count() && j<t->exons.Count()) {
            uint istart=mexons[i].start;
            uint iend=mexons[i].end;
            uint jstart=t->exons[j]->start;
            uint jend=t->exons[j]->end;
            if (iend<jstart) { i++; continue; }
            if (jend<istart) { j++; continue; }
            //exon overlap found if we're here:
            ovlexons.Add(j);
            hasovl=true;
            //extend mexons[i] as needed
            if (jstart<istart) mexons[i].start=jstart;
            if (jend>iend) { //mexon stretch up
                mexons[i].end=jend;
                //now this could overlap the next mexon(s), so we have to merge them all
                while (i<mexons.Count()-1 && mexons[i].end>mexons[i+1].start) {
                    uint nextend=mexons[i+1].end;
                    mexons.Delete(i+1);
                    if (nextend>mexons[i].end) {
                        mexons[i].end=nextend;
                        break; //no need to check next mexons
                    }
                } //while next mexons merge
            } //possible mexons merge

            j++; //check the next t exon
        }//all vs all exon check loop
        if (hasovl) {
            GSeg seg;
             //add the rest of the non-overlapping exons
            for (int i=0;i<t->exons.Count();i++) {
                seg.start=t->exons[i]->start;
                seg.end=t->exons[i]->end;
                if (!ovlexons.Exists(i)) mexons.Add(seg);
                }
            rnas_add(t);
            // add to rnas
            ((GTData*)t->uptr)->locus=this;
            gseq_id=t->gseq_id;
            }
        return hasovl;
    }

    //simpler,basic adding of a mrna
    void rnas_add(GffObj* t) {
      rnas.Add(t);
      // adjust start/end
      //if (start==0 || start>t->start) start=t->start;
      if (start==0) start=t->start;
        else if (start>t->start) {
          start=t->start;
          }
      if (end<t->end) end=t->end;
      if (t_maxcov->covlen<t->covlen) t_maxcov=t;
      if (strand==0) strand=t->strand;
      if (t->ftype_id==gff_fid_mRNA) is_mrna=true;
      }
};

class GenomicSeqData {
  int gseq_id;
 public:
  const char* gseq_name;
  GList<GffObj> gfs; //all non-transcript features -> usually gene features
  GList<GffObj> rnas; //all transcripts on this genomic sequence
  GList<GffLocus> loci; //all loci clusters
  GList<GTData> tdata; //transcript data (uptr holder for all rnas loaded here)
  //GenomicSeqData(int gid=-1):rnas(true,true,false),loci(true,true,true),
  GenomicSeqData(int gid=-1):gfs(true, true, false),rnas((GCompareProc*)gfo_cmpByLoc),loci(true,true,false),
       tdata(false,true,false) {
  gseq_id=gid;
  if (gseq_id>=0) 
    gseq_name=GffObj::names->gseqs.getName(gseq_id);
  
  }
  bool operator==(GenomicSeqData& d){
    return gseq_id==d.gseq_id;
  }
  bool operator<(GenomicSeqData& d){
    return (gseq_id<d.gseq_id);
  }
};

int gseqCmpName(const pointer p1, const pointer p2);

class GSpliceSite {
 public:
  char nt[3];
  GSpliceSite(const char* c, bool revc=false) {
    nt[2]=0;
    if (c==NULL) {
      nt[0]=0;
      nt[1]=0;
      return;
      }
    if (revc) {
      nt[0]=toupper(ntComplement(c[1]));
      nt[1]=toupper(ntComplement(c[0]));
      }
    else {
      nt[0]=toupper(c[0]);
      nt[1]=toupper(c[1]);
      }
    }

  GSpliceSite(const char* intron, int intronlen, bool getAcceptor, bool revc=false) {
    nt[2]=0;
    if (intron==NULL || intronlen==0)
       GError("Error: invalid intron or intron len for GSpliceSite()!\n");
    const char* c=intron;
    if (revc) {
      if (!getAcceptor) c+=intronlen-2;
      nt[0]=toupper(ntComplement(c[1]));
      nt[1]=toupper(ntComplement(c[0]));
      }
    else { //on forward strand
      if (getAcceptor) c+=intronlen-2;
      nt[0]=toupper(c[0]);
      nt[1]=toupper(c[1]);
      }//forward strand
    }

  GSpliceSite(const char n1, const char n2) {
    nt[2]=0;
    nt[0]=toupper(n1);
    nt[1]=toupper(n2);
    }
  bool canonicalDonor() {
    return (nt[0]=='G' && (nt[1]=='C' || nt[1]=='T'));
    }
  bool operator==(GSpliceSite& c) {
    return (c.nt[0]==nt[0] && c.nt[1]==nt[1]);
    }
  bool operator==(GSpliceSite* c) {
    return (c->nt[0]==nt[0] && c->nt[1]==nt[1]);
    }
  bool operator==(const char* c) {
    //return (nt[0]==toupper(c[0]) && nt[1]==toupper(c[1]));
    //assumes given const nucleotides are uppercase already!
    return (nt[0]==c[0] && nt[1]==c[1]);
    }
  bool operator!=(const char* c) {
    //assumes given const nucleotides are uppercase already!
    return (nt[0]!=c[0] || nt[1]!=c[1]);
    }
};

struct GffLoader {
  GStr fname;
  FILE* f;
  bool transcriptsOnly;
  bool fullAttributes;
  bool noExonAttrs;
  bool mergeCloseExons;
  bool showWarnings;
  bool noPseudo;
  void placeGf(GffObj* t, GenomicSeqData* gdata, bool doCluster=true, bool collapseRedundant=true,
                                    bool matchAllIntrons=true, bool fuzzSpan=false);
  void load(GList<GenomicSeqData>&seqdata, GFValidateFunc* gf_validate=NULL, 
                      bool doCluster=true, bool doCollapseRedundant=true, 
                      bool matchAllIntrons=true, bool fuzzSpan=false, bool forceExons=false);
  GffLoader(const char* filename):fname(filename) {
      f=NULL;
      transcriptsOnly=true;
      fullAttributes=false;
      noExonAttrs=false;
      mergeCloseExons=false;
      showWarnings=false;
      noPseudo=false;
      if (fname=="-" || fname=="stdin") {
         f=stdin;
         fname="stdin";
         }
        else {
          if ((f=fopen(fname.chars(), "r"))==NULL) {
            GError("Error: cannot open gff file %s!\n",fname.chars());
            }
          }
      }
  ~GffLoader() {
      if (f!=NULL && f!=stdin) fclose(f);
      }
};

void printFasta(FILE* f, GStr& defline, char* seq, int seqlen=-1);

//"position" a given coordinate x within a list of transcripts sorted by their start (lowest)
//coordinate, using quick-search; the returned int is the list index of the closest *higher*
//GffObj - i.e. starting right *ABOVE* the given coordinate
//Convention: returns -1 if there is no such GffObj (i.e. last GffObj starts below x)
int qsearch_rnas(uint x, GList<GffObj>& rnas);
int qsearch_gloci(uint x, GList<GffLocus>& loci);

GffObj* redundantTranscripts(GffObj& ti, GffObj&  tj, bool matchAllIntrons=true, bool fuzzSpan=false);

//void loadGFF(FILE* f, GList<GenomicSeqData>& seqdata, const char* fname);

void collectLocusData(GList<GenomicSeqData>& ref_data);

#endif
