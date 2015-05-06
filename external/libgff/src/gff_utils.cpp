#include "gff_utils.h"

extern bool verbose;
extern bool debugMode;

//bool debugState=false;

void printFasta(FILE* f, GStr& defline, char* seq, int seqlen) {
 if (seq==NULL) return;
 int len=(seqlen>0)?seqlen:strlen(seq);
 if (len<=0) return;
 if (!defline.is_empty())
     fprintf(f, ">%s\n",defline.chars());
 int ilen=0;
 for (int i=0; i < len; i++, ilen++) {
   if (ilen == 70) {
     fputc('\n', f);
     ilen = 0;
     }
   putc(seq[i], f);
   } //for
 fputc('\n', f);
}

int qsearch_gloci(uint x, GList<GffLocus>& loci) {
  //binary search
  //do the simplest tests first:
  if (loci[0]->start>x) return 0;
  if (loci.Last()->start<x) return -1;
  uint istart=0;
  int i=0;
  int idx=-1;
  int maxh=loci.Count()-1;
  int l=0;
  int h = maxh;
  while (l <= h) {
     i = (l+h)>>1;
     istart=loci[i]->start;
     if (istart < x)  l = i + 1;
          else {
             if (istart == x) { //found matching coordinate here
                  idx=i;
                  while (idx<=maxh && loci[idx]->start==x) {
                     idx++;
                     }
                  return (idx>maxh) ? -1 : idx;
                  }
             h = i - 1;
             }
     } //while
 idx = l;
 while (idx<=maxh && loci[idx]->start<=x) {
    idx++;
    }
 return (idx>maxh) ? -1 : idx;
}

int qsearch_rnas(uint x, GList<GffObj>& rnas) {
  //binary search
  //do the simplest tests first:
  if (rnas[0]->start>x) return 0;
  if (rnas.Last()->start<x) return -1;
  uint istart=0;
  int i=0;
  int idx=-1;
  int maxh=rnas.Count()-1;
  int l=0;
  int h = maxh;
  while (l <= h) {
     i = (l+h)>>1;
     istart=rnas[i]->start;
     if (istart < x)  l = i + 1;
          else {
             if (istart == x) { //found matching coordinate here
                  idx=i;
                  while (idx<=maxh && rnas[idx]->start==x) {
                     idx++;
                     }
                  return (idx>maxh) ? -1 : idx;
                  }
             h = i - 1;
             }
     } //while
 idx = l;
 while (idx<=maxh && rnas[idx]->start<=x) {
    idx++;
    }
 return (idx>maxh) ? -1 : idx;
}

int cmpRedundant(GffObj& a, GffObj& b) {
  if (a.exons.Count()==b.exons.Count()) {
     if (a.covlen==b.covlen) {
       return strcmp(a.getID(), b.getID());
       }
     else return (a.covlen>b.covlen)? 1 : -1;
     }
   else return (a.exons.Count()>b.exons.Count())? 1: -1;
}


bool tMatch(GffObj& a, GffObj& b) {
  //strict intron chain match, or single-exon perfect match
  int imax=a.exons.Count()-1;
  int jmax=b.exons.Count()-1;
  int ovlen=0;
  if (imax!=jmax) return false; //different number of introns

  if (imax==0) { //single-exon mRNAs
    //if (equnspl) {
      //fuzz match for single-exon transfrags: 
      // it's a match if they overlap at least 80% of max len
      ovlen=a.exons[0]->overlapLen(b.exons[0]);
      int maxlen=GMAX(a.covlen,b.covlen);
      return (ovlen>=maxlen*0.8);
    /*}
    else {
      //only exact match
      ovlen=a.covlen;
      return (a.exons[0]->start==b.exons[0]->start &&
          a.exons[0]->end==b.exons[0]->end);
      
       }*/
     }
  //check intron overlaps
  ovlen=a.exons[0]->end-(GMAX(a.start,b.start))+1;
  ovlen+=(GMIN(a.end,b.end))-a.exons.Last()->start;
  for (int i=1;i<=imax;i++) {
    if (i<imax) ovlen+=a.exons[i]->len();
    if ((a.exons[i-1]->end!=b.exons[i-1]->end) ||
      (a.exons[i]->start!=b.exons[i]->start)) {
            return false; //intron mismatch
    }
  }
  return true;
}


bool unsplContained(GffObj& ti, GffObj&  tj, bool fuzzSpan) {
 //returns true only if ti (which MUST be single-exon) is "almost" contained in any of tj's exons
 //but it does not cross any intron-exon boundary of tj
  int imax=ti.exons.Count()-1;
  int jmax=tj.exons.Count()-1;
  if (imax>0) GError("Error: bad unsplContained() call, 1st param must be single-exon transcript!\n");
  int minovl = (int)(0.8 * ti.len()); //minimum overlap for fuzzSpan
  if (fuzzSpan) {
    for (int j=0;j<=jmax;j++) {
       //must NOT overlap the introns
       if ((j>0 && ti.start<tj.exons[j]->start) 
          || (j<jmax && ti.end>tj.exons[j]->end))
         return false;
       if (ti.exons[0]->overlapLen(tj.exons[j])>=minovl)
              return true;
       }
      } else {
    for (int j=0;j<=jmax;j++) {
       //must NOT overlap the introns
       if ((j>0 && ti.start<tj.exons[j]->start) 
          || (j<jmax && ti.end>tj.exons[j]->end))
         return false;
         //strict containment
       if (ti.end<=tj.exons[j]->end && ti.start>=tj.exons[j]->start) 
            return true;
       }
      }
 return false;
}

GffObj* redundantTranscripts(GffObj& ti, GffObj&  tj, bool matchAllIntrons, bool fuzzSpan) {
  // matchAllIntrons==true:  transcripts are considered "redundant" only if
  //                   they have the exact same number of introns and same splice sites (or none)
  //                 (single-exon transcripts can be also fully contained to be considered matching)
  // matchAllIntrons==false: an intron chain could be a subset of a "container" chain, 
  //                   as long as no intron-exon boundaries are violated; also, a single-exon 
  //                   transcript will be collapsed if it's contained in one of the exons of the other
  // fuzzSpan==false: the genomic span of one transcript must be contained in or equal with the genomic 
  //                  span of the other 
  // 
  // fuzzSpan==true: then genomic spans of transcripts are no longer required to be fully contained 
  //                 (i.e. they may extend each-other in opposite directions)
  
  //if redundancy is detected, the "bigger" transcript is returned (otherwise NULL is returned)
 if (ti.start>=tj.end || tj.start>=ti.end || tj.strand!=ti.strand) return NULL; //no span overlap at all
 int imax=ti.exons.Count()-1;
 int jmax=tj.exons.Count()-1;
 GffObj* bigger=NULL;
 GffObj* smaller=NULL;
 if (matchAllIntrons) {
   if (imax!=jmax) return NULL;
   if (ti.covlen>tj.covlen) {
       bigger=&ti;
       if (!fuzzSpan && (ti.start>tj.start || ti.end<tj.end)) return NULL;
       }
     else { //ti.covlen<=tj.covlen
       bigger=&tj;
       if (!fuzzSpan && (tj.start>ti.start || tj.end<ti.end)) return NULL;
       }
   //check that all introns really match
   for (int i=0;i<imax;i++) {
     if (ti.exons[i]->end!=tj.exons[i]->end || 
         ti.exons[i+1]->start!=tj.exons[i+1]->start) return NULL;
     }
   return bigger;
   }
 //--- matchAllIntrons==false: intron-chain containment is also considered redundancy
 //int maxlen=0;
 int minlen=0;
 if (ti.covlen>tj.covlen) {
      if (tj.exons.Count()>ti.exons.Count()) {
          //exon count override
          bigger=&tj;
          smaller=&ti;
          }
        else {
          bigger=&ti;
          smaller=&tj;
          }
      //maxlen=ti.covlen;
      minlen=tj.covlen;
      }
   else { //tj has more bases
      if (ti.exons.Count()>tj.exons.Count()) {
          //exon count override
          bigger=&ti;
          smaller=&tj;
          }
        else {
          bigger=&tj;
          smaller=&ti;
          }
      //maxlen=tj.covlen;
      minlen=ti.covlen;
      }
 if (imax==0 && jmax==0) {
     //single-exon transcripts: if fuzzSpan, at least 80% of the shortest one must be overlapped by the other
     if (fuzzSpan) {
         return (ti.exons[0]->overlapLen(tj.exons[0])>=minlen*0.8) ? bigger : NULL;
         }
       else {
         return (smaller->start>=bigger->start && smaller->end<=bigger->end) ? bigger : NULL;
         }
     }
 //containment is also considered redundancy
 if (smaller->exons.Count()==1) {
   //check if this single exon is contained in any of tj exons
   //without violating any intron-exon boundaries
   return (unsplContained(*smaller, *bigger, fuzzSpan) ? bigger : NULL);
   }

 //--from here on: both are multi-exon transcripts, imax>0 && jmax>0
  if (ti.exons[imax]->start<tj.exons[0]->end ||
     tj.exons[jmax]->start<ti.exons[0]->end )
         return NULL; //intron chains do not overlap at all
 
 
 //checking full intron chain containment
 uint eistart=0, eiend=0, ejstart=0, ejend=0; //exon boundaries
 int i=1; //exon idx to the right of the current intron of ti
 int j=1; //exon idx to the right of the current intron of tj
 //find the first intron overlap:
 while (i<=imax && j<=jmax) {
    eistart=ti.exons[i-1]->end;
    eiend=ti.exons[i]->start;
    ejstart=tj.exons[j-1]->end;
    ejend=tj.exons[j]->start;
    if (ejend<eistart) { j++; continue; }
    if (eiend<ejstart) { i++; continue; }
    //we found an intron overlap
    break;
    }
 if (!fuzzSpan && (bigger->start>smaller->start || bigger->end < smaller->end)) return NULL;
 if ((i>1 && j>1) || i>imax || j>jmax) {
     return NULL; //either no intron overlaps found at all
                  //or it's not the first intron for at least one of the transcripts
     }
 if (eistart!=ejstart || eiend!=ejend) return NULL; //not an exact intron match
 if (j>i) {
   //i==1, ti's start must not conflict with the previous intron of tj
   if (ti.start<tj.exons[j-1]->start) return NULL;
   //so i's first intron starts AFTER j's first intron
   // then j must contain i, so i's last intron must end with or before j's last intron
   if (ti.exons[imax]->start>tj.exons[jmax]->start) return NULL;
      //comment out the line above if you just want "intron compatibility" (i.e. extension of intron chains )
   }
  else if (i>j) {
     //j==1, tj's start must not conflict with the previous intron of ti
     if (tj.start<ti.exons[i-1]->start) return NULL;
     //so j's intron chain starts AFTER i's
     // then i must contain j, so j's last intron must end with or before j's last intron
     if (tj.exons[jmax]->start>ti.exons[imax]->start) return NULL;
        //comment out the line above for just "intronCompatible()" check (allowing extension of intron chain)
     }
 //now check if the rest of the introns overlap, in the same sequence
 i++;
 j++;
 while (i<=imax && j<=jmax) {
  if (ti.exons[i-1]->end!=tj.exons[j-1]->end ||
      ti.exons[i]->start!=tj.exons[j]->start) return NULL;
  i++;
  j++;
  }
 i--;
 j--;
 if (i==imax && j<jmax) {
   // tj has more introns to the right, check if ti's end doesn't conflict with the current tj exon boundary
   if (ti.end>tj.exons[j]->end) return NULL;
   }
 else if (j==jmax && i<imax) {
   if (tj.end>ti.exons[i]->end) return NULL;
   }
 return bigger;
}


int gseqCmpName(const pointer p1, const pointer p2) {
 return strcmp(((GenomicSeqData*)p1)->gseq_name, ((GenomicSeqData*)p2)->gseq_name);
}


void printLocus(GffLocus* loc, const char* pre) {
  if (pre!=NULL) fprintf(stderr, "%s", pre);
  GMessage(" [%d-%d] : ", loc->start, loc->end);
  GMessage("%s",loc->rnas[0]->getID());
  for (int i=1;i<loc->rnas.Count();i++) {
    GMessage(",%s",loc->rnas[i]->getID());
    }
  GMessage("\n");
}

void preserveContainedCDS(GffObj* t, GffObj* tfrom) {
 //transfer CDS info to the container t if it's a larger protein
 if (tfrom->CDstart==0) return;
 if (t->CDstart) {
   if (tfrom->CDstart<t->CDstart && tfrom->CDstart>=t->start)
      t->CDstart=tfrom->CDstart;
   if (tfrom->CDend>t->CDend && tfrom->CDend<=t->end)
      t->CDend=tfrom->CDend;
   }
  else { //no CDS info on container, just copy it from the contained
   t->addCDS(tfrom->CDstart, tfrom->CDend, tfrom->CDphase);
   }
}

bool exonOverlap2Gene(GffObj* t, GffObj& g) {
	if (t->exons.Count()>0) {
		return t->exonOverlap(g.start, g.end);
	}
	else return g.overlap(*t);
}
void GffLoader::placeGf(GffObj* t, GenomicSeqData* gdata, bool doCluster, bool collapseRedundant,
                                               bool matchAllIntrons, bool fuzzSpan) {
  GTData* tdata=new GTData(t); //additional transcript data
  gdata->tdata.Add(tdata);
  //int tidx=-1;
  /*
  if (debug) {
     GMessage(">>Placing transcript %s\n", t->getID());
     debugState=true;
     }
    else debugState=false; 
   */
  //dumb TRNA case for RefSeq: gene parent link missing
  //try to restore it here; BUT this only works if gene feature comes first
  if (t->parent==NULL && t->isTranscript()) {
  	int gidx=gdata->gfs.Count()-1;
  	while (gidx>=0 && gdata->gfs[gidx]->end>=t->start) {
  		GffObj& g = *(gdata->gfs[gidx]);
  		if (g.isGene() && t->strand==g.strand && exonOverlap2Gene(t, g)) {
  			g.children.Add(t);
  			t->parent=&g;
  			//disable printing of gene if transcriptsOnly
  			if (transcriptsOnly) {
  				g.udata|=4; //tag it as non-printable
  			}
  			const char* geneName=g.getAttr("Name");
  			if (t->getAttr("Name")==NULL && geneName) {
  				t->addAttr("Name", geneName);
  				t->addAttr("gene_name", geneName);
  			}
  			t->addAttr("geneID", g.getID());
  			break;
  		}
  		gidx--;
  	}
  }

  /*
	if (t->exons.Count()==0  && t->children.Count()==0 && forceExons) {
		//a non-mRNA feature with no subfeatures
		//just so we get some sequence functions working, add a dummy "exon"-like subfeature here
		//--this could be a single "pseudogene" entry or another genomic region without exons
		//
		t->addExon(t->start,t->end);
	}
  */
  if (t->exons.Count()>0) {
              //tidx=
              gdata->rnas.Add(t); //added it in sorted order
              }
            else {
              if (t->isGene() || !this->transcriptsOnly)
              	  gdata->gfs.Add(t);
              return; //nothing to do with these non-transcript objects
              }
  if (!doCluster) return;
  if (gdata->loci.Count()==0) {
       gdata->loci.Add(new GffLocus(t));
       //GMessage("  <<make it first locus %d-%d \n",t->start, t->end);
       return;
       }
   /*    
  //DEBUG: show available loci:
   if (debug) {
    GMessage("  [%d loci already:\n", gdata->loci.Count());
    for (int l=0;l<gdata->loci.Count();l++) {
       printLocus(gdata->loci[l]);
       }
    }
  */
  int nidx=qsearch_gloci(t->end, gdata->loci); //get index of nearest locus starting just ABOVE t->end
  //GMessage("\tlooking up end coord %d in gdata->loci.. (qsearch got nidx=%d)\n", t->end, nidx);
  if (nidx==0) {
     //cannot have any overlapping loci
     //if (debug) GMessage("  <<no ovls possible, create locus %d-%d \n",t->start, t->end);
     gdata->loci.Add(new GffLocus(t));
     return;
     }
  if (nidx==-1) nidx=gdata->loci.Count();//all loci start below t->end
  int lfound=0; //count of parent loci
  GArray<int> mrgloci(false);
  GList<GffLocus> tloci(true); //candidate parent loci to adopt this
  //if (debug) GMessage("\tchecking all loci from %d to 0\n",nidx-1);
  for (int l=nidx-1;l>=0;l--) {
      GffLocus& loc=*(gdata->loci[l]);
      if (loc.strand!='.' && t->strand!='.'&& loc.strand!=t->strand) continue;
      if (t->start>loc.end) {
           if (t->start-loc.start>GFF_MAX_LOCUS) break; //give up already
           continue;
           }
      if (loc.start>t->end) {
               //this should never be the case if nidx was found correctly
               GMessage("Warning: qsearch_gloci found loc.start>t.end!(t=%s)\n", t->getID());
               continue;
               }
      /*
      if (debug) {
          GMessage(" !range overlap found with locus ");
          printLocus(&loc);
          }
      */
      if (loc.add_RNA(t)) {
         //will add this transcript to loc
         lfound++;
         mrgloci.Add(l);
         if (collapseRedundant) {
           //compare to every single transcript in this locus
           for (int ti=0;ti<loc.rnas.Count();ti++) {
                 if (loc.rnas[ti]==t) continue;
                 GTData* odata=(GTData*)(loc.rnas[ti]->uptr);
                 //GMessage("  ..redundant check vs overlapping transcript %s\n",loc.rnas[ti]->getID());
                 GffObj* container=NULL;
                 if (odata->replaced_by==NULL && 
                      (container=redundantTranscripts(*t, *(loc.rnas[ti]), matchAllIntrons, fuzzSpan))!=NULL) {
                     if (container==t) {
                        odata->replaced_by=t;
                        preserveContainedCDS(t, loc.rnas[ti]);
                        }
                     else {
                        tdata->replaced_by=loc.rnas[ti];
                        preserveContainedCDS(loc.rnas[ti], t);
                        }
                     }
              }//for each transcript in the exon-overlapping locus
          } //if doCollapseRedundant
         } //overlapping locus
      } //for each existing locus
  if (lfound==0) {
      //overlapping loci not found, create a locus with only this mRNA
      /* if (debug) {
        GMessage("  overlapping locus not found, create locus %d-%d \n",t->start, t->end);
        }
      */
      int addidx=gdata->loci.Add(new GffLocus(t));
      if (addidx<0) {
         //should never be the case!
         GMessage("  WARNING: new GffLocus(%s:%d-%d) not added!\n",t->getID(), t->start, t->end);
         }
      }
   else { //found at least one overlapping locus
     lfound--;
     int locidx=mrgloci[lfound];
     GffLocus& loc=*(gdata->loci[locidx]);
     //last locus index found is also the smallest index
     if (lfound>0) {
       //more than one loci found parenting this mRNA, merge loci
       /* if (debug)
          GMessage(" merging %d loci \n",lfound);
       */
       for (int l=0;l<lfound;l++) {
          int mlidx=mrgloci[l]; 
          loc.addMerge(*(gdata->loci[mlidx]), t);
          gdata->loci.Delete(mlidx); //highest indices first, so it's safe to remove
          }
       }
     int i=locidx;  
     while (i>0 && loc<*(gdata->loci[i-1])) {
       //bubble down until it's in the proper order
       i--;
       gdata->loci.Swap(i,i+1);
       }
     }//found at least one overlapping locus
}

void collectLocusData(GList<GenomicSeqData>& ref_data) {
  int locus_num=0;
  for (int g=0;g<ref_data.Count();g++) {
    GenomicSeqData* gdata=ref_data[g];
    for (int l=0;l<gdata->loci.Count();l++) {
      GffLocus& loc=*(gdata->loci[l]);
      GHash<int> gnames(true); //gene names in this locus
      GHash<int> geneids(true); //Entrez GeneID: numbers
      for (int i=0;i<loc.rnas.Count();i++) {
        GffObj& t=*(loc.rnas[i]);
        GStr gname(t.getGeneName());
        if (!gname.is_empty()) {
           gname.upper();
           int* prevg=gnames.Find(gname.chars());
           if (prevg!=NULL) (*prevg)++;
                  else gnames.Add(gname, new int(1));
           }
        //parse GeneID xrefs, if any:
        GStr xrefs(t.getAttr("xrefs"));
        if (!xrefs.is_empty()) {
          xrefs.startTokenize(",");
          GStr token;
          while (xrefs.nextToken(token)) {
            token.upper();
            if (token.startsWith("GENEID:")) {
              token.cut(0,token.index(':')+1);
              int* prevg=geneids.Find(token.chars());
              if (prevg!=NULL) (*prevg)++;
                     else geneids.Add(token, new int(1));
              }
            } //for each xref
          } //xrefs parsing
        }//for each transcript
      locus_num++;
      loc.locus_num=locus_num;
      if (gnames.Count()>0) { //collect all gene names associated to this locus
         gnames.startIterate();
         int* gfreq=NULL;
         char* key=NULL;
         while ((gfreq=gnames.NextData(key))!=NULL) {
            loc.gene_names.AddIfNew(new CGeneSym(key,*gfreq));
            }
         } //added collected gene_names
      if (loc.gene_ids.Count()>0) { //collect all GeneIDs names associated to this locus
         geneids.startIterate();
         int* gfreq=NULL;
         char* key=NULL;
         while ((gfreq=geneids.NextData(key))!=NULL) {
           loc.gene_ids.AddIfNew(new CGeneSym(key,*gfreq));
            }
          }
      } //for each locus
  }//for each genomic sequence
}


void GffLoader::load(GList<GenomicSeqData>& seqdata, GFValidateFunc* gf_validate, 
                          bool doCluster, bool doCollapseRedundant, 
                          bool matchAllIntrons, bool fuzzSpan, bool forceExons) {
   GffReader* gffr=new GffReader(f, this->transcriptsOnly, false); //not only mRNA features, not sorted
   gffr->showWarnings(this->showWarnings);
   //           keepAttrs   mergeCloseExons  noExonAttr
   gffr->readAll(this->fullAttributes,    this->mergeCloseExons,  this->noExonAttrs);
   GVec<int> pseudoAttrIds;
   GVec<int> pseudoFeatureIds;
   if (this->noPseudo) {
   	 GffNameList& fnames = gffr->names->feats;
   	 for (int i=0;i<fnames.Count();i++) {
   		char* n=fnames[i]->name;
   		if (startsWith(n, "pseudo")) {
   			pseudoFeatureIds.Add(fnames[i]->idx);
   		}
   	 }
  	 GffNameList& attrnames = gffr->names->attrs;
  	 for (int i=0;i<attrnames.Count();i++) {
  		char* n=attrnames[i]->name;
  		char* p=strifind(n, "pseudo");
  		if (p==n || (p==n+2 && tolower(n[0])=='i' && tolower(n[1])=='s')) {
  			pseudoAttrIds.Add(attrnames[i]->idx);
  		}
  	}
   }

  //int redundant=0; //redundant annotation discarded
  if (verbose) GMessage("   .. loaded %d genomic features from %s\n", gffr->gflst.Count(), fname.chars());
  //int rna_deleted=0;
  //add to GenomicSeqData, adding to existing loci and identifying intron-chain duplicates
  for (int k=0;k<gffr->gflst.Count();k++) {
     GffObj* m=gffr->gflst[k];
     if (strcmp(m->getFeatureName(), "locus")==0 && 
          m->getAttr("transcripts")!=NULL) {
        continue; //discard locus meta-features
        }
     if (this->noPseudo) {
    	 bool is_pseudo=false;
    	 for (int i=0;i<pseudoFeatureIds.Count();++i) {
    		 if (pseudoFeatureIds[i]==m->ftype_id) {
    			 is_pseudo=true;
    			 break;
    		 }
    	 }
    	 if (is_pseudo) continue;
    	 for (int i=0;i<pseudoAttrIds.Count();++i) {
    		 char* attrv=NULL;
    		 if (m->attrs!=NULL) attrv=m->attrs->getAttr(pseudoAttrIds[i]);
    		 if (attrv!=NULL) {
    			 char fc=tolower(attrv[0]);
    			 if (fc=='t' || fc=='y' || fc=='1') {
    				 is_pseudo=true;
    				 break;
    			 }
    		 }
    	 }
    	 if (is_pseudo) continue;
     }
     char* rloc=m->getAttr("locus");
     if (rloc!=NULL && startsWith(rloc, "RLOC_")) {
        m->removeAttr("locus", rloc);
        }
    /*
     if (m->exons.Count()==0 && m->children.Count()==0) {
       //a non-mRNA feature with no subfeatures
       //add a dummy exon just to have the generic exon checking work
       m->addExon(m->start,m->end);
       }
     */
     if (forceExons) {  // && m->children.Count()==0) {
       m->exon_ftype_id=gff_fid_exon;
       }
     GList<GffObj> gfadd(false,false);
     if (gf_validate!=NULL && !(*gf_validate)(m, &gfadd)) {
       continue;
       }
     m->isUsed(true); //so the gffreader won't destroy it
     int i=-1;
     GenomicSeqData f(m->gseq_id);
     GenomicSeqData* gdata=NULL;
     if (seqdata.Found(&f,i)) gdata=seqdata[i];
         else { //entry not created yet for this genomic seq
           gdata=new GenomicSeqData(m->gseq_id);
           seqdata.Add(gdata);
           }
    for (int k=0;k<gfadd.Count();k++) {
      placeGf(gfadd[k], gdata, doCluster, doCollapseRedundant, matchAllIntrons, fuzzSpan);
      }
    placeGf(m, gdata, doCluster, doCollapseRedundant, matchAllIntrons, fuzzSpan);
    } //for each read gffObj
   //if (verbose) GMessage("  .. %d records from %s clustered into loci.\n", gffr->gflst.Count(), fname.chars());
   if (f!=stdin) { fclose(f); f=NULL; }
   delete gffr;
}
