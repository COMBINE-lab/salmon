//---------------------------------------------------------------------------
/*
Sortable collections of objects and object pointers
*/
#ifndef _GList_HH
#define _GList_HH

#include "GVec.hh"

#define GLIST_SORTED_ERR "Operation not allowed on a sorted list!\n"
#define GLIST_UNSORTED_ERR "Operation not allowed on an unsorted list!\n"

//------ useful macros:
#define BE_UNSORTED if (fCompareProc!=NULL) { GError(GLIST_SORTED_ERR); return; }
#define BE_SORTED if (fCompareProc==NULL) { GError(GLIST_UNSORTED_ERR); return; }

#define SORTED (fCompareProc!=NULL)
#define UNSORTED (fCompareProc==NULL)

// GArray is the sortable array type, requires the comparison operator < to be defined
template <class OBJ> class GArray:public GVec<OBJ> {
  protected:
    bool fUnique;
    static int DefaultCompareProc(const pointer item1, const pointer item2) {
      //operator< MUST be defined for OBJ class!
      if (*((OBJ*)item2) < *((OBJ*)item1)) return 1;
        else if (*((OBJ*)item1) < *((OBJ*)item2)) return -1;
                                             else return  0;
      }
    GCompareProc* fCompareProc;
  public:
    GArray(GCompareProc* cmpFunc=NULL);
    GArray(bool sorted, bool unique=false);
    GArray(int init_capacity, bool sorted, bool unique=false);
    GArray(GArray<OBJ>& array); //copy constructor
    const GArray<OBJ>& operator=(GArray<OBJ>& array);
    //~GArray();
    //assignment operator
    void setSorted(GCompareProc* cmpFunc);
    void setSorted(bool sorted) {
     if (sorted) {
         if (fCompareProc!=&DefaultCompareProc) {
             fCompareProc=&DefaultCompareProc;
             Sort();
             }
          }
      else fCompareProc=NULL;
      }
    //sort the array if cmpFunc not NULL or changes
    int Add(OBJ* item); // specific implementation if sorted
    int Add(OBJ& item) { return Add(&item); } //both will CREATE a new OBJ and COPY to it
                       // using OBJ new operator=
    int AddIfNew(OBJ& item, int* fidx=NULL); //requires == operator
        //if equal item not found, item is added and return the index of it
        //otherwise returns -1 and fidx is set to the equal item location
    int cAdd(OBJ item) { return Add(&item); }
    int cPush(OBJ item) { return Add(&item); }
    int Push(OBJ& item) { return Add(&item); }

    void Add(GArray<OBJ>& list); //add copies of all items from another list
    //this will reject identical items in sorted lists only!
    void setUnique(bool beUnique) { fUnique = beUnique; };
    void Sort(); //explicit sort may be requested
    bool Sorted() { return fCompareProc!=NULL; }
    void Replace(int idx, OBJ& item); //Put, use operator= to copy
    int  Unique() { return fUnique; }
    int IndexOf(OBJ& item);
         //this needs the == operator to have been defined for OBJ
    bool Found(OBJ& item, int& idx); // for sorted arrays only;
         //search by content; if found, returns true and idx will be the index
         //of the first item found matching for which fCompareProc returns 0
    bool Exists(OBJ& item); //same as above without existing index info
    //unsorted only, place item at position idx:
    void Move(int curidx, int newidx);
    void Insert(int idx, OBJ* item);
    void Insert(int idx, OBJ item) { Insert(idx,&item); }
};

//GList is a sortable collection of pointers to objects; requires operator< to be defined, or a custom compare function
template <class OBJ> class GList:public GPVec<OBJ> {
  protected:
    bool fUnique;
    GCompareProc* fCompareProc; //a pointer to a Compare function

    static int DefaultCompareProc(const pointer item1, const pointer item2) {
      //operator< MUST be defined for OBJ class!
      if (*((OBJ*)item2) < *((OBJ*)item1)) return 1;
        else if (*((OBJ*)item1) < *((OBJ*)item2)) return -1;
                                             else return  0;
      }
  public:
    void sortInsert(int idx, OBJ* item); //special insert in sorted lists
         //WARNING: the caller must know the insert index such that the sort order is preserved!
    GList(GCompareProc* compareProc=NULL); //free by default
    GList(GCompareProc* compareProc, //unsorted by default
        GFreeProc *freeProc,
        bool beUnique=false);
    GList(bool sorted, bool free_elements=true, bool beUnique=false);
    GList(int init_capacity, bool sorted, bool free_elements=true, bool beUnique=false);
    GList(GList<OBJ>& list); //copy constructor?
    GList(GList<OBJ>* list); //kind of a copy constructor
    const GList<OBJ>& operator=(GList<OBJ>& list);
    //void Clear();
    //~GList();
    void setSorted(GCompareProc* compareProc);
       //sorted if compareProc not NULL; sort the list if compareProc changes !
    bool Sorted() { return fCompareProc!=NULL; }
    void setSorted(bool sorted) {
     if (sorted) {
         if (fCompareProc!=&DefaultCompareProc) {
             fCompareProc=&DefaultCompareProc;
             Sort();
             }
          }
      else fCompareProc=NULL;
      }
    int Add(OBJ* item); //-- specific implementation if sorted - may become an Insert()
    void Add(GList<OBJ>& list); //add all pointers from another list

    OBJ* AddIfNew(OBJ* item, bool deleteIfFound=true, int* fidx=NULL);
    // default: delete item if Found() (and pointers are not equal)!
    //returns the equal (==) object if it's in the list already
    //or the item itself if it is unique and actually added

    int AddedIfNew(OBJ* item);
    // if Found(item) (and pointers are not equal) delete item and returns -1
    // if added, returns the new item index


    int Unique() { return fUnique; }
    //this will reject identical items in sorted lists only!
    void setUnique(bool beUnique) { fUnique = beUnique; };

    GCompareProc* GetCompareProc() {return fCompareProc;}
    int IndexOf(OBJ* item); //this has a specific implementation for sorted lists
               //if list is sorted, item data is located by binary search
               //based on the Compare function
               //if not, a linear search is performed, but
               //this needs the == operator to have been defined for OBJ

    void Put(int idx, OBJ* item, bool re_sort=false);
    bool Found(OBJ* item, int & idx); // sorted only;
               //search by content; if found, returns true and idx will be the index
               //of the first item found matching for which GTCompareProc returns 0
    bool Exists(OBJ* item); //same as above without existing index info
    bool Exists(OBJ& item); //same as above without existing index info
    void Sort(); //explicit sort may be requested using this function
    int Remove(OBJ* item); //search for pointer, using binary search if sorted
    void Insert(int idx, OBJ* item); //unsorted only, place item at position idx
    void Move(int curidx, int newidx);
}; //GList



//-------------------- TEMPLATE IMPLEMENTATION-------------------------------

template <class OBJ> GArray<OBJ>::GArray(GArray<OBJ>& array):GVec<OBJ>(0) { //copy constructor
 this->fCount=array.fCount;
 this->fCapacity=array.fCapacity;
 this->fArray=NULL;
 if (this->fCapacity>0) {
    //GMALLOC(this->fArray, this->fCapacity*sizeof(OBJ));
    this->fArray=new OBJ[this->fCapacity];
    }
 this->fCount=array.fCount;
 fUnique=array.fUnique;
 fCompareProc=array.fCompareProc;
 // uses OBJ operator=
 for (int i=0;i<this->fCount;i++) this->fArray[i]=array[i];
 }

template <class OBJ> const GArray<OBJ>& GArray<OBJ>::operator=(GArray<OBJ>& array) {
 if (&array==this) return *this;
 GVec<OBJ>::Clear();
 this->fCount=array.fCount;
 this->fUnique=array.fUnique;
 this->fCapacity=array.fCapacity;
 if (this->fCapacity>0) {
    //GMALLOC(this->fArray, this->fCapacity*sizeof(OBJ));
    this->fArray=new OBJ[this->fCapacity];
    }
 this->fCompareProc=array.fCompareProc;
 this->fCount=array.fCount;
 // uses OBJ operator=
 for (int i=0;i<this->fCount;i++) {
   this->fArray[i]=array[i];
   }
 return *this;
}

template <class OBJ> GArray<OBJ>::GArray(GCompareProc* cmpFunc):GVec<OBJ>(0) {
  fCompareProc = cmpFunc;
  fUnique = false; //only affects sorted lists
}

template <class OBJ> GArray<OBJ>::GArray(bool sorted, bool unique):GVec<OBJ>(0) {
  fUnique=unique;
  fCompareProc = sorted ? DefaultCompareProc : NULL;
}

template <class OBJ> GArray<OBJ>::GArray(int init_capacity,
                        bool sorted, bool unique):GVec<OBJ>(init_capacity) {
  fUnique=unique;
  fCompareProc=sorted ? DefaultCompareProc : NULL;
}

template <class OBJ> void GArray<OBJ>::setSorted(GCompareProc* cmpFunc) {
  GCompareProc* old_proc=fCompareProc;
  fCompareProc=cmpFunc;
  if (fCompareProc!=old_proc && fCompareProc!=NULL)
       Sort(); //new compare method
}

template <class OBJ> int GArray<OBJ>::IndexOf(OBJ& item) {
 int result=0;
 if (Found(item, result)) return result;
                     else return -1;
 }

template <class OBJ> bool GArray<OBJ>::Exists(OBJ& item) {
 int result=0;
 if (Found(item, result)) return true;
                     else return false;
 }


template <class OBJ> int GArray<OBJ>::Add(OBJ* item) {
 if (item==NULL) return -1;
 int result;
 if (SORTED) {
   if (Found(*item, result))
      if (fUnique) return -1; //cannot add a duplicate!
   //Found sets result to the position where the item should be!
   GVec<OBJ>::Insert(result, *item);
 }
  else {
   if (fUnique && Found(*item, result)) return -1; //set behaviour
   result = this->fCount;
   if (result==this->fCapacity) GVec<OBJ>::Grow();
   this->fArray[result] = *item; //operator=, copies the item
   this->fCount++;
 }
 return result;
}


template <class OBJ> void GArray<OBJ>::Add(GArray<OBJ>& list) {
  if (list.Count()==0) return;
  if (SORTED) {
    for (int i=0;i<list.fCount;i++) Add(&list[i]);
    }
  else { //simply copy
    this->setCapacity(this->fCapacity+list.fCount);
    int s=this->fCount;
    for (int i=0;i<list.fCount;i++)
           this->fArray[s+i]=list.fArray[i];
    this->fCount+=list.fCount;
    }
}

//returns -1 if existing equal object exists, sets fidx to that equal item index
//or returns the index where the item was added/inserted
template <class OBJ> int GArray<OBJ>::AddIfNew(OBJ& item,
                                     int* fidx) {
 int rpos;
 if (Found(item, rpos)) {
    if (fidx) *fidx=rpos; //the position where the item should be inserted:
    return -1; //found and not added
 }
 //not found, let's insert it
 if (SORTED) {
   //Found() set result to the position where the item should be inserted
   GVec<OBJ>::Insert(rpos, item);
 } else { //simply append
	rpos = this->fCount;
	if (rpos==this->fCapacity) GVec<OBJ>::Grow();
	this->fArray[rpos] = item; //operator= copies the item
	this->fCount++;
 }
 if (fidx!=NULL) *fidx=rpos;
 return rpos;
}

template <class OBJ> bool GArray<OBJ>::Found(OBJ& item, int& idx) {
 //search the list by using fCompareProc (if defined)
 //or == operator for a non-sortable list
 //for sorted lists, even when the result is false, the idx is
 //set to the closest matching object!
 int i;
 idx=-1;
 if (this->fCount==0) { idx=0; return false;}
 if (SORTED) { //binary search based on fCompareProc
   //do the simplest tests first:
   if ((*fCompareProc)(&(this->fArray[0]),&item)>0) {
                       idx=0;
                       return false;
                       }
   if ((*fCompareProc)(&item, &(this->fArray[this->fCount-1]))>0) {
                       idx=this->fCount;
                       return false;
                       }

   int l=0;
   int h = this->fCount - 1;
   int c;
   while (l <= h) {
       i = (l + h) >> 1;
       c = (*fCompareProc)(&(this->fArray[i]), &item);
       if (c < 0)  l = i + 1;
         else {
            h = i - 1;
            if (c == 0) { //found!
                 idx=i;
                 return true;
                }
            }
       } //while
   idx = l;
   return false;
   }
 else {//not sorted: use linear search
   // needs == operator to compare user defined objects !
   i=0;
   while (i<this->fCount) {
      if (this->fArray[i]==item) { //requires operator==
         idx=i;
         return true;
         }
      i++;
      }
   return false;
   }
}

template <class OBJ> void GArray<OBJ>::Insert(int idx, OBJ* item) {
 //idx can be [0..fCount] so an item can be actually added
 BE_UNSORTED; //forbid this operation on sorted data
 GVec<OBJ>::Insert(idx, item);
}


template <class OBJ> void GArray<OBJ>::Move(int curidx, int newidx) {
 BE_UNSORTED; //cannot do this in a sorted list!
 if (curidx!=newidx || newidx>=this->fCount)
     GError(GVEC_INDEX_ERR, newidx);

 OBJ tmp=this->fArray[curidx]; //copy constructor here
 this->fArray[curidx]=this->fArray[newidx];
 this->fArray[newidx]=tmp;
}

template <class OBJ> void GArray<OBJ>::Replace(int idx, OBJ& item) {
 //TEST_INDEX(idx);
 if (idx<0 || idx>=this->fCount) GError(GVEC_INDEX_ERR, __FILE__,__LINE__, idx);
 this->fArray[idx]=item;
 if ( SORTED ) Sort(); //re-sort ! this could be very expensive, don't do it
}

template <class OBJ> void GArray<OBJ>::Sort() {
 if (fCompareProc==NULL) { fCompareProc=DefaultCompareProc; }
 if (this->fArray!=NULL && this->fCount>0)
     this->qSort(0, this->fCount-1, fCompareProc);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//*=> GList implementation -- sortable array of pointers to OBJ

template <class OBJ> GList<OBJ>::GList(GList<OBJ>& list):GPVec<OBJ>(list) { //copy constructor
 fUnique=list.fUnique;
 fCompareProc=list.fCompareProc;
}

template <class OBJ> GList<OBJ>::GList(GList<OBJ>* plist):GPVec<OBJ>(0) { //another copy constructor
 this->fCapacity=plist->fCapacity;
 this->fList=NULL;
 if (this->fCapacity>0) {
     GMALLOC(this->fList, this->fCapacity*sizeof(OBJ*));
 }
 fUnique=plist->fUnique;
 fCompareProc=plist->fCompareProc;
 this->fFreeProc=plist->fFreeProc;
 this->fCount=plist->fCount;
 memcpy(this->fList, plist->fList, this->fCount*sizeof(OBJ*));
 //for (int i=0;i<list->fCount;i++) Add(plist->Get(i));
}

template <class OBJ> void GList<OBJ>::Add(GList<OBJ>& list) {
  if (list.Count()==0) return;
  if (SORTED) {
    for (int i=0;i<list.Count();i++) Add(list[i]);
    }
  else { //simply copy
    this->setCapacity(this->fCapacity+list.fCount);
    memcpy( & (this->fList[this->fCount]), list.fList, list.fCount*sizeof(OBJ*));
    this->fCount+=list.fCount;
    }
}


template <class OBJ> GList<OBJ>::GList(GCompareProc* compareProc,
       GFreeProc* freeProc, bool beUnique) {
  fCompareProc = compareProc;
  this->fFreeProc    = freeProc;
  fUnique = beUnique; //only affects sorted lists
}

template <class OBJ> GList<OBJ>::GList(GCompareProc* compareProc) {
  fCompareProc = compareProc;
  this->fFreeProc = GPVec<OBJ>::DefaultFreeProc;
  fUnique = false; //only affects sorted lists
}

template <class OBJ> GList<OBJ>::GList(bool sorted,
    bool free_elements, bool beUnique) {
  if (sorted) {
     if (free_elements) {
        fCompareProc=&DefaultCompareProc;
        this->fFreeProc = GPVec<OBJ>::DefaultFreeProc;
        fUnique=beUnique;
        }
       else {
        fCompareProc=&DefaultCompareProc;
        this->fFreeProc=NULL;
        fUnique=beUnique;
        }
     }
   else {
     if (free_elements) {
        fCompareProc=NULL;
        this->fFreeProc=GPVec<OBJ>::DefaultFreeProc;
        fUnique=beUnique;
        }
      else {
        fCompareProc=NULL;
        this->fFreeProc=NULL;
        fUnique=beUnique;
        }
     }
}


template <class OBJ> GList<OBJ>::GList(int init_capacity, bool sorted,
    bool free_elements, bool beUnique):GPVec<OBJ>(init_capacity, free_elements) {
  if (sorted) {
      fCompareProc=&DefaultCompareProc;
      fUnique=beUnique;
      }
   else {
      fCompareProc=NULL;
      fUnique=beUnique;
      }
}

template <class OBJ> const GList<OBJ>& GList<OBJ>::operator=(GList& list) {
 if (&list!=this) {
     GPVec<OBJ>::Clear();
     fCompareProc=list.fCompareProc;
     this->fFreeProc=list.fFreeProc;
     //Attention: the object pointers are copied directly,
     //but the actual objects are NOT duplicated
     for (int i=0;i<list.Count();i++) Add(list[i]);
     }
 return *this;
}

template <class OBJ> void GList<OBJ>::setSorted(GCompareProc* compareProc) {
 GCompareProc* old_proc=fCompareProc;
 fCompareProc=compareProc;
 if (fCompareProc!=old_proc && fCompareProc!=NULL)
       Sort(); //new compare method
}

template <class OBJ> int GList<OBJ>::IndexOf(OBJ* item) {
 int result=0;
 if (Found(item, result)) return result;
                     else return -1;
 }

template <class OBJ> bool GList<OBJ>::Exists(OBJ& item) {
 int result=0;
 if (Found(&item, result)) return true;
                      else return false;
 }

template <class OBJ> bool GList<OBJ>::Exists(OBJ* item) {
 int result=0;
 if (Found(item, result)) return true;
                      else return false;
 }

template <class OBJ> int GList<OBJ>::Add(OBJ* item) {
 int result;
 if (item==NULL) return -1;
 if (SORTED) {
   if (Found(item, result))
      if (fUnique) return -1; //duplicates forbidden
   //Found sets result to the position where the item should be!
   sortInsert(result, item);
   }
  else {
   if (fUnique && Found(item,result)) return -1; //set behaviour
   result = this->fCount;
   if (result==this->fCapacity) GPVec<OBJ>::Grow();
   this->fList[result]=item;
   this->fCount++;
   }
 return result;
}

//by default, it deletes item if it an equal is found in the list!
//returns the existing equal (==) object if it's in the list already
//or returns the item itself if it's unique (and adds it)
template <class OBJ> OBJ* GList<OBJ>::AddIfNew(OBJ* item,
                                     bool deleteIfFound, int* fidx) {
 int r;
 if (Found(item, r)) {
    if (deleteIfFound && (pointer)item != (pointer)(this->fList[r])) {
       this->deallocate_item(item);
    }
    if (fidx!=NULL) *fidx=r;
    return this->fList[r]; //found
 }
 //not found:
 if (SORTED) {
   //Found() set result to the position where the item should be inserted:
   sortInsert(r, item);
 } else {
   r = this->fCount;
   if (r==this->fCapacity) GPVec<OBJ>::Grow();
   this->fList[r]=item;
   this->fCount++;
 }
 if (fidx!=NULL) *fidx=r;
 return item;
}

//if item is found already in the list DELETE it and return -1
//otherwise the item is added and its index is returned
template <class OBJ> int GList<OBJ>::AddedIfNew(OBJ* item) {
 int r;
 if (Found(item, r)) {
    if ((pointer)item != (pointer)(this->fList[r])) {
        this->deallocate_item(item);
        }
    return -1;
    }
 //not found:
 if (SORTED) {
   //Found() set r to the position where the item should be inserted:
   sortInsert(r, item);
   }
  else {
   r = this->fCount;
   if (r==this->fCapacity) GPVec<OBJ>::Grow();
   this->fList[r]=item;
   this->fCount++;
   }
 return r;
}


template <class OBJ> bool GList<OBJ>::Found(OBJ* item, int& idx) {
 //search the list by using fCompareProc (if defined)
 //or == operator for a non-sortable list
 //for sorted lists, even when the result is false, the idx is
 //set to the closest matching object!
 int i;
 idx=-1;
 if (this->fCount==0) { idx=0;return false;}
 if (SORTED) { //binary search based on fCompareProc
   //do the simple test first:

   if ((*fCompareProc)(this->fList[0],item)>0) {
                       idx=0;
                       return false;
                       }
   if ((*fCompareProc)(item, this->fList[this->fCount-1])>0) {
                       idx=this->fCount;
                       return false;
                       }

   int l, h, c;
   l = 0;
   h = this->fCount - 1;
   while (l <= h) {
       i = (l + h) >> 1;
       c = (*fCompareProc)(this->fList[i], item);
       if (c < 0)  l = i + 1;
       else {
          h = i - 1;
          if (c == 0) {
               idx=i;
               return true;
          }
       }
   } //while
   idx = l;
   return false;
   }
 else {//not sorted: use linear search
   // needs == operator to compare user defined objects !
   i=0;
   while (i<this->fCount) {
      if (*this->fList[i]==*item) {
         idx=i;
         return true;
         }
      i++;
      }
   return false;
   }
}

template <class OBJ> void GList<OBJ>::sortInsert(int idx, OBJ* item) {
 //idx must be the new position this new item must have
 //so the allowed range is [0..fCount]
 //the current fList[idx] and all the above will be shifted +1
 if (idx<0 || idx>this->fCount) GError(GVEC_INDEX_ERR, idx);
 if (this->fCount==this->fCapacity) {
    GPVec<OBJ>::Grow(idx, item);
    //expand and also copy/move data and insert the new item
    return;
    }
 //room still left, just move data around and insert the new one
 if (idx<this->fCount) //copy/move pointers only!
      memmove(&(this->fList[idx+1]), &(this->fList[idx]), (this->fCount-idx)*sizeof(OBJ*));
 this->fList[idx]=item;
 this->fCount++;
}

template <class OBJ> void GList<OBJ>::Insert(int idx, OBJ* item) {
 //idx can be [0..fCount] so an item can be actually added
 BE_UNSORTED; //cannot do that with a sorted list!
 GPVec<OBJ>::Insert(idx,item);
}

template <class OBJ> void GList<OBJ>::Move(int curidx, int newidx) {
 BE_UNSORTED; //cannot do this in a sorted list!
 GPVec<OBJ>::Move(curidx,newidx);
}

template <class OBJ> void GList<OBJ>::Put(int idx, OBJ* item, bool re_sort) {
 //WARNING: this will never free the replaced item!
 // this may BREAK the sort order unless the "re_sort" parameter is given
 if (idx<0 || idx>this->fCount) GError(GVEC_INDEX_ERR, idx);
 this->fList[idx]=item;
 if (SORTED && item!=NULL && re_sort) Sort(); //re-sort
}

template <class OBJ> int GList<OBJ>::Remove(OBJ* item) {
//removes an item if it's in our list
 int result=IndexOf(item);
 if (result>=0) GPVec<OBJ>::Delete(result);
 return result;
}

template <class OBJ> void GList<OBJ>::Sort() {
 if (fCompareProc==NULL) fCompareProc = DefaultCompareProc;
 if (this->fList!=NULL && this->fCount>0)
     this->qSort(0, this->fCount-1, fCompareProc);
}

//---------------------------------------------------------------------------
#endif
