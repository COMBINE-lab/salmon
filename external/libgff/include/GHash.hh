/********************************************************************************
*                  Hash table class template (char* based)                               *
*********************************************************************************/

#ifndef GHash_HH
#define GHash_HH
#include "GBase.h"

/**
* This class maintains a fast-access hash table of entities
* indexed by a character string (essentially, maps strings to pointers)
*/


template <class OBJ> class GHash {
 protected:
	struct GHashEntry {
	     char*   key;              // Key string
	     bool    keyalloc;         //shared key flag (to not free the key chars)
	     int     hash;             // Hash value of key
	     pointer data;              // Data
	     bool    mark;             // Entry is marked
	     };
  GHashEntry* hash;         // Hash
  int         fCapacity;     // table size
  int         fCount;        // number of valid entries
  int  fCurrentEntry;
  char* lastkeyptr; //pointer to last key string added
    //---------- Raw data retrieval (including empty entries
  // Return key at position pos.
  const char* Key(uint pos) const { return hash[pos].key; }
  // return data OBJ* at given position
  OBJ* Data(uint pos) const { return (OBJ*) hash[pos].data; }
  // Return mark flag of entry at position pos.
  bool Mark(uint pos) const { return hash[pos].mark; }
  // Return position of first filled slot, or >= fCapacity
  int First() const;
  // Return position of last filled slot or -1
  int Last() const;
  // Return position of next filled slot in hash table
  // or a value greater than or equal to fCapacity if no filled
  // slot was found
  int Next(int pos) const;
  //Return position of previous filled slot in hash table
  //or a -1 if no filled slot was found
  int Prev(int pos) const;

private:
  GHash(const GHash&);
  GHash &operator=(const GHash&);
  GFreeProc* fFreeProc; //procedure to free item data
protected:
public:
  static void DefaultFreeProc(pointer item) {
      delete (OBJ*)item;
      }
public:
  GHash(GFreeProc* freeProc); // constructs of an empty hash
  GHash(bool doFree=true); // constructs of an empty hash (free the item objects)
  void setFreeItem(GFreeProc *freeProc) { fFreeProc=freeProc; }
  void setFreeItem(bool doFree) { fFreeProc=(doFree)? &DefaultFreeProc : NULL; }
  int Capacity() const { return fCapacity; } // table's size, including the empty slots.
  void Resize(int m);  // Resize the table to the given size.
  int Count() const { return fCount; }// the total number of entries in the table.
  // Insert a new entry into the table given key and mark.
  // If there is already an entry with that key, leave it unchanged,
  const OBJ* Add(const char* ky, const OBJ* ptr=NULL, bool mrk=false);
  //same as Add, but the key pointer is stored directly, no string duplicate
  //is made (shared-key-Add)
  const OBJ* shkAdd(const char* ky, const OBJ* ptr, bool mrk=false);

  // Replace data at key, if the entry's mark is less than
  // or equal to the given mark.  If there was no existing entry,
  // a new entry is inserted with the given mark.
  OBJ* Replace(const char* ky, const OBJ* ptr, bool mrk=false);
  // Remove a given key and its data
  OBJ* Remove(const char* ky);
  // Find data OBJ* given key.
  OBJ* Find(const char* ky, char** keyptr=NULL);
  bool hasKey(const char* ky);
  char* getLastKey() { return lastkeyptr; }
  OBJ* operator[](const char* ky) { return Find(ky); }
  void startIterate(); //iterator-like initialization
  char* NextKey(); //returns next valid key in the table (NULL if no more)
  OBJ* NextData(); //returns next valid hash[].data
  OBJ* NextData(char*& nextkey); //returns next valid hash[].data
                                //or NULL if no more
                                //nextkey is SET to the corresponding key
  GHashEntry* NextEntry() { //returns a pointer to a GHashEntry
  	 register int pos=fCurrentEntry;
  	 while (pos<fCapacity && hash[pos].hash<0) pos++;
  	 if (pos==fCapacity) {
  	                 fCurrentEntry=fCapacity;
  	                 return NULL;
  	                 }
  	              else {
  	                 fCurrentEntry=pos+1;
  	                 return &hash[pos];
  	                 }
  }
  /// Clear all entries
  void Clear();

  /// Destructor
  virtual ~GHash();
  };
//
//======================== method definitions ========================
//
/*
  Notes:
  - The hash algorithm should yield a fCount in the range [0...GHash::EMPTY)
     GHash::EMPTY and GHash::UNUSED are needed for flag purposes.
  - Since the algorithm doubles the table size when exceeding MAX_LOAD,
    it would be prudent to keep MIN_LOAD less than 1/2 MAX_LOAD;
    otherwise, the algorithm might hip-hop between halving and doubling,
    which would be quite expensive!!
  - Not many people seem to know that hash tables don't have to be prime
    numbers; in fact, a table size of 2**n and odd probe distance are very
    easy to arrange, and this works just as well!
  - We store the hash key, so that 99.999% of the time we can compare hash numbers;
    only when hash numbers match do we need to compare keys.
    Thus, with a good hash function, the fCount of calls to strcmp() should be
    roughly the same as the fCount of successful lookups.
  - The hash table should NEVER get full, or stuff will loop forever!!
*/

// Initial table size (MUST be power of 2)
#define DEF_HASH_SIZE      32
// Maximum hash table load factor (%)
#define MAX_LOAD           80
// Minimum hash table load factor (%)
#define MIN_LOAD           10
// Probe Position [0..n-1]
#define HASH1(x,n) (((unsigned int)(x)*13)%(n))
// Probe Distance [1..n-1]
#define HASH2(x,n) (1|(((unsigned int)(x)*17)%((n)-1)))

#define FREEDATA (fFreeProc!=NULL)

/*******************************************************************************/
// Construct empty hash
template <class OBJ> GHash<OBJ>::GHash(GFreeProc* freeProc) {
  GMALLOC(hash, sizeof(GHashEntry)*DEF_HASH_SIZE);
  fCurrentEntry=-1;
  fFreeProc=freeProc;
  lastkeyptr=NULL;
  for (uint i=0; i<DEF_HASH_SIZE; i++)
         hash[i].hash=-1; //this will be an indicator for 'empty' entries
  fCapacity=DEF_HASH_SIZE;
  fCount=0;
  }

template <class OBJ> GHash<OBJ>::GHash(bool doFree) {
  GMALLOC(hash, sizeof(GHashEntry)*DEF_HASH_SIZE);
  fCurrentEntry=-1;
  lastkeyptr=NULL;
  fFreeProc = (doFree)?&DefaultFreeProc : NULL;
  for (uint i=0; i<DEF_HASH_SIZE; i++)
         hash[i].hash=-1; //this will be an indicator for 'empty' entries
  fCapacity=DEF_HASH_SIZE;
  fCount=0;
  }


// Resize table
template <class OBJ> void GHash<OBJ>::Resize(int m){
  register int i,n,p,x,h;
  GHashEntry *k;
  GASSERT(fCount<=fCapacity);
  if(m<DEF_HASH_SIZE) m=DEF_HASH_SIZE;
  n=fCapacity;
  while((n>>2)>m) n>>=1;            // Shrink until n/4 <= m
  while((n>>1)<m) n<<=1;            // Grow until m <= n/2
  GASSERT(m<=(n>>1));
  GASSERT(DEF_HASH_SIZE<=n);
  if(n!=fCapacity){
    GASSERT(m<=n);
    GMALLOC(k, sizeof(GHashEntry)*n);
    for(i=0; i<n; i++) k[i].hash=-1;
    for(i=0; i<fCapacity; i++){
      h=hash[i].hash;
      if(0<=h){
        p=HASH1(h,n);
        GASSERT(0<=p && p<n);
        x=HASH2(h,n);
        GASSERT(1<=x && x<n);
        while(k[p].hash!=-1) p=(p+x)%n;
        GASSERT(k[p].hash<0);
        k[p]=hash[i];
        }
      }
    GFREE(hash);
    hash=k;
    fCapacity=n;
    }
  }

// add a new entry, or update it if it already exists
template <class OBJ> const OBJ* GHash<OBJ>::Add(const char* ky,
                      const OBJ* pdata,bool mrk){
  register int p,i,x,h,n;
  if(!ky) GError("GHash::insert: NULL key argument.\n");
  GASSERT(fCount<fCapacity);
  h=strhash(ky);
  GASSERT(0<=h);
  p=HASH1(h,fCapacity);
  GASSERT(0<=p && p<fCapacity);
  x=HASH2(h,fCapacity);
  GASSERT(1<=x && x<fCapacity);
  i=-1;
  n=fCapacity;
  while(n && hash[p].hash!=-1){
    if ((i==-1)&&(hash[p].hash==-2)) i=p;
    if (hash[p].hash==h && strcmp(hash[p].key,ky)==0) {
      //replace hash data for this key!
      lastkeyptr=hash[p].key;
      hash[p].data = (void*) pdata;
      return (OBJ*)hash[p].data;
      }
    p=(p+x)%fCapacity;
    n--;
    }
  if(i==-1) i=p;
  GTRACE(("GHash::insert: key=\"%s\"\n",ky));
  //GMessage("GHash::insert: key=\"%s\"\n",ky);
  GASSERT(0<=i && i<fCapacity);
  GASSERT(hash[i].hash<0);
  hash[i].hash=h;
  hash[i].mark=mrk;
  hash[i].key=Gstrdup(ky);
  hash[i].keyalloc=true;
  lastkeyptr=hash[i].key;
  hash[i].data= (void*) pdata;
  fCount++;
  if((100*fCount)>=(MAX_LOAD*fCapacity)) Resize(fCount);
  GASSERT(fCount<fCapacity);
  return pdata;
  }

template <class OBJ> const OBJ* GHash<OBJ>::shkAdd(const char* ky,
                      const OBJ* pdata,bool mrk){
  register int p,i,x,h,n;
  if(!ky) GError("GHash::insert: NULL key argument.\n");
  GASSERT(fCount<fCapacity);
  h=strhash(ky);
  GASSERT(0<=h);
  p=HASH1(h,fCapacity);
  GASSERT(0<=p && p<fCapacity);
  x=HASH2(h,fCapacity);
  GASSERT(1<=x && x<fCapacity);
  i=-1;
  n=fCapacity;
  while(n && hash[p].hash!=-1){
    if((i==-1)&&(hash[p].hash==-2)) i=p;
    if(hash[p].hash==h && strcmp(hash[p].key,ky)==0){
      //replace hash data for this key!
      lastkeyptr=hash[p].key;
      hash[p].data = (void*) pdata;
      return (OBJ*)hash[p].data;
      }
    p=(p+x)%fCapacity;
    n--;
    }
  if(i==-1) i=p;
  GTRACE(("GHash::insert: key=\"%s\"\n",ky));
  //GMessage("GHash::insert: key=\"%s\"\n",ky);
  GASSERT(0<=i && i<fCapacity);
  GASSERT(hash[i].hash<0);
  hash[i].hash=h;
  hash[i].mark=mrk;
  hash[i].key=(char *)ky;
  lastkeyptr=hash[i].key;
  hash[i].keyalloc=false;
  hash[i].data= (void*) pdata;
  fCount++;
  if((100*fCount)>=(MAX_LOAD*fCapacity)) Resize(fCount);
  GASSERT(fCount<fCapacity);
  return pdata;
  }


// Add or replace entry
template <class OBJ>  OBJ* GHash<OBJ>::Replace(const char* ky,const OBJ* pdata, bool mrk){
  register int p,i,x,h,n;
  if(!ky){ GError("GHash::replace: NULL key argument.\n"); }
  GASSERT(fCount<fCapacity);
  h=strhash(ky);
  GASSERT(0<=h);
  p=HASH1(h,fCapacity);
  GASSERT(0<=p && p<fCapacity);
  x=HASH2(h,fCapacity);
  GASSERT(1<=x && x<fCapacity);
  i=-1;
  n=fCapacity;
  while(n && hash[p].hash!=-1){
    if((i==-1)&&(hash[p].hash==-2)) i=p;
    if(hash[p].hash==h && strcmp(hash[p].key,ky)==0){
      if(hash[p].mark<=mrk){
        GTRACE(("GHash::replace: %08x: replacing: \"%s\"\n",this,ky));
        if (FREEDATA) (*fFreeProc)(hash[p].data);
        hash[p].mark=mrk;
        hash[p].data=pdata;
        }
      return hash[p].data;
      }
    p=(p+x)%fCapacity;
    n--;
    }
  if(i==-1) i=p;
  GTRACE(("GHash::replace: %08x: inserting: \"%s\"\n",this,ky));
  GASSERT(0<=i && i<fCapacity);
  GASSERT(hash[i].hash<0);
  hash[i].hash=h;
  hash[i].mark=mrk;
  hash[i].key=Gstrdup(ky);
  hash[i].data=pdata;
  fCount++;
  if((100*fCount)>=(MAX_LOAD*fCapacity)) Resize(fCount);
  GASSERT(fCount<fCapacity);
  return pdata;
  }


// Remove entry
template <class OBJ> OBJ* GHash<OBJ>::Remove(const char* ky){
  register int p,x,h,n;
  if(!ky){ GError("GHash::remove: NULL key argument.\n"); }
  if(0<fCount){
    h=strhash(ky);
    GASSERT(0<=h);
    p=HASH1(h,fCapacity);
    GASSERT(0<=p && p<fCapacity);
    x=HASH2(h,fCapacity);
    GASSERT(1<=x && x<fCapacity);
    GASSERT(fCount<fCapacity);
    n=fCapacity;
    while(n && hash[p].hash!=-1){
      if(hash[p].hash==h && strcmp(hash[p].key,ky)==0){
        GTRACE(("GHash::remove: %08x removing: \"%s\"\n",this,ky));
        hash[p].hash=-2;
        hash[p].mark=false;
        if (hash[p].keyalloc) GFREE((hash[p].key));
        if (FREEDATA) (*fFreeProc)(hash[p].data);
        hash[p].key=NULL;
        hash[p].data=NULL;
        fCount--;
        if((100*fCount)<=(MIN_LOAD*fCapacity)) Resize(fCount);
        GASSERT(fCount<fCapacity);
        return NULL;
        }
      p=(p+x)%fCapacity;
      n--;
      }
    }
  return NULL;
  }


// Find entry
template <class OBJ> bool GHash<OBJ>::hasKey(const char* ky) {
  register int p,x,h,n;
  if(!ky){ GError("GHash::find: NULL key argument.\n"); }
  if(0<fCount){
    h=strhash(ky);
    GASSERT(0<=h);
    p=HASH1(h,fCapacity);
    GASSERT(0<=p && p<fCapacity);
    x=HASH2(h,fCapacity);
    GASSERT(1<=x && x<fCapacity);
    GASSERT(fCount<fCapacity);
    n=fCapacity;
    while(n && hash[p].hash!=-1){
      if(hash[p].hash==h && strcmp(hash[p].key,ky)==0){
        return true;
        }
      p=(p+x)%fCapacity;
      n--;
      }
    }
  return false;
}

template <class OBJ> OBJ* GHash<OBJ>::Find(const char* ky, char** keyptr){
  register int p,x,h,n;
  if(!ky){ GError("GHash::find: NULL key argument.\n"); }
  if(0<fCount){
    h=strhash(ky);
    GASSERT(0<=h);
    p=HASH1(h,fCapacity);
    GASSERT(0<=p && p<fCapacity);
    x=HASH2(h,fCapacity);
    GASSERT(1<=x && x<fCapacity);
    GASSERT(fCount<fCapacity);
    n=fCapacity;
    while(n && hash[p].hash!=-1){
      if(hash[p].hash==h && strcmp(hash[p].key,ky)==0){
        if (keyptr!=NULL) *keyptr = hash[p].key;
        return (OBJ*)hash[p].data;
        }
      p=(p+x)%fCapacity;
      n--;
      }
    }
  return NULL;
  }


template <class OBJ> void GHash<OBJ>::startIterate() {// initialize a key iterator; call
 fCurrentEntry=0;
}

template <class OBJ> char* GHash<OBJ>::NextKey() {
 register int pos=fCurrentEntry;
 while (pos<fCapacity && hash[pos].hash<0) pos++;
 if (pos==fCapacity) {
                 fCurrentEntry=fCapacity;
                 return NULL;
                 }
              else {
                 fCurrentEntry=pos+1;
                 return hash[pos].key;
                 }
}

template <class OBJ> OBJ* GHash<OBJ>::NextData() {
 register int pos=fCurrentEntry;
 while (pos<fCapacity && hash[pos].hash<0) pos++;
 if (pos==fCapacity) {
                 fCurrentEntry=fCapacity;
                 return NULL;
                 }
              else {
                 fCurrentEntry=pos+1;
                 return (OBJ*)hash[pos].data;
                 }

}

template <class OBJ> OBJ* GHash<OBJ>::NextData(char* &nextkey) {
 register int pos=fCurrentEntry;
 while (pos<fCapacity && hash[pos].hash<0) pos++;
 if (pos==fCapacity) {
                 fCurrentEntry=fCapacity;
                 nextkey=NULL;
                 return NULL;
                 }
              else {
                 fCurrentEntry=pos+1;
                 nextkey=hash[pos].key;
                 return (OBJ*)hash[pos].data;
                 }

}


// Get first non-empty entry
template <class OBJ> int GHash<OBJ>::First() const {
  register int pos=0;
  while(pos<fCapacity){ if(0<=hash[pos].hash) break; pos++; }
  GASSERT(fCapacity<=pos || 0<=hash[pos].hash);
  return pos;
  }

// Get last non-empty entry
template <class OBJ> int GHash<OBJ>::Last() const {
  register int pos=fCapacity-1;
  while(0<=pos){ if(0<=hash[pos].hash) break; pos--; }
  GASSERT(pos<0 || 0<=hash[pos].hash);
  return pos;
  }


// Find next valid entry
template <class OBJ> int GHash<OBJ>::Next(int pos) const {
  GASSERT(0<=pos && pos<fCapacity);
  while(++pos <= fCapacity-1){ if(0<=hash[pos].hash) break; }
  GASSERT(fCapacity<=pos || 0<=hash[pos].hash);
  return pos;
  }


// Find previous valid entry
template <class OBJ> int GHash<OBJ>::Prev(int pos) const {
  GASSERT(0<=pos && pos<fCapacity);
  while(--pos >= 0){ if(0<=hash[pos].hash) break; }
  GASSERT(pos<0 || 0<=hash[pos].hash);
  return pos;
  }


// Remove all
template <class OBJ> void GHash<OBJ>::Clear(){
  register int i;
  for(i=0; i<fCapacity; i++){
    if(hash[i].hash>=0){
      if (hash[i].keyalloc) GFREE((hash[i].key));
      if (FREEDATA)
            (*fFreeProc)(hash[i].data);
      }
    }
  GFREE(hash);
  GMALLOC(hash, sizeof(GHashEntry)*DEF_HASH_SIZE);
  //reinitialize it
  for (i=0; i<DEF_HASH_SIZE; i++)
         hash[i].hash=-1; //this will be an indicator for 'empty' entries
  fCapacity=DEF_HASH_SIZE;
  fCount=0;
  }


// Save data
/*
void GHash::Save(Stream& store) const {
  Object::save(store);
  store << fCapacity;
  store << fCount;
  for(int i=0; i<fCapacity; i++){
    store << hash[i].hash;
    if(hash[i].hash>=0){
      uint len=strlen(hash[i].key);
      store << len;
      store << hash[i].mark;
      store.save(hash[i].key,len);
      }
    }
  }


// Load data
void GHash::Load(Stream& store){
  Object::load(store);
  store >> fCapacity;
  store >> fCount;
  for(int i=0; i<fCapacity; i++){
    store >> hash[i].hash;
    if(hash[i].hash>=0){
      uint len;
      store >> len;
      store >> hash[i].mark;
      GMALLOC(hash[i].key,len+1);
      store.load(hash[i].key,len);
      hash[i].key[len]='\0';
      }
    }
  }
*/

// Destroy table
template <class OBJ> GHash<OBJ>::~GHash(){
  register int i;
  for(i=0; i<fCapacity; i++){
    if(hash[i].hash>=0){
      if (hash[i].keyalloc) GFREE((hash[i].key));
      if (FREEDATA) (*fFreeProc)(hash[i].data);
      }
    }
  GFREE(hash);
  }

#endif
