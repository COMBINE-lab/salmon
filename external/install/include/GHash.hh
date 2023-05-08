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

//#define HASH_DBG_PRINT 1

#define GSTR_HASH(s) strhash(s)
//#define GSTR_HASH(s) djb_hash(s)
//#define GSTR_HASH(s) fnv1a_hash(s)
//#define GSTR_HASH(s) murmur3(s)

template <class OBJ> class GHash {
protected:
	struct GHashEntry {
		char*   key;              // Key string
		bool    keyalloc;         // shared key flag (to free/not the key)
		int     hash;             // Hash value of key
		pointer data;             // Data
	};
	GHashEntry* hash;          // Hash
	int         fCapacity;     // table size
	int         fCount;        // number of valid entries
	int  fCurrentEntry;
	char* lastkeyptr; //pointer to last key string added
	//---------- Raw data retrieval (including empty entries)
	// Return key at position pos.
	const char* Key(uint pos) const { return hash[pos].key; }
	// return data OBJ* at given position
	OBJ* Data(uint pos) const { return (OBJ*) hash[pos].data; }
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

	// Insert a new entry into the table given key.
	// If there is already an entry with that key, leave it unchanged
	OBJ* Add(const char* ky, OBJ* ptr=NULL);

	//same with Add, but frees the old element if it's a replacement
	OBJ* fAdd(const char* ky, OBJ* ptr=NULL);

	//same as Add, but the key pointer is stored directly, no string copy needed
	//(shared-key-Add)
	OBJ* shkAdd(const char* ky, OBJ* ptr);

	// Replace data at key. If there was no existing entry,
	// a new entry is inserted.
	OBJ* Replace(const char* ky, OBJ* ptr);
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
		int pos=fCurrentEntry;
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
  - Since the algorithm doubles the table size when exceeding MAX_LOAD,
    it would be prudent to keep MIN_LOAD less than 1/2 MAX_LOAD;
    otherwise, the algorithm might flip between halving and doubling!
  - We store the key hash value so that 99.999% of the time we can compare hash numbers;
    only when hash numbers match we need to compare keys.
  - Thus with a good hash function the fCount of calls to strcmp() should be
    roughly the same as the fCount of successful lookups.
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
template <class OBJ> void GHash<OBJ>::Resize(int m) {
	int i,n,p,x,h;
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
			if(h>=0){
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


template <class OBJ> OBJ* GHash<OBJ>::Add(const char* ky, OBJ* pdata) {
	int p,i,x,h,n;
	if(!ky) GError("GHash::insert: NULL key argument.\n");
	GASSERT(fCount<fCapacity);
	h=GSTR_HASH(ky);
	GASSERT(0<=h);
	p=HASH1(h,fCapacity);
	GASSERT(0<=p && p<fCapacity);
	x=HASH2(h,fCapacity);
	GASSERT(1<=x && x<fCapacity);
	i=-1;
	n=fCapacity;
#ifdef HASH_DBG_PRINT
	int iterations=0;
	int init_p=p;
	int init_x=x;
#endif
	while(n && hash[p].hash!=-1) {
		if ((i==-1)&&(hash[p].hash==-2)) i=p;
		if (hash[p].hash==h && strcmp(hash[p].key,ky)==0) {
			//replace hash data for this key!
			lastkeyptr=hash[p].key;
			OBJ* r = (OBJ*) hash[p].data;
			hash[p].data = (void*) pdata;
#ifdef HASH_DBG_PRINT
			GMessage("Add.R\t%s\t%d,%d,%d\t%d\t%d\t%d\n",
					ky, h,init_p,init_x, iterations,  fCount, fCapacity);
#endif
			return r;
		}
		p=(p+x)%fCapacity;
		n--;
	}
	if(i==-1) i=p;
#ifdef HASH_DBG_PRINT
	GMessage("Add.N\t%s\t%d,%d,%d\t%d\t%d\t%d\n",
			ky, h,init_p,init_x, iterations,  fCount, fCapacity);
#endif
	GTRACE(("GHash::insert: key=\"%s\"\n",ky));
	//GMessage("GHash::insert: key=\"%s\"\n",ky);
	GASSERT(0<=i && i<fCapacity);
	GASSERT(hash[i].hash<0);
	hash[i].hash=h;
	hash[i].key=Gstrdup(ky);
	hash[i].keyalloc=true;
	lastkeyptr=hash[i].key;
	hash[i].data= (void*) pdata;
	fCount++;
	if((100*fCount)>=(MAX_LOAD*fCapacity)) Resize(fCount);
	GASSERT(fCount<fCapacity);
	return pdata;
}

template <class OBJ> OBJ* GHash<OBJ>::fAdd(const char* ky, OBJ* pdata) {
	int p,i,x,h,n;
	if(!ky) GError("GHash::insert: NULL key argument.\n");
	GASSERT(fCount<fCapacity);
	h=GSTR_HASH(ky);
	GASSERT(0<=h);
	p=HASH1(h,fCapacity);
	GASSERT(0<=p && p<fCapacity);
	x=HASH2(h,fCapacity);
	GASSERT(1<=x && x<fCapacity);
	i=-1;
	n=fCapacity;
#ifdef HASH_DBG_PRINT
	int iterations=0;
	int init_p=p;
	int init_x=x;
#endif
	while(n && hash[p].hash!=-1) {
		if ((i==-1)&&(hash[p].hash==-2)) i=p;
		if (hash[p].hash==h && strcmp(hash[p].key,ky)==0) {
			//replace hash data for this key!
			lastkeyptr=hash[p].key;
			if (FREEDATA) (*fFreeProc)(hash[p].data);
			hash[p].data = (void*) pdata;
#ifdef HASH_DBG_PRINT
			GMessage("Add.R\t%s\t%d,%d,%d\t%d\t%d\t%d\n",
					ky, h,init_p,init_x, iterations,  fCount, fCapacity);
#endif
			return pdata;
		}
		p=(p+x)%fCapacity;
#ifdef HASH_DBG_PRINT
		++iterations;
#endif
		n--;
	}
	if(i==-1) i=p;
#ifdef HASH_DBG_PRINT
	GMessage("Add.N\t%s\t%d,%d,%d\t%d\t%d\t%d\n",
			ky, h,init_p,init_x, iterations,  fCount, fCapacity);
#endif
	GTRACE(("GHash::insert: key=\"%s\"\n",ky));
	//GMessage("GHash::insert: key=\"%s\"\n",ky);
	GASSERT(0<=i && i<fCapacity);
	GASSERT(hash[i].hash<0);
	hash[i].hash=h;
	hash[i].key=Gstrdup(ky);
	hash[i].keyalloc=true;
	lastkeyptr=hash[i].key;
	hash[i].data= (void*) pdata;
	fCount++;
	if((100*fCount)>=(MAX_LOAD*fCapacity)) Resize(fCount);
	GASSERT(fCount<fCapacity);
	return pdata;
}

template <class OBJ> OBJ* GHash<OBJ>::shkAdd(const char* ky, OBJ* pdata) {
	int p,i,x,h,n;
	if(!ky) GError("GHash::insert: NULL key argument.\n");
	GASSERT(fCount<fCapacity);
	h=GSTR_HASH(ky);
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
template <class OBJ>  OBJ* GHash<OBJ>::Replace(const char* ky, OBJ* pdata){
	int p,i,x,h,n;
	if(!ky){ GError("GHash::replace: NULL key argument.\n"); }
	GASSERT(fCount<fCapacity);
	h=GSTR_HASH(ky);
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
			GTRACE(("GHash::replace: %08x: replacing: \"%s\"\n",this,ky));
			if (FREEDATA) (*fFreeProc)(hash[p].data);
			hash[p].data=pdata;
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
	hash[i].key=Gstrdup(ky);
	hash[i].data=pdata;
	fCount++;
	if((100*fCount)>=(MAX_LOAD*fCapacity)) Resize(fCount);
	GASSERT(fCount<fCapacity);
	return pdata;
}


// Remove entry
template <class OBJ> OBJ* GHash<OBJ>::Remove(const char* ky){
	int p,x,h,n;
	if(!ky){ GError("GHash::remove: NULL key argument.\n"); }
	OBJ* removed=NULL;
	if(0<fCount){
		h=GSTR_HASH(ky);
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
				if (hash[p].keyalloc) GFREE((hash[p].key));
				if (FREEDATA) (*fFreeProc)(hash[p].data);
				else removed=(OBJ*)hash[p].data;
				hash[p].key=NULL;
				hash[p].data=NULL;
				fCount--;
				if((100*fCount)<=(MIN_LOAD*fCapacity)) Resize(fCount);
				GASSERT(fCount<fCapacity);
				return removed;
			}
			p=(p+x)%fCapacity;
			n--;
		}
	}
	return removed;
}


// Find entry
template <class OBJ> bool GHash<OBJ>::hasKey(const char* ky) {
	int p,x,h,n;
	if(!ky){ GError("GHash::find: NULL key argument.\n"); }
	if(0<fCount){
		h=GSTR_HASH(ky);
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
	int p,x,h,n;
	if(!ky){ GError("GHash::find: NULL key argument.\n"); }
	if (fCount==0) return NULL;
	h=GSTR_HASH(ky);
	GASSERT(0<=h);
	p=HASH1(h,fCapacity);
	GASSERT(0<=p && p<fCapacity);
	x=HASH2(h,fCapacity);
	GASSERT(1<=x && x<fCapacity);
	GASSERT(fCount<fCapacity);
	n=fCapacity;
#ifdef HASH_DBG_PRINT
	int iterations=0;
	int init_p=p;
	int init_x=x;
#endif
	while(n && hash[p].hash!=-1){
		if(hash[p].hash==h && strcmp(hash[p].key,ky)==0){
			if (keyptr!=NULL) *keyptr = hash[p].key;
#ifdef HASH_DBG_PRINT
			GMessage("Found \t%s\t%d,%d,%d\t%d\t%d\t%d\n",
					ky, h,init_p,init_x, iterations,  fCount, fCapacity);
#endif
			return (OBJ*)hash[p].data;
		}
		p=(p+x)%fCapacity;
		n--;
#ifdef HASH_DBG_PRINT
		++iterations;
#endif
	}
#ifdef HASH_DBG_PRINT
	GMessage("Nfound\t%s\t%d,%d,%d\t%d\t%d\t%d\n",
			ky, h,init_p,init_x, iterations,  fCount, fCapacity);
#endif
	return NULL;
}

template <class OBJ> void GHash<OBJ>::startIterate() {// initialize a key iterator; call
	fCurrentEntry=0;
}

template <class OBJ> char* GHash<OBJ>::NextKey() {
	int pos=fCurrentEntry;
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
	int pos=fCurrentEntry;
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
	int pos=fCurrentEntry;
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
	int pos=0;
	while(pos<fCapacity){ if(0<=hash[pos].hash) break; pos++; }
	GASSERT(fCapacity<=pos || 0<=hash[pos].hash);
	return pos;
}

// Get last non-empty entry
template <class OBJ> int GHash<OBJ>::Last() const {
	int pos=fCapacity-1;
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
	int i;
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

// Destroy table
template <class OBJ> GHash<OBJ>::~GHash(){
	for(int i=0; i<fCapacity; i++){
		if(hash[i].hash>=0){
			if (hash[i].keyalloc) GFREE((hash[i].key));
			if (FREEDATA) (*fFreeProc)(hash[i].data);
		}
	}
	GFREE(hash);
}

class GStrSet:public GHash<int> {
protected:
	bool free_keys;
public:
	GStrSet(bool shared_keys=false):GHash<int>(false), free_keys(!shared_keys) {
	}
	void Add(const char* str) {
		if (free_keys) {
			//allocates a copy of str
			GHash<int>::Add(str, NULL);
		}
		else this->shkAdd(str, NULL);
	}
	void add(const char* str) { this->Add(str); }
	void push(const char* str) { this->Add(str); }
	bool has(const char* str) {
		return hasKey(str);
	}

};

#endif
