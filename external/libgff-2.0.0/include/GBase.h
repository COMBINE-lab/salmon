#ifndef G_BASE_DEFINED
#define G_BASE_DEFINED
#define GCLIB_VERSION "0.11.9"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if defined(__WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_WIN64) || defined(__MINGW64__) || defined(__WINDOWS__)
  #ifndef _WIN32
    #define _WIN32
  #endif
  #ifndef _WIN64
    #define _WIN64
  #endif
  #define __USE_MINGW_ANSI_STDIO 1
  //#define __ISO_C_VISIBLE 1999
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdint.h>
#include <stdarg.h>

#ifdef _WIN32
  #include <windows.h>
  #include <io.h>
  #define CHPATHSEP '\\'
  #undef off_t
  #define off_t int64_t
  #ifndef popen
   #define popen _popen
  #endif
  #ifndef fseeko
		#ifdef _fseeki64
			#define fseeko(stream, offset, origin) _fseeki64(stream, offset, origin)
		#else
			#define fseeko fseek
		#endif
  #endif
 #ifndef ftello
  #ifdef _ftelli64
    #define ftello(stream) _ftelli64(stream)
  #else
    #define ftello ftell
  #endif
 #endif
 #else
  #define CHPATHSEP '/'
  #include <unistd.h>
#endif

#ifndef fseeko
 #define fseeko fseek
#endif
#ifndef ftello
 #define ftello ftell
#endif

#ifdef DEBUG
#undef NDEBUG
#define _DEBUG 1
#define _DEBUG_ 1
#endif

typedef int32_t int32;
typedef uint32_t uint32;
typedef int16_t int16;
typedef uint16_t uint16;

typedef unsigned char uchar;
typedef unsigned char byte;

#ifndef MAXUINT
#define MAXUINT ((unsigned int)-1)
#endif

#ifndef MAXINT
#define MAXINT INT_MAX
#endif

#ifndef MAX_UINT
#define MAX_UINT ((unsigned int)-1)
#endif

#ifndef MAX_INT
#define MAX_INT INT_MAX
#endif

typedef int64_t int64;
typedef uint64_t uint64;

/****************************************************************************/

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

/****************************************************************************/
#define ERR_ALLOC "Error allocating memory.\n"

//-------------------

#define GEXIT(a) { \
fprintf(stderr, "Error: "); fprintf(stderr, a); \
GError("Exiting from line %i in file %s\n",__LINE__,__FILE__); \
}

// Debug helpers
#ifndef NDEBUG
 #define GASSERT(exp) ((exp)?((void)0):(void)GAssert(#exp,__FILE__,__LINE__))
 #define GVERIFY(condition) \
if (!(condition)) { \
fprintf(stderr, "Assumption \"%s\"\nFailed in file %s: at line:%i\n", \
#condition,__FILE__,__LINE__); \
GEXIT(#condition);}
 #ifdef TRACE
  #define GTRACE(exp)  (GMessage(exp))
 #else
  #define GTRACE(exp)
 #endif
#else
 #define GASSERT(exp)
 #define GTRACE(exp)
 #define GVERIFY(condition)
#endif

#define GERROR(exp) (GError(exp))

// Abolute value
#define GABS(val) (((val)>=0)?(val):-(val))

// Min and Max
#define GMAX(a,b) (((a)>(b))?(a):(b))
#define GMIN(a,b) (((a)>(b))?(b):(a))

// Min of three
#define GMIN3(x,y,z) ((x)<(y)?GMIN(x,z):GMIN(y,z))

// Max of three
#define GMAX3(x,y,z) ((x)>(y)?GMAX(x,z):GMAX(y,z))

// Return minimum and maximum of a, b
#define GMINMAX(lo,hi,a,b) ((a)<(b)?((lo)=(a),(hi)=(b)):((lo)=(b),(hi)=(a)))

// Clamp value x to range [lo..hi]
#define GCLAMP(lo,x,hi) ((x)<(lo)?(lo):((x)>(hi)?(hi):(x)))

typedef void* pointer;
typedef unsigned int uint;

typedef int GCompareProc(const pointer item1, const pointer item2);
typedef long GFStoreProc(const pointer item1, FILE* fstorage); //for serialization
typedef pointer GFLoadProc(FILE* fstorage); //for deserialization

typedef void GFreeProc(pointer item); //usually just delete,
      //but may also support structures with embedded dynamic members

#define GMALLOC(ptr,size)  if (!GMalloc((pointer*)(&ptr),size)) \
                                     GError(ERR_ALLOC)
#define GCALLOC(ptr,size)  if (!GCalloc((pointer*)(&ptr),size)) \
                                     GError(ERR_ALLOC)
#define GREALLOC(ptr,size) if (!GRealloc((pointer*)(&ptr),size)) \
                                     GError(ERR_ALLOC)
#define GFREE(ptr)       GFree((pointer*)(&ptr))

inline char* strMin(char *arg1, char *arg2) {
    return (strcmp(arg1, arg2) < 0)? arg1 : arg2;
}

inline char* strMax(char *arg1, char *arg2) {
    return (strcmp(arg2, arg1) < 0)? arg1 : arg2;
}

inline int iround(double x) {
   return (int)floor(x + 0.5);
}

int Gmkdir(const char *path, bool recursive=true, int perms = (S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH));
void Gmktempdir(char* templ);

/****************************************************************************/

inline int Gintcmp(int a, int b) {
 //return (a>b)? 1 : ((a==b)?0:-1);
  return a-b;
}

int Gstrcmp(const char* a, const char* b, int n=-1);
//same as strcmp but doesn't crash on NULL pointers

int Gstricmp(const char* a, const char* b, int n=-1);
bool GstrEq(const char* a, const char* b);
bool GstriEq(const char* a, const char* b);

//basic swap template function
template<class T> void Gswap(T& lhs, T& rhs) {
 T tmp=lhs; //requires copy operator
 lhs=rhs;
 rhs=tmp;
}


/**************** Memory management ***************************/

bool GMalloc(pointer* ptr, unsigned long size); // Allocate memory
bool GCalloc(pointer* ptr, unsigned long size); // Allocate and initialize memory
bool GRealloc(pointer* ptr,unsigned long size); // Resize memory
void GFree(pointer* ptr); // Free memory, resets ptr to NULL

//int saprintf(char **retp, const char *fmt, ...);

void GError(const char* format,...); // Error routine (aborts program)
void GMessage(const char* format,...);// Log message to stderr
// Assert failed routine:- usually not called directly but through GASSERT
void GAssert(const char* expression, const char* filename, unsigned int lineno);

// ****************** basic string manipulation *************************
char *Gstrdup(const char* str, int xtracap=0); //string duplication with extra capacity added
//duplicate a string by allocating a copy for it (+xtracap heap room) and returning the new pointer
//caller is responsible for deallocating the returned pointer!

char* Gstrdup(const char* sfrom, const char* sto);
//same as GStrdup, but with an early termination (e.g. on delimiter)

char* Gsubstr(const char* str, char* from, char* to=NULL);
//extracts a substring, allocating it, including boundaries (from/to)

char* replaceStr(char* &str, char* newvalue);

//conversion: to Lower/Upper case
// creating a new string:
char* upCase(const char* str);
char* loCase(const char* str);
// changing string in place:
char* strlower(char * str);
char* strupper(char * str);


//strstr but for memory zones: scans a memory region
//for a substring:
void* Gmemscan(void *mem, unsigned int len,
                  void *part, unsigned int partlen);


FILE* Gfopen(const char *path, char *mode=NULL);

// test if a char is in a string:
bool chrInStr(char c, const char* str);

char* rstrchr(char* str, char ch);
/* returns a pointer to the rightmost
  occurence of ch in str - like rindex for platforms missing it*/

char* strchrs(const char* s, const char* chrs);
//strchr but with a set of chars instead of only one

char* rstrfind(const char* str, const char *substr);
// like rindex() but for strings;  right side version of strstr()

char* reverseChars(char* str, int slen=0); //in place reversal of string

char* rstrstr(const char* rstart, const char *lend, const char* substr);
/*the reversed, rightside equivalent of strstr: starts searching
 from right end (rstart), going back to left end (lend) and returns
 a pointer to the last (right) matching character in str */

char* strifind(const char* str,  const char* substr);
// case insensitive version of strstr -- finding a string within another string
// returns NULL if not found

//Determines if a string begins with a given prefix
//(returns false when any of the params is NULL,
// but true when prefix is '' (empty string)!)
bool startsWith(const char* s, const char* prefix);

bool startsiWith(const char* s, const char* prefix); //case insensitive


bool endsWith(const char* s, const char* suffix);
//Note: returns true if suffix is empty string, but false if it's NULL
bool endsiWith(const char* s, const char* suffix); //case insensitive version

//like endsWith but also remove the suffix if found
//returns true if the given suffix was found and removed
bool trimSuffix(char* s, const char* suffix);
//case insensitive version:
bool trimiSuffix(char* s, const char* suffix);

// ELF hash function for strings
int strhash(const char* str);

//alternate hash functions:
int fnv1a_hash(const char* cp);
int djb_hash(const char* cp);

//---- generic base GSeg : genomic segment (interval) --
// coordinates are considered 1-based (so 0 is invalid)
class GSeg {
 public:
  uint start; //start<end always!
  uint end;
  GSeg(uint s=0,uint e=0) {
    if (s>e) { start=e;end=s; }
        else { start=s;end=e; }
  }
  //check for overlap with other segment
  uint len() { return end-start+1; }
  bool overlap(GSeg* d) {
     //return start<d->start ? (d->start<=end) : (start<=d->end);
     return (start<=d->end && end>=d->start);
  }

  bool overlap(GSeg& d) {
     //return start<d.start ? (d.start<=end) : (start<=d.end);
     return (start<=d.end && end>=d.start);
  }

  bool overlap(GSeg& d, int fuzz) {
     //return start<d.start ? (d.start<=end+fuzz) : (start<=d.end+fuzz);
     return (start<=d.end+fuzz && end+fuzz>=d.start);
  }

  bool overlap(uint x) {
	return (start<=x && x<=end);
  }

  bool overlap(uint s, uint e) {
     if (s>e) { Gswap(s,e); }
     return (start<=e && end>=s);
  }

  //return the length of overlap between two segments
  int overlapLen(GSeg* r) {
     if (start<r->start) {
        if (r->start>end) return 0;
        return (r->end>end) ? end-r->start+1 : r->end-r->start+1;
        }
       else { //r->start<=start
        if (start>r->end) return 0;
        return (r->end<end)? r->end-start+1 : end-start+1;
        }
     }
  int overlapLen(uint rstart, uint rend) {
     if (rstart>rend) { Gswap(rstart,rend); }
     if (start<rstart) {
        if (rstart>end) return 0;
        return (rend>end) ? end-rstart+1 : rend-rstart+1;
        }
       else { //rstart<=start
        if (start>rend) return 0;
        return (rend<end)? rend-start+1 : end-start+1;
        }
  }

  bool contains(GSeg* s) {
	  return (start<=s->start && end>=s->end);
  }
  bool contained(GSeg* s) {
	  return (s->start<=start && s->end>=end);
  }

  bool equals(GSeg& d){
      return (start==d.start && end==d.end);
  }
  bool equals(GSeg* d){
      return (start==d->start && end==d->end);
  }

  //fuzzy coordinate matching:
  bool coordMatch(GSeg* s, uint fuzz=0) { //caller must check for s!=NULL
    if (fuzz==0) return (start==s->start && end==s->end);
    uint sd = (start>s->start) ? start-s->start : s->start-start;
    uint ed = (end>s->end) ? end-s->end : s->end-end;
    return (sd<=fuzz && ed<=fuzz);
  }
  void expand(int by) { //expand in both directions
	  start-=by;
	  end+=by;
  }
  void expandInclude(uint rstart, uint rend) { //expand to include given coordinates
	 if (rstart>rend) { Gswap(rstart,rend); }
	 if (rstart<start) start=rstart;
	 if (rend>end) end=rend;
  }
  //comparison operators required for sorting
  bool operator==(GSeg& d){
      return (start==d.start && end==d.end);
      }
  bool operator<(GSeg& d){
     return (start==d.start)?(end<d.end):(start<d.start);
     }
};

//basic dynamic array template for primitive types
//which can only grow (reallocate) as needed

//optimize index test
#define GDynArray_INDEX_ERR "Error: use of index (%d) in GDynArray of size %d!\n"
 #if defined(NDEBUG) || defined(NODEBUG) || defined(_NDEBUG) || defined(NO_DEBUG)
 #define GDynArray_TEST_INDEX(x)
#else
 #define GDynArray_TEST_INDEX(x) \
 if (fCount==0 || x>=fCount) GError(GDynArray_INDEX_ERR, x, fCount)
#endif

#define GDynArray_MAXCOUNT UINT_MAX-1
#define GDynArray_NOIDX UINT_MAX

//basic dynamic array (vector) template for simple/primitive types or structs
//Warning: uses malloc so it will never call the item's default constructor when growing

template<class OBJ> class GDynArray {
 protected:
	bool byptr; //in-place copy (pointer) takeover of existing OBJ[]
    OBJ *fArray;
    uint fCount;
    uint fCapacity; // size of allocated memory
	const static uint dyn_array_defcap = 16; // initial capacity (in elements)
 public:
    GDynArray(int initcap=dyn_array_defcap):byptr(false), fArray(NULL), fCount(0),
	     fCapacity(initcap) { // constructor
    	  GMALLOC(fArray, fCapacity*sizeof(OBJ));
    }

    GDynArray(const GDynArray &a):byptr(false), fArray(NULL),
    		fCount(a.fCount), fCapacity(a.fCapacity) { // copy constructor
        GMALLOC(fArray, sizeof(OBJ)*a.fCapacity);
        memcpy(fArray, a.fArray, sizeof(OBJ)* a.fCapacity);
    }
    GDynArray(OBJ* ptr, uint pcap):byptr(true), fArray(ptr), fCount(0), fCapacity(pcap) {
    	//this will never deallocate the passed pointer
    }

    virtual ~GDynArray() { if (!byptr) { GFREE(fArray); } }

    GDynArray& operator = (const GDynArray &a) { // assignment operator
        if (this == &a) return *this;
    	if (a.fCount == 0) {
    		Clear();
    		return *this;
    	}
    	growTo(a.fCapacity); //set size
        memcpy(fArray, a.fArray, sizeof(OBJ)*a.fCount);
        return *this;
    }

    OBJ& operator[] (uint idx) {// get array item
    	GDynArray_TEST_INDEX(idx);
    	return fArray[idx];
    }

    void Grow() {
    	int delta = (fCapacity>16) ? (fCapacity>>2) : 2;
    	if (GDynArray_MAXCOUNT-delta<=fCapacity)
    		delta=GDynArray_MAXCOUNT-fCapacity;
    	if (delta<=1) GError("Error at GDynArray::Grow(): max capacity reached!\n");
    	growTo(fCapacity + delta);
    }
#define GDynArray_ADD(item) \
    	if (fCount==MAX_UINT-1) GError("Error at GDynArray: cannot add item, maximum count reached!\n"); \
    	if ((++fCount) > fCapacity) Grow(); \
    	fArray[fCount-1] = item;

    uint Add(OBJ* item) { // Add item to the end of array
      //element given by pointer
      if (item==NULL) return GDynArray_NOIDX;
  	  GDynArray_ADD( (*item) );
  	  return (fCount-1);
    }

    uint Add(OBJ item) { // Add OBJ copy to the end of array
	  GDynArray_ADD(item);
	  return (fCount-1);
    }

    uint Push(OBJ item) { //same as Add
    	GDynArray_ADD(item);
    	return (fCount-1);
    }

    OBJ Pop() { //shoddy.. Do NOT call this for an empty array!
    	if (fCount==0) return (OBJ)NULL; //a NULL cast operator is required
    	--fCount;
    	return fArray[fCount];
    }

    uint Count() { return fCount; } // get size of array (elements)
    uint Capacity() { return fCapacity; }
    void growTo(uint newcap) {
    	if (newcap==0) { Clear(); return; }
    	if (newcap <= fCapacity) return; //never shrink! (use Pack() for shrinking)
    	GREALLOC(fArray, newcap*sizeof(OBJ));
    	fCapacity=newcap;
    }

    void append(OBJ* arr, uint count) {
    	//fast adding of a series of objects
    	growTo(fCount+count);
    	memcpy(fArray+fCount, arr, count*sizeof(OBJ));
    	fCount+=count;
    }

    void append(GDynArray<OBJ> arr) {
    	//fast adding of a series of objects
    	growTo(fCount+arr.fCount);
    	memcpy(fArray+fCount, arr.fArray, arr.fCount*sizeof(OBJ));
    	fCount+=arr.fCount;
    }

    void Trim(int tcount=1) {
    	//simply cut (discard) the last tcount items
    	//new Count is now fCount-tcount
    	//does NOT shrink capacity accordingly!
    	if (fCount>=tcount) fCount-=tcount;
    }

    void Pack() { //shrink capacity to fCount+dyn_array_defcap
    	if (fCapacity-fCount<=dyn_array_defcap) return;
    	int newcap=fCount+dyn_array_defcap;
    	GREALLOC(fArray, newcap*sizeof(OBJ));
    	fCapacity=newcap;
    }

    void zPack(OBJ z) { //shrink capacity to fCount+1 and adds a z terminator
    	if (fCapacity-fCount<=1) { fArray[fCount]=z; return; }
    	int newcap=fCount+1;
    	GREALLOC(fArray, newcap*sizeof(OBJ));
    	fCapacity=newcap;
    	fArray[fCount]=z;
    }


    inline void Shrink() { Pack(); }

    void Delete(uint idx) {
	  GDynArray_TEST_INDEX(idx);
	  --fCount;
	  if (idx<fCount)
		  memmove(&fArray[idx], &fArray[idx+1], (fCount-idx)*sizeof(OBJ));
    }

    inline void Remove(uint idx) { Delete(idx); }

    void Clear() { // clear array, shrinking its allocated memory
    	fCount = 0;
    	GREALLOC(fArray, sizeof(OBJ)*dyn_array_defcap);
    	// set initial memory size again
    	fCapacity = dyn_array_defcap;
    }

    void Reset() {// fast clear array WITHOUT deallocating it
    	fCount = 0;
    }

    OBJ* operator()() { return fArray; }

    //use methods below in order to prevent deallocation of fArray pointer on destruct
    //could be handy for adopting stack objects (e.g. cheap dynamic strings)
    void ForgetPtr() { byptr=true;  }
    void DetachPtr() { byptr=true;  }

};


int strsplit(char* str, GDynArray<char*>& fields, const char* delim, int maxfields=MAX_INT);
//splits a string by placing 0 where any of the delim chars are found, setting fields[] to the beginning
//of each field (stopping after maxfields); returns number of fields parsed

int strsplit(char* str, GDynArray<char*>& fields, const char delim, int maxfields=MAX_INT);
//splits a string by placing 0 where the delim char is found, setting fields[] to the beginning
//of each field (stopping after maxfields); returns number of fields parsed

int strsplit(char* str, GDynArray<char*>& fields, int maxfields=MAX_INT); //splits by tab or space
//splits a string by placing 0 where tab or space is found, setting fields[] to the beginning
//of each field (stopping after maxfields); returns number of fields parsed

// ************** simple line reading class for text files
//GLineReader -- text line reading/buffering class
class GLineReader {
   bool closeFile;
   //int len;
   //int allocated;
   GDynArray<char> buf;
   int textlen; //length of actual text, without newline character(s)
   bool isEOF;
   FILE* file;
   off_t filepos; //current position
   bool pushed; //pushed back
   int lcount; //line counter (read lines)
 public:
   char* chars() { return buf(); }
   char* line() { return buf(); }
   int readcount() { return lcount; } //number of lines read
   void setFile(FILE* stream) { file=stream; }
   int blength() { return buf.Count(); } //binary/buffer length, including newline character(s)
   int charcount() { return buf.Count(); } //line length, including newline character(s)
   int tlength() { return textlen; } //text length excluding newline character(s)
   int linelen() { return textlen; } //line length, excluding newline character(s)
   //int size() { return buf.Count(); } //same as size();
   bool isEof() {return isEOF; }
   bool eof() { return isEOF; }
   off_t getfpos() { return filepos; }
   off_t getFpos() { return filepos; }
   char* nextLine() { return getLine(); }
   char* getLine() { if (pushed) { pushed=false; return buf(); }
                            else return getLine(file);  }
   char* getLine(FILE* stream) {
                 if (pushed) { pushed=false; return buf(); }
                          else return getLine(stream, filepos); }
   char* getLine(FILE* stream, off_t& f_pos); //read a line from a stream and update
                           // the given file position
   void pushBack() { if (lcount>0) pushed=true; } // "undo" the last getLine request
            // so the next call will in fact return the same line
   GLineReader(const char* fname):closeFile(false),buf(1024), textlen(0),
		   isEOF(false),file(NULL),filepos(0), pushed(false), lcount(0) {
      FILE* f=fopen(fname, "rb");
      if (f==NULL) GError("Error opening file '%s'!\n",fname);
      closeFile=true;
      file=f;
      }
   GLineReader(FILE* stream=NULL, off_t fpos=0):closeFile(false),buf(1024),
		   textlen(0), isEOF(false),file(stream),
		   filepos(fpos), pushed(false), lcount(0) {
     }
   ~GLineReader() {
     if (closeFile) fclose(file);
     }
};


/* extended fgets() -  to read one full line from a file and
  update the file position correctly !
  buf will be reallocated as necessary, to fit the whole line
  */
char* fgetline(char* & buf, int& buflen, FILE* stream, off_t* f_pos=NULL, int* linelen=NULL);


//print int/values nicely formatted in 3-digit groups
char* commaprintnum(uint64 n);

/*********************** File management functions *********************/

// removes the last part (file or directory name) of a full path
// WARNING: this is a destructive operation for the given string!
void delFileName(char* filepath);

// returns a pointer to the last file or directory name in a full path
const char* getFileName(const char* filepath);
// returns a pointer to the file "extension" part in a filename
const char* getFileExt(const char* filepath);


int fileExists(const char* fname);
//returns 0 if path doesn't exist
//        1 if it's a directory
//        2 if it's a regular file
//        3 something else (but entry exists)

int64 fileSize(const char* fpath);

//write a formatted fasta record, fasta formatted
void writeFasta(FILE *fw, const char* seqid, const char* descr,
        const char* seq, int linelen=60, int seqlen=0);

//parses the next number found in a string at the current position
//until a non-digit (and not a '.', 'e','E','-','+') is encountered;
//updates the char* pointer to be after the last digit parsed
bool parseNumber(char* &p, double& v);
bool parseDouble(char* &p, double& v); //just an alias for parseNumber
bool parseFloat(char* &p, float& v);

bool strToInt(char* p, int& i);
bool strToUInt(char* p, uint& i);
bool parseInt(char* &p, int& i); //advance pointer p after the number
bool parseUInt(char* &p, uint& i); //advance pointer p after the number
bool parseHex(char* &p,  uint& i);

#endif /* G_BASE_DEFINED */
