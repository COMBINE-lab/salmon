#ifndef G_BASE_DEFINED
#define G_BASE_DEFINED
#ifndef _POSIX_SOURCE
//mostly for MinGW
#define _POSIX_SOURCE
#endif
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdint.h>

#if defined __WIN32__ || defined WIN32 || defined _WIN32 || defined _WIN32_
  #ifndef __WIN32__
    #define __WIN32__
  #endif
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
			/*
			#define _DEFINE_WIN32_FSEEKO
			int fseeko(FILE *stream, off_t offset, int whence);
			*/
			#define fseeko fseek
		#endif
  #endif
 #ifndef ftello
  #ifdef _ftelli64
    #define ftello(stream) _ftelli64(stream)
  #else
    /*
    #define _DEFINE_WIN32_FTELLO
    off_t ftello(FILE *stream);
    */
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

// Debug helpers
#ifndef NDEBUG
 #define GASSERT(exp) ((exp)?((void)0):(void)GAssert(#exp,__FILE__,__LINE__))
 #ifdef TRACE
  #define GTRACE(exp)  (GMessage exp)
 #else
  #define GTRACE(exp)  ((void)0)
 #endif
#else
 #define GASSERT(exp) ((void)0)
 #define GTRACE(exp)  ((void)0)
#endif

#define GERROR(exp) (GError exp)
/**********************************  Macros  ***********************************/
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

/****************************************************************************/

inline int Gintcmp(int a, int b) {
 //return (a>b)? 1 : ((a==b)?0:-1);
  return a-b;
}

int Gstrcmp(const char* a, const char* b, int n=-1);
//same as strcmp but doesn't crash on NULL pointers

int Gstricmp(const char* a, const char* b, int n=-1);

//basic swap template function
template<class T> void Gswap(T& lhs, T& rhs) {
 //register T tmp=lhs;
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

// ****************** string manipulation *************************
char *Gstrdup(const char* str);
//duplicate a string by allocating a copy for it and returning it
char* Gstrdup(const char* sfrom, const char* sto);
//same as GStrdup, but with an early termination (e.g. on delimiter)

char* Gsubstr(const char* str, char* from, char* to=NULL);
//extracts a substring, allocating it, including boundaries (from/to)

int strsplit(char* str, char** fields, int maxfields, const char* delim);
int strsplit(char* str, char** fields, int maxfields, const char delim);
int strsplit(char* str, char** fields, int maxfields); //splits by tab or space

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
// the case insensitive version of strstr -- finding a string within a strin


//Determines if a string begins with a given prefix
//(returns false when any of the params is NULL,
// but true when prefix is '' (empty string)!)
bool startsWith(const char* s, const char* prefix);

bool endsWith(const char* s, const char* suffix);
//Note: returns true if suffix is empty string, but false if it's NULL


// ELF hash function for strings
int strhash(const char* str);



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

  bool overlap(uint s, uint e) {
     if (s>e) { Gswap(s,e); }
     //return start<s ? (s<=end) : (start<=e);
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

  //fuzzy coordinate matching:
  bool coordMatch(GSeg* s, uint fuzz=0) {
    if (fuzz==0) return (start==s->start && end==s->end);
    uint sd = (start>s->start) ? start-s->start : s->start-start;
    uint ed = (end>s->end) ? end-s->end : s->end-end;
    return (sd<=fuzz && ed<=fuzz);
    }
  //comparison operators required for sorting
  bool operator==(GSeg& d){
      return (start==d.start && end==d.end);
      }
  bool operator<(GSeg& d){
     return (start==d.start)?(end<d.end):(start<d.start);
     }
};



//--------------------------------------------------------
// ************** simple line reading class for text files

//GLineReader -- text line reading/buffering class
class GLineReader {
   bool closeFile;
   int len;
   int allocated;
   char* buf;
   bool isEOF;
   FILE* file;
   off_t filepos; //current position
   bool pushed; //pushed back
   int lcount; //line counter (read lines)
 public:
   char* chars() { return buf; }
   char* line() { return buf; }
   int readcount() { return lcount; } //number of lines read
   void setFile(FILE* stream) { file=stream; }
   int length() { return len; }
   int size() { return len; } //same as size();
   bool isEof() {return isEOF; }
   bool eof() { return isEOF; }
   off_t getfpos() { return filepos; }
   off_t getFpos() { return filepos; }
   char* nextLine() { return getLine(); }
   char* getLine() { if (pushed) { pushed=false; return buf; }
                            else return getLine(file);  }
   char* getLine(FILE* stream) {
                 if (pushed) { pushed=false; return buf; }
                          else return getLine(stream, filepos); }
   char* getLine(FILE* stream, off_t& f_pos); //read a line from a stream and update
                           // the given file position
   void pushBack() { if (lcount>0) pushed=true; } // "undo" the last getLine request
            // so the next call will in fact return the same line
   GLineReader(const char* fname) {
      FILE* f=fopen(fname, "rb");
      if (f==NULL) GError("Error opening file '%s'!\n",fname);
      closeFile=true;
      init(f);
      }
   GLineReader(FILE* stream=NULL, off_t fpos=0) {
     closeFile=false;
     init(stream,fpos);
     }
   void init(FILE* stream, off_t fpos=0) {
     len=0;
     isEOF=false;
     allocated=1024;
     GMALLOC(buf,allocated);
     lcount=0;
     buf[0]=0;
     file=stream;
     filepos=fpos;
     pushed=false;
     }
   ~GLineReader() {
     GFREE(buf);
     if (closeFile) fclose(file);
     }
};


/* extended fgets() -  to read one full line from a file and
  update the file position correctly !
  buf will be reallocated as necessary, to fit the whole line
  */
char* fgetline(char* & buf, int& buflen, FILE* stream, off_t* f_pos=NULL, int* linelen=NULL);


//print int/values nicely formatted in 3-digit groups
char* commaprint(uint64 n);

/*********************** File management functions *********************/

// removes the last part (file or directory name) of a full path
// WARNING: this is a destructive operation for the given string!
void delFileName(char* filepath);

// returns a pointer to the last file or directory name in a full path
const char* getFileName(const char* filepath);
// returns a pointer to the file "extension" part in a filename
const char* getFileExt(const char* filepath);


int fileExists(const char* fname);
//returns 0 if file entry doesn't exist
//        1 if it's a directory
//        2 if it's a regular file
//        3 otherwise (?)

int64 fileSize(const char* fpath);

//write a formatted fasta record, fasta formatted
void writeFasta(FILE *fw, const char* seqid, const char* descr,
        const char* seq, int linelen=60, int seqlen=0);

//parses the next number found in a string at the current position
//until a non-digit (and not a '.', 'e','E','-','+') is encountered;
//updates the char* pointer to be after the last digit parsed
bool parseNumber(char* &p, double& v);
bool parseDouble(char* &p, double& v); //just an alias for parseNumber

bool parseInt(char* &p, int& i);
bool parseUInt(char* &p, uint& i);
bool parseHex(char* &p,  uint& i);

#endif /* G_BASE_DEFINED */
