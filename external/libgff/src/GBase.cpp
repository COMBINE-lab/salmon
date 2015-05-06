#include "GBase.h"
#include <stdarg.h>
#include <ctype.h>
#include <sys/stat.h>

#ifndef S_ISDIR
#define S_ISDIR(mode)  (((mode) & S_IFMT) == S_IFDIR)
#endif

#ifndef S_ISREG
#define S_ISREG(mode)  (((mode) & S_IFMT) == S_IFREG)
#endif

/*
#ifdef _DEFINE_WIN32_FSEEKO
 int fseeko(FILE *stream, off_t offset, int whence) {
   
   }
#endif

#ifdef _DEFINE_WIN32_FTELLO
 off_t ftello(FILE *stream) {
  
  }
#endif
*/

/*
int saprintf(char **retp, const char *fmt, ...) {
  va_list argp;
  int len;
  char *buf;

  va_start(argp, fmt);
  len = vsnprintf(NULL, 0, fmt, argp);
  va_end(argp);
  GMALLOC(buf, (len + 1));
  if(buf == NULL)
    {
    *retp = NULL;
    return -1;
    }

  va_start(argp, fmt);
  vsnprintf(buf, len+1, fmt, argp);
  va_end(argp);

  *retp = buf;
  return len;
}
*/

//************************* Debug helpers **************************
// Assert failed routine
void GAssert(const char* expression, const char* filename, unsigned int lineno){
  char msg[4096];
  sprintf(msg,"%s(%d): ASSERT(%s) failed.\n",filename,lineno,expression);
  fprintf(stderr,"%s",msg);
  //abort();
  }
// Error routine (prints error message and exits!)
void GError(const char* format,...){
  #ifdef __WIN32__
    char msg[4096];
    va_list arguments;
    va_start(arguments,format);
    vsprintf(msg,format,arguments);
    va_end(arguments);
    OutputDebugString(msg);
    fprintf(stderr,"%s",msg); // if a console is available
    MessageBox(NULL,msg,NULL,MB_OK|MB_ICONEXCLAMATION|MB_APPLMODAL);
  #else
    va_list arguments;
    va_start(arguments,format);
    vfprintf(stderr,format,arguments);
    va_end(arguments);
    #ifdef DEBUG
     // modify here if you want a core dump
     abort();
    #endif
  #endif
    exit(1);
  }
  
// Warning routine (just print message without exiting)
void GMessage(const char* format,...){
  char msg[4096];
  va_list arguments;
  va_start(arguments,format);
  vsprintf(msg,format,arguments);
  va_end(arguments);
  #ifdef __WIN32__
    OutputDebugString(msg);
  #endif
  fprintf(stderr,"%s",msg);fflush(stderr);
  }

/*************** Memory management routines *****************/
// Allocate memory
bool GMalloc(pointer* ptr,unsigned long size){
  //GASSERT(ptr);
  if (size!=0) *ptr=malloc(size);
  return *ptr!=NULL;
  }

// Allocate cleaned memory (0 filled)
bool GCalloc(pointer* ptr,unsigned long size){
  GASSERT(ptr);
  *ptr=calloc(size,1);
  return *ptr!=NULL;
  }

// Resize memory
bool GRealloc(pointer* ptr,unsigned long size){
  //GASSERT(ptr);
  if (size==0) {
    GFree(ptr);
    return true;
    }
  if (*ptr==NULL) {//simple malloc
   void *p=malloc(size);
   if (p != NULL) {
     *ptr=p;
     return true;
     }
    else return false;
   }//malloc
  else {//realloc
   void *p=realloc(*ptr,size);
   if (p) {
       *ptr=p;
       return true;
       }
   return false;
   }
 }
// Free memory, resets ptr to NULL afterward
void GFree(pointer* ptr){
  GASSERT(ptr);
  if (*ptr) free(*ptr);
  *ptr=NULL;
  }

char* Gstrdup(const char* str) {
  if (str==NULL) return NULL;
  char *copy=NULL;
  GMALLOC(copy, strlen(str)+1);
  strcpy(copy,str);
  return copy;
  }

char* newEmptyStr() {
  char* zs=NULL;
  GMALLOC(zs,1);
  zs[0]=0;
  return zs;
}

char* Gstrdup(const char* sfrom, const char* sto) {
  if (sfrom==NULL || sto==NULL) return NULL;
  char *copy=NULL;
  if (sfrom[0]==0) return newEmptyStr();
  GMALLOC(copy, sto-sfrom+2);
  strncpy(copy, sfrom, sto-sfrom+1);
  copy[sto-sfrom+1]=0;
  return copy;
  }

int Gstrcmp(const char* a, const char* b, int n) {
 if (a==NULL || b==NULL) {
   return a==NULL ? -1 : 1;
   }
 else {
   if (n<0) return strcmp(a,b);
       else return strncmp(a,b,n);
 }

}

int Gstricmp(const char* a, const char* b, int n) {
 if (a==NULL || b==NULL) return a==NULL ? -1 : 1;
 register int ua, ub;
 if (n<0) {
   while ((*a!=0) && (*b!=0)) {
    ua=tolower((unsigned char)*a);
    ub=tolower((unsigned char)*b);
    a++;b++;
    if (ua!=ub) return ua < ub ? -1 : 1;
    }
    return (*a == 0) ? ( (*b == 0) ? 0 : -1 ) : 1 ;
  }
 else {
   while (n && (*a!=0) && (*b!=0)) {
    ua=tolower((unsigned char)*a);
    ub=tolower((unsigned char)*b);
    a++;b++;n--;
    if (ua!=ub) return ua < ub ? -1 : 1;
    }
    //return (*a == 0) ? ( (*b == 0) ? 0 : -1 ) : 1 ;
   if (n==0) return 0;
   else { return (*a == 0) ? ( (*b == 0) ? 0 : -1 ) : 1 ; }
  }
}

int strsplit(char* str, char** fields, int maxfields, const char* delim) {
 //splits by placing 0 where delim chars are found, setting fields[] to the beginning
 //of each field (stopping after maxfields); returns number of fields parsed
 int tidx=0;
 bool afterdelim=true;
 int i=0;
 while (str[i]!=0 && tidx<maxfields) {
    if (afterdelim) {
        fields[tidx]=str+i;
        tidx++;
        }
    afterdelim=false;
    if (chrInStr(str[i],(char*)delim)) {
        str[i]=0;
        i++;
        while (str[i]!=0 && chrInStr(str[i], (char*)delim)) i++;
        afterdelim=true;
        continue;
        }
    i++;
    }
 return tidx;
}

int strsplit(char* str, char** fields, int maxfields, const char delim) {
  //splits by placing 0 where delim is found, setting fields[] to the beginning
  //of each field (stopping after maxfields); returns number of fields parsed
  int tidx=0;
  bool afterdelim=true;
  int i=0;
  while (str[i]!=0 && tidx<maxfields) {
     if (afterdelim) {
         fields[tidx]=str+i;
         tidx++;
         }
     afterdelim=false;
     if (str[i]==delim) {
         str[i]=0;
         i++;
         while (str[i]!=0 && str[i]==delim) i++;
         afterdelim=true;
         continue;
         }
     i++;
     }
  return tidx;
}

int strsplit(char* str, char** fields, int maxfields) {
  //splits by placing 0 where delim is found, setting fields[] to the beginning
  //of each field (stopping after maxfields); returns number of fields parsed
  int tidx=0;
  bool afterdelim=true;
  int i=0;
  while (str[i]!=0 && tidx<maxfields) {
     if (afterdelim) {
         fields[tidx]=str+i;
         tidx++;
         }
     afterdelim=false;
     if (str[i]==' ' || str[i]=='\t') {
         str[i]=0;
         i++;
         while (str[i]!=0 && (str[i]=='\t' || str[i]==' ')) i++;
         afterdelim=true;
         continue;
         }
     i++;
     }
  return tidx;
}


char* Gsubstr(const char* str, char* from, char* to) {
 //extract (and allocate) a substring, including boundaries (from/to)
 if (str==NULL || from==NULL) return NULL;
 if (from[0]==0 || str[0]==0) return newEmptyStr();
 if (from<str) return NULL;
 if (to==NULL) {
    to=from;
    while (to[1]) to++;
    }
 if (to<from) return newEmptyStr();
 int newlen=to-from+1;
 char* subs;
 GMALLOC(subs, newlen);
 memcpy(subs, str, newlen-1);
 subs[newlen]='\0';
 return subs;
 }

char* replaceStr(char* &str, char* newvalue) {
 if (str!=NULL) GFREE(str);
 if (newvalue==NULL) { return NULL; }
 GMALLOC(str, strlen(newvalue)+1);
 strcpy(str,newvalue);
 return str;
 }

void* Gmemscan(void *mem, unsigned int len,
                   void *part, unsigned int partlen) {
char* p;
unsigned int restlen=len-partlen+1;
void* oldp=mem;
while ( (p=(char*)memchr(oldp, ((char*)part)[0], restlen))!=NULL) {
  //located first char, try to match the rest:
  p++;
  if (memcmp(p, &((char*)part)[1], partlen-1)==0) return p-1;
  //no string match, prepare next iteration
  restlen-=(p-(char*)oldp);
  oldp=p;
  }//while
return NULL;
}

//rindex function is missing on some platforms ?
char* rstrchr(char* str, char ch) {  /* returns a pointer to the rightmost
  occurence of ch in str  */
 char *p;
 if (str==NULL) return NULL;
 p=str+strlen(str)-1;
 while (p>=str) {
    if (*p==ch) return p;
    p--;
    }
 return NULL;
 }


/* DOS/UNIX safer fgets : reads a text line from a (binary) file and
  update the file position accordingly and the buffer capacity accordingly.
  The given buf is resized to read the entire line in memory
    -- even when it's abnormally long
  */
char* fgetline(char* & buf, int& buf_cap, FILE *stream, off_t* f_pos, int* linelen) {
  //reads a char at a time until \n and/or \r are encountered
  int i=0;
  int c=0;
  off_t fpos=(f_pos!=NULL) ? *f_pos : 0;
  while ((c=getc(stream))!=EOF) {
    if (i>=buf_cap-1) {
       buf_cap+=1024;
       GREALLOC(buf, buf_cap);
       }
    if (c=='\n' || c=='\r') {
       if (c=='\r') {
         if ((c=getc(stream))!='\n') ungetc(c,stream);
                                else fpos++;
         }
       fpos++;
       break;
       }
    fpos++;
    buf[i]=(char)c;
    i++;
    } //while i<buf_cap-1
  if (linelen!=NULL) *linelen=i;
  if (f_pos!=NULL) *f_pos=fpos;
  if (c==EOF && i==0) return NULL;
  buf[i]='\0';
  return buf;
  }

char* GLineReader::getLine(FILE* stream, off_t& f_pos) {
   if (pushed) { pushed=false; return buf; }
   //reads a char at a time until \n and/or \r are encountered
   len=0;
   int c=0;
   while ((c=getc(stream))!=EOF) {
     if (len>=allocated-1) {
        allocated+=1024;
        GREALLOC(buf, allocated);
     }
     if (c=='\n' || c=='\r') {
       buf[len]='\0';
       if (c=='\r') { //DOS file -- special case
         if ((c=getc(stream))!='\n') ungetc(c,stream);
                                else f_pos++;
         }
       f_pos++;
       lcount++;
       return buf;
       }
     f_pos++;
     buf[len]=(char)c;
     len++;
     } //while i<buf_cap-1
   if (c==EOF) {
     isEOF=true;
     if (len==0) return NULL;
     }
   buf[len]='\0';
   lcount++;
   return buf;
}


//strchr but with a set of chars instead of only one
char* strchrs(const char* s, const char* chrs) {
  if (s==NULL || chrs==NULL || *chrs=='\0' || *s=='\0')
         return NULL;
  unsigned int l=strlen(s);
  unsigned int r=strcspn(s, chrs);
  if (r==l) return NULL;
  return ((char*)s+r);
}

char* upCase(const char* str) {
 if (str==NULL) return NULL;
 int len=strlen(str);
 char* upstr;
 GMALLOC(upstr, len+1);
 upstr[len]='\0';
 for (int i=0;i<len;i++) upstr[i]=toupper(str[i]);
 return upstr;
 }

char* loCase(const char* str) {
 if (str==NULL) return NULL;
 int len=strlen(str);
 char* lostr;
 GMALLOC(lostr, len+1);
 lostr[len]='\0';
 for (int i=0;i<len;i++) lostr[i]=tolower(str[i]);
 return lostr;
 }

char* strlower(char * str) {//changes string in place
  if (str==NULL) return NULL;
  int i=0;
  while (str[i]!=0) { str[i]=tolower(str[i]); i++; }
  return str;
}

char* strupper(char * str) {//changes string in place
  if (str==NULL) return NULL;
  int i=0;
  while (str[i]!=0) { str[i]=toupper(str[i]); i++; }
  return str;
}



//test if a char is in a given string (set)
bool chrInStr(char c, const char* str) {
 if (str==NULL || *str=='\0') return false;
 for (const char* p=str; (*p)!='\0'; p++) {
   if ((*p)==c) return true;
   }
 return false;
 }



char* rstrfind(const char* str, const char* substr) {
/* like rindex() for a string */
 int l,i;
 if (str==NULL || *str=='\0') return NULL;
 if (substr==NULL || *substr=='\0') return NULL;
 l=strlen(substr);
 char* p=(char*)str+strlen(str)-l;
   //rightmost position that could match

 while (p>=str) {
    for (i=0; i<l && *(p+i) == *(substr+i); i++) ;
    if (i==l) return p; //found!
    p--;
    }
 return NULL;
}


char* strifind(const char* str,  const char* substr) {
 // the case insensitive version of strstr -- finding a string within a strin
  int l,i;
  if (str==NULL || *str==0) return NULL;
  if (substr==NULL || *substr==0) return NULL;
  l=strlen(substr);
  char* smax=(char*)str+strlen(str)-l;
  //rightmost position that could match
  char* p=(char*)str;
  while (p<=smax) {
     for (i=0; i<l && tolower(*(p+i))==tolower(*(substr+i)); i++) ;
     if (i==l) return p; //found!
     p++;
     }
  return NULL;
}



// tests if string s has the given prefix
bool startsWith(const char* s, const char* prefix) {
 if (prefix==NULL || s==NULL) return false;
 int i=0;
 while (prefix[i]!='\0' && prefix[i]==s[i]) i++;
 return (prefix[i]=='\0');
 }

// tests if string s ends with given suffix
bool endsWith(const char* s, const char* suffix) {
 if (suffix==NULL || s==NULL) return false;
 if (suffix[0]==0) return true; //special case: empty suffix
 int j=strlen(suffix)-1;
 int i=strlen(s)-1;
 if (i<j) return false;
 while (j>=0 && s[i]==suffix[j]) { i--; j--; }
 return (j==-1);
 }


char* reverseChars(char* str, int slen) {
  if (slen==0) slen=strlen(str);
  int l=0;
  int r=slen-1;
  char c;
  while (l<r) {
     c=str[l];str[l]=str[r];
     str[r]=c;
     l++;r--;
     }
  return str;
}


char* rstrstr(const char* rstart, const char *lend, const char* substr) {  /*like strstr, but starts searching
 from right end, going up to lend and returns a pointer to the last (right)
 matching character in str */
 char *p;
 int l,i;
 l=strlen(substr);
 p=(char*)rstart-l+1;
 while (p>=lend) {
    for (i=0;i<l;i++) if (*(p+i) != *(substr+i)) break;
    if (i==l) return p+l-1;
    p--;
    }
 return NULL;
 }


//hash function used for strings in GHash
int strhash(const char* str){
  register int h=0;
  register int g;
  while (*str) {
    h=(h<<4)+*str++;
    g=h&0xF0000000;
    if(g) h^=g>>24;
    h&=0x0fffffff;
    }
  GASSERT(h<=0x0fffffff);
  return h;
  }

// removes the last part (file or directory name) of a full path
// this is a destructive operation for the given string!!!
// the trailing '/' is guaranteed to be there
void delFileName(char* filepath) {
 char *p, *sep;
 if (filepath==NULL) return;
 for (p=filepath, sep=filepath;*p!='\0';p++)
     if (*p=='/' || *p=='\\') sep=p+1;
 *sep='\0'; // truncate filepath
}

// returns a pointer to the last file or directory name in a full path
const char* getFileName(const char* filepath) {
 const char *p, *sep;
 if (filepath==NULL) return NULL;
 for (p=filepath, sep=filepath;*p!='\0';p++)
     if (*p=='/' || *p=='\\') sep=p+1;
 return sep;
}

// returns a pointer to the file "extension" part in a filename
const char* getFileExt(const char* filepath) {
 const char *p, *dp, *sep;
 if (filepath==NULL) return NULL;
 for (p=filepath, dp=filepath, sep=filepath;*p!='\0';p++) {
     if (*p=='.') dp=p+1;
       else if (*p=='/' || *p=='\\') 
                  sep=p+1;
     }
 return (dp>sep) ? dp : NULL ;
}

int fileExists(const char* fname) {
  struct stat stFileInfo;
  int r=0;
  // Attempt to get the file attributes
  int fs = stat(fname,&stFileInfo);
  if (fs == 0) {
      r=3;
      // We were able to get the file attributes
      // so the file obviously exists.
      if (S_ISREG (stFileInfo.st_mode)) {
         r=2;
         }
      if (S_ISDIR (stFileInfo.st_mode)) {
          r=1;
          }
      }
  return r;
}

/*bool fileExists(const char* filepath) {
  if (filepath==NULL) return false;
  FILE* ft=fopen(filepath, "rb");
  if (ft==NULL) return false;
  fclose(ft);
  return true;
}
*/
int64 fileSize(const char* fpath) {
  struct stat results;
  if (stat(fpath, &results) == 0)
      // The size of the file in bytes is in
      return (int64)results.st_size;
  else
      // An error occurred
    //GMessage("Error at stat(%s)!\n", fpath);
    return 0;
}

bool parseNumber(char* &p, double& v) {
 //skip any spaces..
 while (*p==' ' || *p=='\t') p++;
 char* start=p;
 /*if (*p=='-') p++;
       else if (*p=='+') { p++;start++; }*/

 /* while ((*p>='1' && *p<='9') || *p=='0' ||
          *p=='.' || *p=='-' || tolower(*p)=='e') p++; */
 int numlen=strspn(start, "0123456789eE.-+");
 p=start+numlen;
 //now p is on a non-digit;
 if (*start=='-' && p==start+1) return false;
 char saved=*p;
 *p='\0';
 char* endptr=p;
 v=strtod(start,&endptr);
 *p=saved;
 if (endptr!=p) return false;
 return true;
}


bool parseDouble(char* &p, double& v) {
 return parseNumber(p,v);
}

bool parseInt(char* &p, int& i) {
 while (*p==' ' || *p=='\t') p++;
 char* start=p;
 if (*p=='-') p++;
       else if (*p=='+') { p++;start++; }
 while ((*p>='1' && *p<='9') || *p=='0') p++;
 //now p is on a non-digit;
 if (*start=='-' && p==start+1) return false;
 char saved=*p;
 *p='\0';
 char* endptr=p;
 long l=strtol(start,&endptr,10);
 i=(int)l;
 *p=saved;
 if (endptr!=p || i!=l) return false;
 return true;
}

bool parseUInt(char* &p, uint& i) {
 while (*p==' ' || *p=='\t') p++;
 char* start=p;
 if (*p=='-') return false;
       else if (*p=='+') { p++;start++; }
 while ((*p>='1' && *p<='9') || *p=='0') p++;
 //now p is on a non-digit;
 if (*start=='-' && p==start+1) return false;
 char saved=*p;
 *p='\0';
 char* endptr=p;
 unsigned long l=strtoul(start,&endptr,10);
 i=(uint) l;
 *p=saved;
 if (endptr!=p || i!=l) return false;
 return true;
}

bool parseHex(char* &p, uint& i) {
 //skip initial spaces/prefix
 while (*p==' ' || *p=='\t' || *p=='0' || *p=='x') p++;
 char* start=p;
 if (*p=='-') return false;
       else if (*p=='+') { p++;start++; }
 while (isxdigit(*p)) p++;
 //now p is on a non-hexdigit;
 if (p==start+1) return false;
 char saved=*p;
 *p='\0';
 char* endptr=p;
 unsigned long l=strtoul(start,&endptr,16);
 i=(uint) l;
 *p=saved;
 if (endptr!=p || i!=l) return false;
 return true;
}

//write a formatted fasta record, fasta formatted
void writeFasta(FILE *fw, const char* seqid, const char* descr,
        const char* seq, int linelen, int seqlen) {
  fflush(fw);
  // write header line only if given!
  if (seqid!=NULL) {
    if (descr==NULL || descr[0]==0)
             fprintf(fw,">%s\n",seqid);
        else fprintf(fw,">%s %s\n",seqid, descr);
    }
  fflush(fw);
  if (seq==NULL || *seq==0) return; //nothing to print
  if (linelen==0) { //unlimited line length: write the whole sequence on a line
     if (seqlen>0)
             fwrite((const void*)seq, 1, seqlen,fw);
        else fprintf(fw,"%s",seq);
     fprintf(fw,"\n");
     fflush(fw);
     return;
     }
  int ilen=0;
  if (seqlen>0) { //seq length given, so we know when to stop
    for (int i=0; i < seqlen; i++, ilen++) {
            if (ilen == linelen) {
                 fputc('\n', fw);
                 ilen = 0;
                 }
            fputc(seq[i], fw);
            }
    fputc('\n', fw);
    }
  else { //seq length not given, stop when 0 encountered
    for (int i=0; seq[i]!=0; i++, ilen++) {
            if (ilen == linelen) {
                 fputc('\n', fw);
                 ilen = 0;
                 }
            fputc(seq[i], fw);
            } //for
    fputc('\n', fw);
    }
  fflush(fw);
 }

char* commaprint(uint64 n) {
  int comma = '\0';
  char retbuf[48];
  char *p = &retbuf[sizeof(retbuf)-1];
  int i = 0;
  if(comma == '\0') {
    /* struct lconv *lcp = localeconv();
    if(lcp != NULL) {
      if(lcp->thousands_sep != NULL &&
        *lcp->thousands_sep != '\0')
        comma = *lcp->thousands_sep;
      else  */
                          comma = ',';
     // }
    }
  *p = '\0';
  do {
    if(i%3 == 0 && i != 0)
      *--p = comma;
    *--p = '0' + n % 10;
    n /= 10;
    i++;
  } while(n != 0);
  return p;
}
