//---------------------------------------------------------------------------
#ifndef GSTR_H
#define GSTR_H
//---------------------------------------------------------------------------
#include "GBase.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

// This class uses reference counting and copy-on-write semantics

// All indexes are zero-based.  For all functions that accept an index, a
// negative index specifies an index from the right of the string.  Also,
// for all functions that accept a length, a length of -1 specifies the rest
// of the string.
enum enTokenizeMode {
 tkFullString,
 tkCharSet
 };

class GStr {
        friend GStr operator+(const char* s1, const GStr& s2);
        friend bool operator==(const char* s1, const GStr& s2);
        friend bool operator<(const char* s1, const GStr& s2);
        friend bool operator<=(const char* s1, const GStr& s2);
        friend bool operator>(const char* s1, const GStr& s2);
        friend bool operator>=(const char* s1, const GStr& s2);
        friend bool operator!=(const char* s1, const GStr& s2);
        friend void Gswap(GStr& s1, GStr& s2);
    public:
        GStr();
        GStr(const GStr& s);
        GStr(const char* s, uint addcap=4);
        GStr(const int i);
        GStr(const double f);
        GStr(char c, int n = 1);
        ~GStr();
        operator const char* () const { return my_data->chars;} //inline here
        char& operator[](int index);
        char operator[](int index) const;
        GStr& operator=(const GStr& s);
        GStr& operator=(const char* s);
        GStr& operator=(const int i);
        GStr& operator=(const double f);
        GStr operator+(const GStr& s) const;
        GStr operator+(const char* s) const;
        GStr operator+(const char c) const;
        GStr operator+(const int i) const;
        GStr operator+(const double f) const;
        bool operator==(const GStr& s) const;
        bool operator==(const char* s) const;
        bool operator<(const GStr& s) const;
        bool operator<(const char* s) const;
        bool operator<=(const GStr& s) const;
        bool operator<=(const char* s) const;
        bool operator>(const GStr& s) const;
        bool operator>(const char* s) const;
        bool operator>=(const GStr& s) const;
        bool operator>=(const char* s) const;
        bool operator!=(const GStr& s) const;
        bool operator!=(const char* s) const;
        GStr& operator+=(const GStr& s) { return append(s.chars()); }
        GStr& operator+=(const char* s) { return append(s); }
        GStr& operator+=(char c) { return append(c); }
        GStr& operator+=(int i) { return append(i); }
        GStr& operator+=(uint i) { return append(i); }
        GStr& operator+=(long l) { return append(l); }
        GStr& operator+=(unsigned long l) { return append(l); }
        GStr& operator+=(double f);
      //interface:
        int length() const;
        bool is_empty() const;
        bool is_space() const;
        GStr substr(int index = 0, int len = -1) const;
        GStr to(char c); //return the first part up to first occurence of c
                           //or whole string if c not found
        GStr from(char c); //same as to, but starting from the right side
        GStr copy() const;
        GStr& format(const char *fmt,...);
        GStr& reverse();
        GStr& appendfmt(const char *fmt,...);
        GStr& cut(int index = 0, int len = -1); //delete a specified length
        GStr& remove(int from, int to) {
            return cut(from, to-from+1);
        }

        //paste a string at the specified position
        GStr& paste(const GStr& s, int index = 0, int len=-1);
        GStr& paste(const char* s, int index = 0, int len = -1);
        GStr& replace(const char* from, const char* to=NULL);
        GStr& insert(const GStr& s, int index = 0);
        GStr& insert(const char* s, int index = 0);
        GStr& append(const char* s);
        GStr& appendQuoted(const char* s, char q='"', bool onlyIfSpaced=false);
        GStr& appendmem(const char* m, int len);
        GStr& append(const char* m, int len); //same as appendmem but stops at '\0'

        GStr& append(const GStr& s);
        GStr& append(char c);
        GStr& append(int i);
        GStr& append(long l);
        GStr& append(double f);
        GStr& append(uint i);
        GStr& append(unsigned long l);

        GStr& upper();
        GStr& lower();
        GStr& clear(int init_cap=0);//make empty, but can specify initial capacity
        //character translation or removal:
        GStr& tr(const char* from, const char* to=NULL);
        //number of occurences of a char in the string:
        int count(char c);
        void startTokenize(const char* delimiter=" \t\n", enTokenizeMode tokenizemode=tkCharSet);
        bool nextToken(GStr& token);
        int asInt(int base=10);
        double asReal();
        double asDouble() { return asReal(); }
        bool asReal(double& r);
        bool asDouble(double& r) { return asReal(r); }
        bool asInt(int& r, int base=10);
        int index(const GStr& s, int start_index = 0) const;
        int index(const char* s, int start_index = 0) const;
        int index(char c, int start_index = 0) const;
        int rindex(char c, int end_index = -1) const;
        int rindex(const char* str, int end_index = -1) const;
        bool contains(const GStr& s) const;
        bool contains(const char* s) const;
        bool contains(char c) const;
        bool startsWith(const char* s) const;
        bool startsWith(const GStr& s) const;
        bool endsWith(const char* s) const;
        bool endsWith(const GStr& s) const;
        GStr split(const char* delim);
        GStr split(char c);
           /* splits "this" in two parts, at the first (leftmost)
                 encounter of delim:
                 1st would stay in "this"
                 (which this way is truncated)
                 2nd will go to the returned string
           */
        GStr splitr(const char* delim);
        GStr splitr(char c);
           /* splits "this" in two parts, at the last (rightmost)
                 encounter of delim:
                 1st would stay in "this"
                 2nd will be returned
           */

        int peelInt() const; //extract an integer, (left to right), from a
                //mixed alphanumeric string, e.g. 'T24HC1234b'=> 2
        int peelIntR() const; //same as above, but starts from the right side
        //e.g. 'T2HC1234b'=> 1234
        GStr& trim(char c);
        GStr& trim(const char* c=" \t\n\r"); //trim both ends of characters in given set
        GStr& trimR(const char* c=" \t\n\r"); //trim only right end
        GStr& trimR(char c=' ');
        GStr& chomp(char c='\n') { return trimR(c); }
        GStr& chomp(const char* cstr); //like trimR, but given string is taken as a whole
        GStr& trimL(const char* c=" \t\n\r"); //trim only left end
        GStr& trimL(char c=' ');
        GStr& padR(uint len, char c=' '); //align it in len spaces to the right
        GStr& padL(uint len, char c=' '); //align it in len spaces to the left
        GStr& padC(uint len, char c=' '); //center it
        size_t read(FILE* stream, const char* delimiter="\n", size_t bufsize=4096);
          //read next token from stream, using the given string as
          //a marker where the block should stop
        const char* chars() const;
        const char* text() const;
        char* detach(); //returns pointer to the string, giving up on its memory management
    protected:
        char* fTokenDelimiter;
        int fLastTokenStart;
        enTokenizeMode fTokenizeMode;
        void* readbuf; //file read buffer for the read() function
        size_t readbufsize; //last setting for the readbuf
        static void invalid_args_error(const char* fname);
        static void invalid_index_error(const char* fname);
        struct Data {//structure holding actual
                     //string data and reference count information
               Data():ref_count(0), cap(0),length(0) { chars[0] = 0; }
               uint ref_count; //reference count
               uint cap; //allocated string capacity (excluding \0 end char)
               uint length; //actual string length (excluding \0 end char)
               char chars[1];
              };
        static Data* new_data(uint len, uint addcap=0); //alloc a specified length string's Data
        static Data* new_data(const char* str, uint addcap=0); //alloc a copy of a specified string, with an additional cap
        void prep_data(uint len, uint addcap=0); //allocates memory for the string, if needed
        void replace_data(Data* data);
        //WARNING (dangerous): direct access to pointer; string editing cannot change the length!
        char* chrs();
        void make_unique();
        static Data null_data; //a null (empty) string Data is available here
        Data* my_data; //pointer to a Data object holding actual string data
};

/***************************************************************************/

inline int GStr::length() const {
 return my_data->length;
 }


inline const char *GStr::chars() const {
 return my_data->chars;
 }

inline char *GStr::chrs() { //allows direct modification of the chars !
 return my_data->chars;
 }

inline const char *GStr::text() const {
 return my_data->chars;
 }

inline bool operator>=(const char *s1, const GStr& s2) {
 return (strcmp(s1, s2.chars()) >= 0);
 }

inline bool operator!=(const char *s1, const GStr& s2) {
 return (strcmp(s1, s2.chars()) != 0);
 }

inline void Gswap(GStr& s1, GStr& s2) {
 GStr::Data *tmp = s1.my_data; s1.my_data = s2.my_data;
 s2.my_data = tmp;
 }

#endif
