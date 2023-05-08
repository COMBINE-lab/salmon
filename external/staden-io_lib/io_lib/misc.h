/*
 * Copyright (c) 2007, 2010, 2013 Genome Research Ltd.
 * Author(s): James Bonfield
 * 
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 * 
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 * 
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *    Institute nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific
 *    prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Copyright (c) 1994-1997, 2001-2002 MEDICAL RESEARCH COUNCIL
 * All rights reserved
 * 
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 * 
 *    1 Redistributions of source code must retain the above copyright notice, 
 *      this list of conditions and the following disclaimer.
 * 
 *    2 Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 * 
 *    3 Neither the name of the MEDICAL RESEARCH COUNCIL, THE LABORATORY OF
 *      MOLECULAR BIOLOGY nor the names of its contributors may be used
 *      to endorse or promote products derived from this software without
 *      specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _misc_h
#define _misc_h

#include "io_lib/os.h"

#include <stdio.h>
#include <stdarg.h>  /* varargs needed for v*printf() prototypes */
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * This informs gcc that crash() doesn't return, so it doesn't need to
 * concern itself that code paths going via crash could mean some variables
 * being undefined and then issuing uninitialised variable warnings.
 * This particularly affected convert.
 */
#ifdef __GNUC__
#    define __NORETURN__ __attribute__ ((__noreturn__))
#else
#    define __NORETURN__
#endif

/*
 * Used for printf style argument checking. We can request a function such
 * as vTcl_SetResult does argument checking, avoiding bugs with using
 * %d and passing in a 64-bit record.
 */
#ifdef __GNUC__
#    define __PRINTF_FORMAT__(a,b) __attribute__ ((format (printf, a, b)))
#else
#    define __PRINTF_FORMAT__(a,b)
#endif

/*
 * Stop gcc from warning about unused variables.
 */

#ifdef __GNUC__
#    define __UNUSED__  __attribute__ ((__unused__))
#else
#    define __UNUSED__
#endif


/* from files.c */
#ifdef _WIN32
FILE *tmpfile_win(void);
#endif

extern int is_directory(char * fn);
extern int is_file(char * fn);
extern int file_exists(char * fn);
extern int compressed_file_exists(char *fname);
extern int file_size(char * fn);
extern FILE *open_fofn(char *files);
extern char *read_fofn(FILE *fp);
extern void close_fofn(FILE *fp);
extern int fstrlen(char *f, int max_f);
extern void f2cstr(char *f, int max_f, char *c, int max_c);
extern void c2fstr(char *c, int max_c, char *f, int max_f);
extern char *mystrtok(char *s, char *ct);
extern char *myfind(char *file, char* searchpath, int (*found) (char *) );
extern void crash (char* format,...);
extern void str_tolower (char *s);
extern void str_toupper (char *s);
extern char *fn_tail (char *s);
extern void fn_tolower (char *s);
extern void fn_toupper (char *s);
extern void shell_call(char *command, char *output, int len);
extern char *date_str(void);
#ifdef NOSTRDUP
extern char *strdup(const char *s);
#endif
#ifdef NOSTRSTR
extern char *strstr(char *cs, char *ct);
#endif

#ifdef NOMEMMOVE
#define memmove(d,s,l) bcopy(s,d,l)
#endif
extern int myusleep(unsigned int useconds);

extern void errout(char *fmt, ...);
extern void messout(char *fmt, ...);

/*
 * Useful macros
 */
#define findfile(F,S) myfind((F),(S),file_exists)
/*is_file fails for symbolic links*/
/*#define findfile(F,S) myfind((F),(S),is_file)*/

#ifdef MIN
#undef MIN
#endif
#define MIN(A,B) ( ( (A) < (B) ) ? (A) : (B) )
#ifdef MAX
#undef MAX
#endif
#define MAX(A,B) ( ( (A) > (B) ) ? (A) : (B) )
#define SGN(A) ( (A) ? ( ( (A) < 0 ) ? -1 : 1 ) : 0 )
#define ABS(A) ( (A) < 0 ? -(A) : (A) )

/* Number of elements in array */
#define Number(A) ( sizeof(A) / sizeof((A)[0]) )

/*
 * Things taken from the new gap text_output.h. They'll be used globally
 * across all the programs in the end.
 */

/*
 * Usage: verror(priority, format, args...);
 * NB: don't pass more than 8K per call
 */
#define ERR_WARN 0
#define ERR_FATAL 1
void verror(int priority, char *name, char *fmt, ...);

/*
 * Usage: vmessage(format, args...);
 * NB: don't pass more than 8K per call
 */
void vmessage(char *fmt, ...);

/*
 * Adds a new header to the text output window.
 */
void vfuncheader(char *fmt, ...);

/*
 * As vfuncheader, but only outputting when necessary.
 */
void vfuncgroup(int group, char *fmt, ...);

#ifdef __cplusplus
}
#endif

#endif /*_misc_h*/
