/* MRC local changes made by James Bonfield. Derived from tar.h: */

/*-
 * Copyright (c) 1992 Keith Muller.
 * Copyright (c) 1992, 1993
 * The Regents of the University of California.  All rights reserved.
 *
 * This code is derived from software contributed to Berkeley by
 * Keith Muller of the University of California, San Diego.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 *@(#)tar.h8.2 (Berkeley) 4/18/94
 * $FreeBSD$
 */

#ifndef _TAR_FORMAT_H
#define _TAR_FORMAT_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Our own tar block defines - we cannot rely on UNIX to provide these for us
 * as the sun tar.h is minimal and Alliant's does not even exist.
 */
#define TBLOCK  512
#define NAMSIZ  100

/* Values used in typeflag field. */
#define REGTYPE         '0'             /* Regular File */
#define AREGTYPE        '\0'            /* Regular File */
#define LNKTYPE         '1'             /* Hard Link */
#define SYMTYPE         '2'             /* Symbolic Link */
#define CHRTYPE         '3'             /* Character Special File */
#define BLKTYPE         '4'             /* Block Special File */
#define DIRTYPE         '5'             /* Directory */
#define FIFOTYPE        '6'             /* FIFO */
#define CONTTYPE        '7'             /* Reserved */

/*
 * There will usually be more data than this in a tar header - but we don't
 * need to concern ourselves with it.
 */
typedef union hblock {
    char data[TBLOCK];
    struct header {
	char name[NAMSIZ];
	char mode[8];
	char uid[8];
	char gid[8];
	char size[12];
	char mtime[12];
	char chksum[8];
	char typeflag;
	char linkname[NAMSIZ];
	char magic[6];
	char version[2];
	char uname[32];
	char gname[32];
	char devmajor[8];
	char devminor[8];
	char prefix[155];
    } header;
} tar_block;

#ifdef __cplusplus
}
#endif

#endif
