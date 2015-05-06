/*
 * Copyright (c) 2005, 2007 Genome Research Ltd.
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
 * Author(s): James Bonfield
 * 
 * Copyright (c) 1994-1995 MEDICAL RESEARCH COUNCIL
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

/*
 * Copyright (c) Medical Research Council 1994. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * as part of the Staden Package at the MRC Laboratory of Molecular
 * Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

#ifndef _mach_io_h
#define _mach_io_h
/*
 * Machine independant io
 * For reading and writing to big-endian and little-endian files
 *
 * Routines available:
 *     be_write_int_1()
 *     be_write_int_2()
 *     be_write_int_4()
 *     be_read_int_1()
 *     be_read_int_2()
 *     be_read_int_4()
 *     le_write_int_1()
 *     le_write_int_2()
 *     le_write_int_4()
 *     le_read_int_1()
 *     le_read_int_2()
 *     le_read_int_4()
 *
 * All routine return:
 *    0 - an error has occurred during io operation
 *    1 - value successfully read or written
 */

#include <stdio.h>
#include "io_lib/os.h"
#include "io_lib/mFILE.h"

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/
/* IO for big-endian files                                            */
/**********************************************************************/

/*
 * Write a big-endian int1
 */
extern int be_write_int_1(mFILE *fp, uint1 *i1);

/*
 * Write a big-endian int2
 */
extern int be_write_int_2(mFILE *fp, uint2 *i2);

/*
 * Write a big-endian int4
 */
extern int be_write_int_4(mFILE *fp, uint4 *i4);

/*
 * Write a big-endian int8
 */
extern int be_write_int_8(mFILE *fp, uint8 *i8);

/*
 * Read a big-endian int1
 */
extern int be_read_int_1(mFILE *fp, uint1 *i1);

/*
 * Read a big-endian int2
 */
extern int be_read_int_2(mFILE *fp, uint2 *i2);

/*
 * Read a big-endian int4
 */
extern int be_read_int_4(mFILE *fp, uint4 *i4);

/*
 * Read a big-endian int8
 */
extern int be_read_int_8(mFILE *fp, uint8 *i8);

/**********************************************************************/
/* IO for little-endian files                                         */
/**********************************************************************/

/*
 * Write a little-endian int1
 */
extern int le_write_int_1(mFILE *fp, uint1 *i1);

/*
 * Write a little-endian int2
 */
extern int le_write_int_2(mFILE *fp, uint2 *i2);

/*
 * Write a little-endian int4
 */
extern int le_write_int_4(mFILE *fp, uint4 *i4);

/*
 * Write a little-endian int8
 */
extern int le_write_int_8(mFILE *fp, uint8 *i8);

/*
 * Read a little-endian int1
 */
extern int le_read_int_1(mFILE *fp, uint1 *i1);

/*
 * Read a little-endian int2
 */
extern int le_read_int_2(mFILE *fp, uint2 *i2);

/*
 * Read a little-endian int4
 */
extern int le_read_int_4(mFILE *fp, uint4 *i4);

/*
 * Read a little-endian int8
 */
extern int le_read_int_8(mFILE *fp, uint8 *i8);

#ifdef __cplusplus
}
#endif

#endif /* _mach_io_h */
