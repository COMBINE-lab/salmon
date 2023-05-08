/*
 * Copyright (c) 2005, 2007, 2010 Genome Research Ltd.
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
 * Author(s): James Bonfield, Simon Dear
 * 
 * Copyright (c) 1992, 1995 MEDICAL RESEARCH COUNCIL
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

/*
 * Machine independant io:
 * For reading and writing to big-endian and little-endian files
 *
 * Routines available:
 *     be_write_int_1()
 *     be_write_int_2()
 *     be_write_int_4()
 *     be_write_int_8()
 *     be_read_int_1()
 *     be_read_int_2()
 *     be_read_int_4()
 *     be_read_int_8()
 *     le_write_int_1()
 *     le_write_int_2()
 *     le_write_int_4()
 *     le_write_int_8()
 *     le_read_int_1()
 *     le_read_int_2()
 *     le_read_int_4()
 *     le_read_int_8()
 *
 * All routine return:
 *    0 - an error has occurred during io operation
 *    1 - value suggessfully read or written
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include "io_lib/stdio_hack.h"
#include "io_lib/mach-io.h"




/**********************************************************************/
/* IO for big-endian files                                            */
/**********************************************************************/

/*
 * Write a big-endian int1
 */
int be_write_int_1(FILE *fp, uint1 *i1)
{
    if (fwrite(i1, sizeof(uint1), 1, fp) != 1) return (0);
    return (1);
}


/*
* Write a big-endian int2
*/
int be_write_int_2(FILE *fp, uint2 *i2)
{
    uint2 i = be_int2(*i2);
    
    if (fwrite(&i, 2, 1, fp) != 1) return (0);
    return (1);
}

/*
 * Write a big-endian int4
 */
int be_write_int_4(FILE *fp, uint4 *i4)
{
    uint4 i = be_int4(*i4);
    
    if (fwrite(&i, 4, 1, fp) != 1) return (0);

    return (1);
}

/*
 * Write a big-endian int8
 */
int be_write_int_8(FILE *fp, uint8 *i8)
{
    uint8 i = be_int8(*i8);
    
    if (fwrite(&i, 8, 1, fp) != 1) return (0);

    return (1);
}


/*
 * Read a big-endian int1
 */
int be_read_int_1(FILE *fp, uint1 *i1)
{
    if (fread(i1, sizeof(uint1), 1, fp) != 1) return (0);
    return (1);
}


/*
 * Read a big-endian int2
 */
int be_read_int_2(FILE *fp, uint2 *i2)
{
    uint2 i;
    
    if (fread(&i, 2, 1, fp) != 1) return (0);
    *i2 = be_int2(i);

    return (1);
}


/*
 * Read a big-endian int4
 */
int be_read_int_4(FILE *fp, uint4 *i4)
{
    uint4 i;
    
    if (fread(&i, 4, 1, fp) != 1) return (0);
    *i4 = be_int4(i);

    return (1);
}


/*
 * Read a big-endian int8
 */
int be_read_int_8(FILE *fp, uint8 *i8)
{
    uint8 i;
    
    if (fread(&i, 8, 1, fp) != 1) return (0);
    *i8 = be_int8(i);

    return (1);
}





/**********************************************************************/
/* IO for little-endian files                                         */
/**********************************************************************/

/*
 * Write a little-endian int1
 */
int le_write_int_1(FILE *fp, uint1 *i1)
{
    if (fwrite(i1, sizeof(uint1), 1, fp) != 1) return (0);
    return (1);
}


/*
 * Write a little-endian int2
 */
int le_write_int_2(FILE *fp, uint2 *i2)
{
    uint2 i = le_int2(*i2);
    
    if (fwrite(&i, 2, 1, fp) != 1) return (0);

    return (1);
}


/*
 * Write a little-endian int4
 */
int le_write_int_4(FILE *fp, uint4 *i4)
{
    uint4 i = le_int4(*i4);
    
    if (fwrite(&i, 4, 1, fp) != 1) return (0);

    return (1);
}


/*
 * Write a little-endian int8
 */
int le_write_int_8(FILE *fp, uint8 *i8)
{
    uint8 i = le_int8(*i8);
    
    if (fwrite(&i, 8, 1, fp) != 1) return (0);

    return (1);
}


/*
 * Read a little-endian int1
 */
int le_read_int_1(FILE *fp, uint1 *i1)
{
    if (fread(i1, sizeof(uint1), 1, fp) != 1) return (0);
    return (1);
}


/*
 * Read a little-endian int2
 */
int le_read_int_2(FILE *fp, uint2 *i2)
{
    uint2 i;
    
    if (fread(&i, 2, 1, fp) != 1) return (0);
    *i2 = le_int2(i);

    return (1);
}

/*
 * Read a little-endian int4
 */
int le_read_int_4(FILE *fp, uint4 *i4)
{
    uint4 i;
    
    if (fread(&i, 4, 1, fp) != 1) return (0);
    *i4 = le_int4(i);

    return (1);
}


/*
 * Read a little-endian int8
 */
int le_read_int_8(FILE *fp, uint8 *i8)
{
    uint8 i;
    
    if (fread(&i, 8, 1, fp) != 1) return (0);
    *i8 = le_int8(i);

    return (1);
}
