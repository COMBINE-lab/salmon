/*
 * Copyright (c) 1994 MEDICAL RESEARCH COUNCIL
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
 * File: array.h
 * Version:
 *
 * Description:
 *
 * Created:
 * Updated:
 *
 */

#ifndef _ARRAY_H_
#define _ARRAY_H_

/* 12/1/99 johnt - use stddef.h not sys/types.h for size_t */
#include <stddef.h>		/* IMPORT: size_t */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t size;		/* element size */
    size_t dim;			/* allocated number of elements */
    size_t max;			/* elements accessed */
    void *base;			/* base address of array */
} ArrayStruct, *Array;



extern Array ArrayCreate(size_t size, size_t dim);

extern int ArrayExtend(Array a, size_t dim);

extern void *ArrayRef(Array a, size_t i);

extern int ArrayDestroy(Array a);

#define ArrayMax(a) ( (a)->max )

#define ArrayBase(t,a) ( (t *)((a)->base) )

/*
#define arr(t,a,n) \
    (*(t*)((a)->base + (a)->size*(n)))

#define arrp(t,a,n) \
    ((t*)((a)->base + (a)->size*(n)))
*/




#define arr(t,a,n) \
    ((t*)((a)->base))[n]

#define ARR(t,a,n) \
    (*((t*)ArrayRef((a),(n))))

#define arrp(t,a,n) \
    &((t*)((a)->base))[n]
#define ARRP(t,a,n) \
    ((t*)ArrayRef(a,n))

#define ARRAY_NO_ERROR			 0
#define ARRAY_FULL			-1
#define ARRAY_INVALID_ARGUMENTS		-2
#define ARRAY_OUT_OF_MEMORY		-3

extern int ArrayError;

extern char *ArrayErrorString(int error);


#ifdef __cplusplus
}
#endif

#endif /*_ARRAY_H_*/
