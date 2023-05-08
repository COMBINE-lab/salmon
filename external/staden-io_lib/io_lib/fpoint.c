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

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <math.h>
#include "io_lib/fpoint.h"
/*
extern double log ( double x ) ;
extern double exp ( double x ) ;
*/
#define IEEE

float int_to_float(int in)
/*
** interpret the integer in as a
** floating point number in IEEE format
*/
{
   /*
  Assume `in' is stored as a float according to the 
  ANSI IEEE 754-1985 standard. See the tables below:

  s = sign ( 1 bit)
  e = biased exponent (8 bits)
  f = fraction (23 bits)

  floating point number =  (-1)^s 2^(e-127) 1.f

     Bits  Name      Content
      31   Sign      1 iff number is negative
    23-30  Exponent  Eight-Bit exponent, biased by 127
     0-22  Fraction  23-bit fraction component of normalised significant.
		     The "one" bit is "hidden"

  If IEEE floating point format is supported on your machine...
  ensure there is a #define IEEE somewhere. 
  */

#ifdef IEEE
  union {
    int i;
    float f;
  } cvt;
  cvt.i = in;
  return cvt.f;
#else
  int fraction;
  int exponent;
  int sign;

  fraction = in & ( (1<<23)-1 );
  exponent = (in >> 23) & ( (1<<8)-1 );
  sign = (in >> 31);

  return
    (float) (
      (sign?-1.0:1.0) *
      exp ( log ( (double) 2.0) * (double) (exponent - 127 - 23) ) *
      (double) ((1<<23)+fraction)) ;
#endif
}
