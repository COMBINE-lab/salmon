/*
 * Copyright (c) 2014 Genome Research Ltd.
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
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2014
 */

#include <io_lib/binning.h>

/* See http://res.illumina.com/documents/products/whitepapers/whitepaper_datacompression.pdf */
unsigned int illumina_bin[256] = {
     0,                                     /* 0 reserved for N */
     1,                                     /* Unused, but for completeness */
     6,  6,  6,  6,  6,  6,  6, 6,          /* 2-9 */
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, /* 10-19 */
    22, 22, 22, 22, 22,                     /* 20-24 */
    27, 27, 27, 27, 27,                     /* 25-29 */
    33, 33, 33, 33, 33,                     /* 30-34 */
    37, 37, 37, 37, 37,                     /* 35-39 */
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40, /* 40+ */
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40,
};

unsigned int illumina_bin_33[256] = {
     0+33,
     1+33,
     6+33,  6+33,  6+33,  6+33,  6+33,  6+33,  6+33,  6+33,
    15+33, 15+33, 15+33, 15+33, 15+33, 15+33, 15+33, 15+33, 15+33, 15+33,
    22+33, 22+33, 22+33, 22+33, 22+33,
    27+33, 27+33, 27+33, 27+33, 27+33,
    33+33, 33+33, 33+33, 33+33, 33+33,
    37+33, 37+33, 37+33, 37+33, 37+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
    40+33, 40+33, 40+33, 40+33, 40+33, 40+33,
};
