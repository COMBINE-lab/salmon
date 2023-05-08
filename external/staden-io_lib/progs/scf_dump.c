/*
 * Copyright (c) 2007, 2010 Genome Research Ltd.
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

#include <stdio.h>
#include <io_lib/scf.h>

int main(int argc, char **argv) {
    Scf *scf;
    int i;

    if (argc != 2) {
	fprintf(stderr, "Usage: scf_dump scf_filename\n");
	return 1;
    }

    scf = read_scf(argv[1]);

    if (scf == NULL) {
	fprintf(stderr, "read_scf failed\n");
	return 1;
    }

    printf("[Header]\n");
    printf("%d\t# magic_number\n",	scf->header.magic_number);
    printf("%d\t\t# samples\n",		scf->header.samples);
    printf("%d\t\t# samples_offset\n",	scf->header.samples_offset);
    printf("%d\t\t# bases\n",		scf->header.bases);
    printf("%d\t\t# bases_left_clip\n",	scf->header.bases_left_clip);
    printf("%d\t\t# bases_right_clip\n",	scf->header.bases_right_clip);
    printf("%d\t\t# bases_offset\n",	scf->header.bases_offset);
    printf("%d\t\t# comments_size\n",	scf->header.comments_size);
    printf("%d\t\t# comments_offset\n",	scf->header.comments_offset);
    printf("%.4s\t\t# version\n",	        scf->header.version);
    printf("%d\t\t# sample_size\n",	scf->header.sample_size);
    printf("%d\t\t# code_set\n",		scf->header.code_set);
    printf("%d\t\t# private_size\n",	scf->header.private_size);
    printf("%d\t\t# private_offset\n",	scf->header.private_offset);
    printf("%d\t\t# spare[0]\n",		scf->header.spare[0]);
    printf("%d\t\t# spare[1]\n",		scf->header.spare[1]);
    printf("%d\t\t# spare[2]\n",		scf->header.spare[2]);
    printf("%d\t\t# spare[3]\n",		scf->header.spare[3]);
    printf("%d\t\t# spare[4]\n",		scf->header.spare[4]);
    printf("%d\t\t# spare[5]\n",		scf->header.spare[5]);
    printf("%d\t\t# spare[6]\n",		scf->header.spare[6]);
    printf("%d\t\t# spare[7]\n",		scf->header.spare[7]);
    printf("%d\t\t# spare[8]\n",		scf->header.spare[8]);
    printf("%d\t\t# spare[9]\n",		scf->header.spare[9]);
    printf("%d\t\t# spare[10]\n",		scf->header.spare[10]);
    printf("%d\t\t# spare[11]\n",		scf->header.spare[11]);
    printf("%d\t\t# spare[12]\n",		scf->header.spare[12]);
    printf("%d\t\t# spare[13]\n",		scf->header.spare[13]);
    printf("%d\t\t# spare[14]\n",		scf->header.spare[14]);
    printf("%d\t\t# spare[15]\n",		scf->header.spare[15]);
    printf("%d\t\t# spare[16]\n",		scf->header.spare[16]);
    printf("%d\t\t# spare[17]\n",		scf->header.spare[17]);

    puts("\n[Bases]");
    for (i = 0; i < scf->header.bases; i++) {
	printf("%c %05d %03d %03d %03d %03d   %03d %03d %03d  #%3d\n",
	       scf->bases[i].base,
	       scf->bases[i].peak_index,
	       scf->bases[i].prob_A,
	       scf->bases[i].prob_C,
	       scf->bases[i].prob_G,
	       scf->bases[i].prob_T,
	       scf->bases[i].spare[0],
	       scf->bases[i].spare[1],
	       scf->bases[i].spare[2],
	       i);
    }

    puts("\n[A_Trace]");
    if (scf->header.sample_size == 1) {
	for (i = 0; i < scf->header.samples; i++)
	    printf("%d\t#%5d\n", scf->samples.samples1[i].sample_A, i);
    } else {
	for (i = 0; i < scf->header.samples; i++)
	    printf("%d\t#%5d\n", scf->samples.samples2[i].sample_A, i);
    }
 
    puts("\n[C_Trace]");
    if (scf->header.sample_size == 1) {
	for (i = 0; i < scf->header.samples; i++)
	    printf("%d\t#%5d\n", scf->samples.samples1[i].sample_C, i);
    } else {
	for (i = 0; i < scf->header.samples; i++)
	    printf("%d\t#%5d\n", scf->samples.samples2[i].sample_C, i);
    }
 
    puts("\n[G_Trace]");
    if (scf->header.sample_size == 1) {
	for (i = 0; i < scf->header.samples; i++)
	    printf("%d\t#%5d\n", scf->samples.samples1[i].sample_G, i);
    } else {
	for (i = 0; i < scf->header.samples; i++)
	    printf("%d\t#%5d\n", scf->samples.samples2[i].sample_G, i);
    }
 
    puts("\n[T_Trace]");
    if (scf->header.sample_size == 1) {
	for (i = 0; i < scf->header.samples; i++)
	    printf("%d\t#%5d\n", scf->samples.samples1[i].sample_T, i);
    } else {
	for (i = 0; i < scf->header.samples; i++)
	    printf("%d\t#%5d\n", scf->samples.samples2[i].sample_T, i);
    }

    puts("\n[Comments]");
    printf("%.*s\n", (int)scf->header.comments_size, scf->comments);

    scf_deallocate(scf);

    return 0;
}
