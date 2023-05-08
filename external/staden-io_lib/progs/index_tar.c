/*
 * Copyright (c) 2004, 2007, 2009-2010, 2013 Genome Research Ltd.
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
 * Copyright (c) 2001 MEDICAL RESEARCH COUNCIL
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

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <io_lib/tar_format.h>

int main(int argc, char **argv) {
    int directories = 0;
    FILE *fp;
    tar_block blk;
    char member[257];
    size_t size, extra;
    int LongLink = 0;
    size_t offset = 0;
    
    /* process command line arguments of the form -arg */
    for (argc--, argv++; argc > 0; argc--, argv++) {
	if (**argv != '-' || strcmp(*argv, "--") == 0)
	    break;

	if (strcmp(*argv, "-d") == 0)
	    directories = 1;
    }

    if (argc != 1) {
	fprintf(stderr, "Usage: index_tar [-d] tarfile > tarfile.index\n");
	return 1;
    }

    /* open the tarfile */
    if (NULL == (fp = fopen(argv[0], "rb"))) {
	perror(argv[0]);
	return 1;
    }

    while(fread(&blk, sizeof(blk), 1, fp) == 1) {
	/*
	 * If a directory is too large to fit in the name (>100) but short
	 * enough to fit in the prefix the name field will be empty, this is
	 * not the case for ordinary files where the name field is always
	 * non-empty
	 */
	if (!blk.header.name[0] && !blk.header.prefix[0])
	    break;

        /* get size of member, rounded to a multiple of TBLOCK */
	size = strtoul(blk.header.size, NULL, 8);
        extra = TBLOCK*((size+TBLOCK-1)/TBLOCK) - size;

        /* skip directories unless requested */
        if (directories || blk.header.typeflag != DIRTYPE || LongLink) {

            /*
	     * extract member name (prefix + name), unless last member
	     * was ././@LongLink
	     */
            if (LongLink == 0) {
		char *cp;
                (void) strncpy(member, blk.header.prefix, 155);
		member[155] = 0;
	        if (strlen(blk.header.prefix) > 0 && blk.header.name[0])
		    (void) strcat(member, "/");
		cp = member + strlen(member);
    	        (void) strncpy(cp, blk.header.name, 100);
		cp[100] = 0;
            }
            
            /* account for gtar ././@LongLink */
            if (strcmp(member, "././@LongLink") == 0) {
                /* still expect filenames to fit into 256 bytes */
                if (size > 256) {
                    int r = fread(member, 1, 256, fp);
                    fprintf(stderr,"././@LongLink too long size=%ld\n",
			    (long)size);
		    if (r > 0)
			fprintf(stderr,"%s...\n", member);
                    exit(1);
                }
                /*
		 * extract full name of next member then rewind to start
		 * of header
		 */
                if (size != fread(member, 1, size > 256 ? 256 : size, fp))
		    exit(1);
                fseek(fp, -size, SEEK_CUR);
                LongLink = 1;
            } else {
                /* output offset, member name */
                printf("%lu %.256s\n", (long)offset, member);
                LongLink = 0;
            }
        }

        /* increment offset */
        size += extra;
        fseek(fp, size, SEEK_CUR);
        offset += sizeof(blk) + size;
    }

    fclose(fp);
    return 0;
}
