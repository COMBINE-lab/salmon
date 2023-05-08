/*
 * Copyright (c) 2007-2010 Genome Research Ltd.
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

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <io_lib/hash_table.h>
#include <io_lib/srf.h>

static char qlookup[256];
void init_qlookup(void) {
    int i;
    for (i = -128; i < 128; i++) {
        qlookup[i+128] = '!' + (int)((10*log(1+pow(10, i/10.0))/log(10)+.499));
    }
}

/* ------------------------------------------------------------------------ */

#define MAX_READ_LEN 10000
void ztr2fastq(ztr_t *z, char *name, int calibrated) {
    int i, nc, seq_len;
    char buf[MAX_READ_LEN*2 + 512 + 6];
    char *seq, *qual, *sdata, *qdata;
    ztr_chunk_t **chunks;

    /* Extract the sequence only */
    chunks = ztr_find_chunks(z, ZTR_TYPE_BASE, &nc);
    if (nc != 1) {
	fprintf(stderr, "Zero or greater than one BASE chunks found.\n");
	return;
    }
    uncompress_chunk(z, chunks[0]);
    sdata = chunks[0]->data+1;
    seq_len = chunks[0]->dlength-1;

    /* Extract the quality */
    if (calibrated)
	chunks = ztr_find_chunks(z, ZTR_TYPE_CNF1, &nc);
    else
	chunks = ztr_find_chunks(z, ZTR_TYPE_CNF4, &nc);

    if (nc != 1) {
	fprintf(stderr, "Zero or greater than one CNF chunks found.\n");
	return;
    }
    uncompress_chunk(z, chunks[0]);
    qdata = chunks[0]->data+1;

    /* Construct fastq entry */
    seq = buf;
    *seq++ = '@';
    while (*name)
	*seq++ = *name++;
    *seq++ = '\n';

    qual = seq + seq_len;
    *qual++ = '\n';
    *qual++ = '+';
    *qual++ = '\n';

    for (i = 0; i < seq_len; i++) {
	if (*sdata != '.') {
	    *seq++ = *sdata++;
	} else {
	    *seq++ = 'N';
	    sdata++;
	}
	*qual++ = calibrated
	    ? *qdata++ + '!'
	    : qlookup[*qdata++ + 128];
    }
    *qual++ = '\n';

    fwrite(buf, 1, qual - buf, stdout);
}

/* ------------------------------------------------------------------------ */
void usage(void) {
    fprintf(stderr, "Usage: srf_extract [-fastq] [-c] archive_name trace_name ...\n");
    exit(1);
}

int main(int argc, char **argv) {
    FILE *fp;
    srf_t *srf;
    char *archive, *trace;
    uint64_t cpos, hpos, dpos;
    int fastq = 0, calibrated = 0, i;

    /* Parse args */
    for (i = 1; i < argc && argv[i][0] == '-'; i++) {
	if (!strcmp(argv[i], "-")) {
	    break;
	} else if (!strcmp(argv[i], "-c")) {
	    calibrated = 1;
	} else if (!strcmp(argv[i], "-fastq")) {
	    fastq = 1;
	} else {
	    usage();
	}
    }    

    if (argc < (i+2)) {
      usage();
    }

    /* the SRF archive */
    archive = argv[i++];
    fp = fopen(archive, "rb");
    if (NULL == fp) {
        perror(archive);
        return 1;
    }
    srf = srf_create(fp);

    /* the trace */
    trace = argv[i];

    if( fastq ){
        read_sections(READ_BASES);
        init_qlookup();
    }else{
#ifdef _WIN32
        _setmode(_fileno(stdout), _O_BINARY);
#endif
    }

    /* Search index */
    switch (srf_find_trace(srf, trace, &cpos, &hpos, &dpos)) {
    case -1:
        fprintf(stderr, "Malformed or missing index hash. "
                "Consider running srf_index_hash\n");
        return 1;

    case -2:
        fprintf(stderr, "%s: not found\n", trace);
        break;

    default:
        /* The srf object holds the latest data and trace header blocks */
        if( fastq ){
            mFILE *mf = mfcreate(NULL, 0);
            mfwrite(srf->th.trace_hdr, 1, srf->th.trace_hdr_size, mf);
            mfwrite(srf->tb.trace,     1, srf->tb.trace_size,     mf);
            mfseek(mf, 0, SEEK_SET);
            ztr_t *ztr = partial_decode_ztr(srf, mf, NULL);
            ztr2fastq(ztr, trace, calibrated);
            delete_ztr(ztr);
            mfdestroy(mf);
        } else {
            fwrite(srf->th.trace_hdr, 1, srf->th.trace_hdr_size, stdout);
            fwrite(srf->tb.trace,     1, srf->tb.trace_size,     stdout);
        }
        break;
    }
	
    return 0;
}
