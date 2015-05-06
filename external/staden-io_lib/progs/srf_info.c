/*
 * Copyright (c) 2008-2010, 2013 Genome Research Ltd.
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
 * This performs a linear (non-indexed) search for a trace in an SRF archive.
 *
 * It's not intended as a suitable production program or as a library of code
 * to use, but as a test and benchmark statistic.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <io_lib/Read.h>
#include <io_lib/misc.h>
#include <io_lib/ztr.h>
#include <io_lib/srf.h>
#include <io_lib/hash_table.h>
#include <io_lib/xalloc.h>

#define LEVEL_READ  (1 << 0)
#define LEVEL_CHUNK (1 << 1)
#define LEVEL_NAME  (1 << 2)
#define LEVEL_BASE  (1 << 3)
#define LEVEL_ALL   (LEVEL_READ | LEVEL_CHUNK | LEVEL_NAME | LEVEL_BASE );

/* only checks the first 10 traces */
#define LEVEL_CHECK 255 

#define READ_GOOD  0
#define READ_BAD   1
#define READ_TOTAL 2

#define NREADS 3

#define READ_GOOD_STR   "GOOD"
#define READ_BAD_STR    "BAD"
#define READ_TOTAL_STR  "TOTAL"

/* see ztr.h for a list of all possible ztr chunk types */

#define CHUNK_BASE 0
#define CHUNK_CNF1 1
#define CHUNK_CNF4 2
#define CHUNK_SAMP 3
#define CHUNK_SMP4 4
#define CHUNK_REGN 5

#define NCHUNKS 6

#define CHUNK_BASE_TYPE ZTR_TYPE_BASE
#define CHUNK_CNF1_TYPE ZTR_TYPE_CNF1
#define CHUNK_CNF4_TYPE ZTR_TYPE_CNF4
#define CHUNK_SAMP_TYPE ZTR_TYPE_SAMP
#define CHUNK_SMP4_TYPE ZTR_TYPE_SMP4
#define CHUNK_REGN_TYPE ZTR_TYPE_REGN

#define KEY_TYPE    0
#define KEY_VALTYPE 1
#define KEY_GROUP   2
#define KEY_OFFS    3
#define KEY_SCALE   4
#define KEY_COORD   5
#define KEY_NAME    6

#define NKEYS 7

#define KEY_TYPE_STR    "TYPE"
#define KEY_VALTYPE_STR "VALTYPE"
#define KEY_GROUP_STR   "GROUP"
#define KEY_OFFS_STR    "OFFS"
#define KEY_SCALE_STR   "SCALE"
#define KEY_COORD_STR   "COORD"
#define KEY_NAME_STR    "NAME"

#define TYPE_PROC 0
#define TYPE_SLXI 1
#define TYPE_SLXN 2
#define TYPE_0FAM 3
#define TYPE_1CY3 4
#define TYPE_2TXR 5
#define TYPE_3CY5 6

#define NTYPES 7

#define TYPE_PROC_STR "PROC"
#define TYPE_SLXI_STR "SLXI"
#define TYPE_SLXN_STR "SLXN"
#define TYPE_0FAM_STR "0FAM"
#define TYPE_1CY3_STR "1CY3"
#define TYPE_2TXR_STR "2TXR"
#define TYPE_3CY5_STR "3CY5"

#define MAX_REGIONS   4

/* regn chunk */
typedef struct {
    char coord;
    char *region_names;
    int nregions;
    char *name[MAX_REGIONS];
    char code[MAX_REGIONS];
    int start[MAX_REGIONS];
    int length[MAX_REGIONS];
    int count;
} regn_t;

/* ------------------------------------------------------------------------ */
/*
 * Print usage message to stderr and exit with the given \"code\".
 */
void usage(int code) {
    fprintf(stderr, "Usage: srf_info [-level level_bitmap] input(s)\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -l level_bitmap \n");
    fprintf(stderr, "              1\tCount of good/bad reads.\n");
    fprintf(stderr, "              2\tCounts for selected chunk types.\n");
    fprintf(stderr, "              4\tTrace count and trace name prefix for each trace_header.\n");
    fprintf(stderr, "              8\tBase count.\n");
    fprintf(stderr, "\n");

    exit(code);
}

/*
 * Parse the REGN chunk, add to regn HASH
 *
 * Returns corresponding HashItem * from regn Hash
 */
HashItem *parse_regn(ztr_t *z, ztr_chunk_t *chunk, HashTable *regn_hash) {
    char key[1024];
    char *name;
    HashItem *hi;
    regn_t *regn;
    size_t l;
    
    uncompress_chunk(z, chunk);

    /* the hash key is a combination of the region names and boundaries */
    name = ztr_lookup_mdata_value(z, chunk, "NAME");
    l = snprintf(key, sizeof(key), "names=%s", name);
    if( chunk->dlength ){
        int nbndy = (chunk->dlength-1)/4;
        uint4 *bndy = (uint4 *)(chunk->data+1);
        int ibndy;
	for (ibndy=0; ibndy<nbndy; ibndy++) {
            if( ibndy )
                l += snprintf(key + l, sizeof(key) - l,
			      ";%d", be_int4(bndy[ibndy]));
            else
                l += snprintf(key + l, sizeof(key) - l,
			      " boundaries=%d", be_int4(bndy[ibndy]));
        }
    }

    if (NULL == (hi = (HashTableSearch(regn_hash, key, strlen(key))))) {
        int iregion, nregions = 0;
        char *coord;
	char *cp1;
        uint4 bndy[MAX_REGIONS];
        int ibndy, nbndy = 0;
        HashData hd;

        if( NULL == (regn = (regn_t *)malloc(sizeof(regn_t)))) {
	    return NULL;
	}

	coord = ztr_lookup_mdata_value(z, chunk, "COORD");
	regn->coord = (NULL == coord ? 'B' : *coord );

	regn->region_names = strdup(name);

        cp1 = strtok (regn->region_names,";");
        while(cp1) {
            char *cp2;
            if(NULL == (cp2 = strchr(cp1,':'))) {
                fprintf(stderr, "Invalid region name/code pair %s\n", cp1);
                return NULL;
            }
            *cp2++ = '\0';
            regn->name[nregions] = cp1;
            regn->code[nregions] = *cp2;
            nregions++;
            cp1 = strtok (NULL, ";");
        }

        regn->nregions = nregions;

	if( chunk->dlength ) {
            nbndy = (chunk->dlength-1)/4;
            memcpy(bndy, chunk->data+1, chunk->dlength-1);
	}

        for( iregion=0, ibndy=0; iregion<nregions; iregion++) {
            /* start = (start + length of previous region) or 0 if no previous region */
            /* length = (next boundary - start of region) or -1 if no next boundary */
            if( regn->code[iregion] == 'E' ){
                /* no sequence, length = 0 */
                regn->start[iregion] = (iregion ? (regn->start[iregion-1] + regn->length[iregion-1]) : 0);
                regn->length[iregion] = 0;
            }else{
                if( ibndy > nbndy ){
                    fprintf(stderr, "More name/code pairs than boundaries\n");
                    return NULL;
                }
                regn->start[iregion] = (iregion ? (regn->start[iregion-1] + regn->length[iregion-1]) : 0);
                regn->length[iregion] = (ibndy == nbndy ? -1 : (be_int4(bndy[ibndy])-regn->start[iregion]));
                ibndy++;
            }
        }

        regn->count = 1;
            
	hd.p = regn;
	if (NULL == (hi = HashTableAdd(regn_hash, key, strlen(key), hd, NULL))) {
	    free(regn->region_names);
	    free(regn);
	    return NULL;
	}
    } else {
	regn = (regn_t *)(hi->data.p);
	regn->count++;
    }

    return hi;
}

/*
 * Parse the sequence
 *
 * Returns 0 on success.
 */
int parse_base(ztr_t *z, ztr_chunk_t *chunk, uint64_t *base_count) {
    int i;
    
    uncompress_chunk(z, chunk);

    for (i = 1; i < chunk->dlength; i++) {
	char base = chunk->data[i];
        uint1 key;
	switch (base) {
	case 'A': case 'a':
            key = 0;
            break;
	case 'C': case 'c':
            key = 1;
            break;
	case 'G': case 'g':
            key = 2;
            break;
	case 'T': case 't':
            key = 3;
            break;
	default:
            key = 4;
            break;
	}
        base_count[key]++;
    }

    return 0;
}

/*
 * count the mdata keys
 *
 * Returns 0 on success.
 */
int count_mdata_keys(ztr_t *z, ztr_chunk_t *chunk, int ichunk, long key_count[NCHUNKS][NKEYS], long type_count[NCHUNKS][NTYPES]) {
    char *keys_str[] = {KEY_TYPE_STR, KEY_VALTYPE_STR, KEY_GROUP_STR, KEY_OFFS_STR, KEY_SCALE_STR, KEY_COORD_STR, KEY_NAME_STR};
    char *types_str[] = {TYPE_PROC_STR, TYPE_SLXI_STR, TYPE_SLXN_STR, TYPE_0FAM_STR, TYPE_1CY3_STR, TYPE_2TXR_STR, TYPE_3CY5_STR};
    int ikey, itype;

    if (z->header.version_major > 1 ||
	z->header.version_minor >= 2) {
	/* ZTR format 1.2 onwards */

	char *cp = chunk->mdata;
	int32_t dlen = chunk->mdlength;

	/*
	 * NB: we may wish to rewrite this using a dedicated state machine
	 * instead of strlen/strcmp as this currently assumes the meta-
	 * data is correctly formatted, which we cannot assume as the 
	 * metadata is external and outside of our control.
	 * Passing in non-nul terminated strings could crash this code.
	 */
	while (dlen > 0) {
	    size_t l;

	    /* key */
	    l = strlen(cp);
	    for (ikey=0; ikey<NKEYS; ikey++)
		if(0 == strcmp(cp, keys_str[ikey]))
		    break;

	    cp += l+1;
	    dlen -= l+1;

	    /* value */
	    if (ikey < NKEYS)
		key_count[ichunk][ikey]++;

	    /* for the type key check the value */
	    if (ikey == KEY_TYPE && (ichunk == CHUNK_SAMP || ichunk == CHUNK_SMP4)) {
		for (itype=0; itype<NTYPES; itype++)
		    if(0 == strcmp(cp, types_str[itype]))
			break;
		if(itype < NTYPES)
		    type_count[ichunk][itype]++;
	    }

	    l = strlen(cp);
	    cp += l+1;
	    dlen -= l+1;
	}

    } else {
	/* v1.1 and before only supported a few types, specifically coded
	 * per chunk type.
	 */

	switch (chunk->type) {
	case ZTR_TYPE_SAMP:
	case ZTR_TYPE_SMP4:
	    key_count[ichunk][KEY_TYPE]++;
	    for (itype=0; itype<NTYPES; itype++)
		if(0 == strcmp(chunk->mdata, types_str[itype])) {
		    type_count[ichunk][itype]++;
		}
	    break;

	default:
	    break;
	}
    }

    return 0;
}

/*
 * As per partial_decode_ztr in srf.c, but without uncompress_ztr.
 */
ztr_t *partial_decode_ztr2(srf_t *srf, mFILE *mf, ztr_t *z) {
    ztr_t *ztr;
    ztr_chunk_t *chunk;
    long pos = 0;

    if (z) {
	/* Use existing ZTR object => already loaded header */
	ztr = z;

    } else {
	/* Allocate or use existing ztr */
	if (NULL == (ztr = new_ztr()))
	    return NULL;

	/* Read the header */
	if (-1 == ztr_read_header(mf, &ztr->header)) {
	    if (!z)
		delete_ztr(ztr);
	    mrewind(mf);
	    return NULL;
	}

	/* Check magic number and version */
	if (memcmp(ztr->header.magic, ZTR_MAGIC, 8) != 0) {
	    if (!z)
		delete_ztr(ztr);
	    mrewind(mf);
	    return NULL;
	}

	if (ztr->header.version_major != ZTR_VERSION_MAJOR) {
	    if (!z)
		delete_ztr(ztr);
	    mrewind(mf);
	    return NULL;
	}
    }

    /* Load chunks */
    pos = mftell(mf);
    while ((chunk = ztr_read_chunk_hdr(mf))) {
	chunk->data = (char *)xmalloc(chunk->dlength);
	if (chunk->dlength != mfread(chunk->data, 1, chunk->dlength, mf))
	    break;

	ztr->nchunks++;
	ztr->chunk = (ztr_chunk_t *)xrealloc(ztr->chunk, ztr->nchunks *
					     sizeof(ztr_chunk_t));
	memcpy(&ztr->chunk[ztr->nchunks-1], chunk, sizeof(*chunk));
	xfree(chunk);
	pos = mftell(mf);
    }

    /*
     * At this stage we're 'pos' into the mFILE mf with any remainder being
     * a partial block.
     */
    if (0 == ztr->nchunks) {
	if (!z)
	    delete_ztr(ztr);
	mrewind(mf);
	return NULL;
    }

    /* Ensure we exit at the start of a ztr CHUNK */
    mfseek(mf, pos, SEEK_SET);

    /* If this is the header part, ensure we uncompress and init. data */
    if (!z) {
	/* Force caching of huffman code_sets */
	ztr_find_hcode(ztr, CODE_USER);

	/* And uncompress the rest */
	uncompress_ztr(ztr);
    }

    return ztr;
}

/*
 * Given the archive name and the level_mode
 * generate information about the archive
 *
 * Note the generated srf file is NOT indexed
 *
 * Returns 0 on success.
 */
int srf_info(char *input, int level_mode, long *read_count, long *chunk_count,
	     uint64_t *chunk_size, long key_count[NCHUNKS][NKEYS],
	     long type_count[NCHUNKS][NTYPES], HashTable *regn_hash,
	     uint64_t *base_count) {
    srf_t *srf;
    off_t pos;
    int type;
    int count = 0;
    long trace_body_count = 0;
    char name[1024];

    if (NULL == (srf = srf_open(input, "rb"))) {
	perror(input);
	return 1;
    }

    while ((type = srf_next_block_type(srf)) >= 0) {

      switch (type) {
	case SRFB_CONTAINER:
	    if( trace_body_count ){
		if( level_mode & LEVEL_NAME )
		    printf( " ... %s x%ld\n", name+strlen(srf->th.id_prefix), trace_body_count);
		trace_body_count = 0;
	    }
	    if (0 != srf_read_cont_hdr(srf, &srf->ch)) {
		fprintf(stderr, "Error reading container header.\nExiting.\n");
		exit(1);
	    }
	    break;

        case SRFB_XML:
	    if( trace_body_count ){
		if( level_mode & LEVEL_NAME )
		    printf( " ... %s x%ld\n", name+strlen(srf->th.id_prefix), trace_body_count);
		trace_body_count = 0;
	    }
	    if (0 != srf_read_xml(srf, &srf->xml)) {
		fprintf(stderr, "Error reading XML.\nExiting.\n");
		exit(1);
	    }
	    break;

	case SRFB_TRACE_HEADER:
	    if( trace_body_count ){
		if( level_mode & LEVEL_NAME )
		    printf( " ... %s x%ld\n", name+strlen(srf->th.id_prefix), trace_body_count);
		trace_body_count = 0;
	    }
	    if (0 != srf_read_trace_hdr(srf, &srf->th)) {
		fprintf(stderr, "Error reading trace header.\nExiting.\n");
		exit(1);
	    }

	    if( 0 == (level_mode & (LEVEL_CHUNK | LEVEL_BASE)) )
		break;

	    /* Decode ZTR chunks in the header */
	    if (srf->mf)
		mfdestroy(srf->mf);

	    srf->mf = mfcreate(NULL, 0);
	    if (srf->th.trace_hdr_size)
		mfwrite(srf->th.trace_hdr, 1, srf->th.trace_hdr_size, srf->mf);
	    if (srf->ztr)
		delete_ztr(srf->ztr);
	    mrewind(srf->mf);

	    if (NULL != (srf->ztr = partial_decode_ztr(srf, srf->mf, NULL))) {
		srf->mf_pos = mftell(srf->mf);
	    } else {
		/* Maybe not enough to decode or no headerBlob. */
		/* So delay until decoding the body. */
		srf->mf_pos = 0;
	    }
	    mfseek(srf->mf, 0, SEEK_END);
	    srf->mf_end = mftell(srf->mf);

	    break;

	case SRFB_TRACE_BODY: {
	    srf_trace_body_t old_tb;
	    ztr_t *ztr_tmp;
            int no_trace = (level_mode & (LEVEL_CHUNK | LEVEL_BASE) ? 0 : 1);

	    if (0 != srf_read_trace_body(srf, &old_tb, no_trace)) {
		fprintf(stderr, "Error reading trace body.\nExiting.\n");
		exit(1);
	    }

	    if (-1 == construct_trace_name(srf->th.id_prefix,
					   (unsigned char *)old_tb.read_id,
					   old_tb.read_id_length,
					   name, 512)) {
		fprintf(stderr, "Error constructing trace name.\nExiting.\n");
		exit(1);
	    }

	    trace_body_count++;
	    if( 1 == trace_body_count ){
		if( level_mode & LEVEL_NAME )
		    printf( "trace_name: %s + %s", srf->th.id_prefix, name+strlen(srf->th.id_prefix));
	    }
          
	    read_count[READ_TOTAL]++;

	    if (old_tb.flags & SRF_READ_FLAG_BAD_MASK ){
		read_count[READ_BAD]++;
	    } else {
		read_count[READ_GOOD]++;
	    }
          
	    if( 0 == (level_mode & (LEVEL_CHUNK | LEVEL_BASE)) )
		break;

	    if (!srf->mf) {
		fprintf(stderr, "Error reading trace body.\nExiting.\n");
		exit(1);
	    }

	    mfseek(srf->mf, srf->mf_end, SEEK_SET);
	    if (old_tb.trace_size) {
		mfwrite(old_tb.trace, 1, old_tb.trace_size, srf->mf);
		free(old_tb.trace);
		old_tb.trace = NULL;
	    }
          
	    mftruncate(srf->mf, mftell(srf->mf));
	    mfseek(srf->mf, srf->mf_pos, SEEK_SET);

	    if (srf->ztr)
		ztr_tmp = ztr_dup(srf->ztr); /* inefficient, but simple */
	    else
		ztr_tmp = NULL;

	    if ((ztr_tmp = partial_decode_ztr(srf, srf->mf, ztr_tmp))) {
		int i;
		for (i=0; i<ztr_tmp->nchunks; i++) {
		    int ichunk = -1;

		    switch (ztr_tmp->chunk[i].type) {
		    case ZTR_TYPE_BASE:
			ichunk = CHUNK_BASE;
			chunk_size[ichunk] += ztr_tmp->chunk[i].dlength;
			if( parse_base(ztr_tmp, &ztr_tmp->chunk[i], base_count) ){
			    delete_ztr(ztr_tmp);
			    return 1;
			}
			break;
		    case ZTR_TYPE_CNF1:
			ichunk = CHUNK_CNF1;
			chunk_size[ichunk] += ztr_tmp->chunk[i].dlength;
			break;
		    case ZTR_TYPE_CNF4:
			ichunk = CHUNK_CNF4;
			chunk_size[ichunk] += ztr_tmp->chunk[i].dlength;
			break;
		    case ZTR_TYPE_SAMP:
			ichunk = CHUNK_SAMP;
			chunk_size[ichunk] += ztr_tmp->chunk[i].dlength;
			break;
		    case ZTR_TYPE_SMP4:
			ichunk = CHUNK_SMP4;
			chunk_size[ichunk] += ztr_tmp->chunk[i].dlength;
			break;
		    case ZTR_TYPE_REGN:
			ichunk = CHUNK_REGN;
			chunk_size[ichunk] += ztr_tmp->chunk[i].dlength;
			if( NULL == parse_regn(ztr_tmp, &ztr_tmp->chunk[i], regn_hash) ){
			    delete_ztr(ztr_tmp);
			    return 1;
			}
			break;
		    default:
			break;
		    }

		    if( ichunk > -1 ) {
			chunk_count[ichunk]++;
			count_mdata_keys(ztr_tmp, &ztr_tmp->chunk[i], ichunk, key_count, type_count);
		    }
		}

	    }

	    if( ztr_tmp )
		delete_ztr(ztr_tmp);

	    count++;
	    if( (level_mode == LEVEL_CHECK) && (count == 10) ){
		printf( " ... %s x%ld\n", name+strlen(srf->th.id_prefix), trace_body_count);
		srf_destroy(srf, 1);
		return 0;
	    }
          
	    break;
        }

	case SRFB_INDEX: {
            off_t pos = ftell(srf->fp);
	    if( trace_body_count ){
		if( level_mode & LEVEL_NAME )
		    printf( " ... %s x%ld\n", name+strlen(srf->th.id_prefix), trace_body_count);
		trace_body_count = 0;
	    }
            printf( "Reading srf index block\n");
	    if (0 != srf_read_index_hdr(srf, &srf->hdr, 1)) {
		srf_destroy(srf, 1);
		fprintf(stderr, "Error reading srf index block header.\nExiting.\n");
		exit(1);
	    }
            /* Skip the index body */
	    fseeko(srf->fp, pos + srf->hdr.size, SEEK_SET);

            break;
        }
        
	case SRFB_NULL_INDEX: {
            uint64_t ilen;
	    if( trace_body_count ){
		if( level_mode & LEVEL_NAME )
		    printf( " ... %s x%ld\n", name+strlen(srf->th.id_prefix), trace_body_count);
		trace_body_count = 0;
	    }
            printf( "Reading srf null index block\n");
	    /*
	     * Maybe the last 8 bytes of a the file (or previously was
	     * last 8 bytes prior to concatenating SRF files together).
	     * If so it's the index length and should always be 8 zeros.
	     */
            if (1 != fread(&ilen, 8, 1, srf->fp)) {
		srf_destroy(srf, 1);
		fprintf(stderr, "Error reading srf null index block.\nExiting.\n");
		exit(1);
            }
            if (ilen != 0) {
		srf_destroy(srf, 1);
		fprintf(stderr, "Invalid srf null index block.\nExiting.\n");
		exit(1);
            }

            break;
        }
        
        default:
            srf_destroy(srf, 1);
	    fprintf(stderr, "Block of unknown type '%c'\nExiting.\n", type);
	    exit(1);
	}

    }

    if( trace_body_count ){
        if( level_mode & LEVEL_NAME )
            printf( " ... %s x%ld\n", name+strlen(srf->th.id_prefix), trace_body_count);
        trace_body_count = 0;
    }

    /* the type should be -1 (EOF) */
    if( type != -1 ) {
        fprintf(stderr, "Block of unknown type '%c'\nExiting.\n", type);
	exit(1);
    }

    /* are we really at the end of the srf file */
    pos = ftell(srf->fp);
    fseek(srf->fp, 0, SEEK_END);
    if( pos != ftell(srf->fp) ){
        fprintf(stderr, "srf file is corrupt\n");
	exit(1);
    }
    
    srf_destroy(srf, 1);
    return 0;
}

/* ------------------------------------------------------------------------ */

/*
 * Main method.
 */
int main(int argc, char **argv) {
    int ifile, nfiles;
    char *input = NULL;

    int c;
    int errflg = 0;
    extern char *optarg;
    extern int optind, optopt;

    int level_mode = LEVEL_ALL;

    long read_count[NREADS];
    char *read_str[] = {READ_GOOD_STR, READ_BAD_STR, READ_TOTAL_STR};
    long chunk_count[NCHUNKS];
    uint64_t chunk_size[NCHUNKS];
    uint4 chunk_type[] = {CHUNK_BASE_TYPE, CHUNK_CNF1_TYPE, CHUNK_CNF4_TYPE, CHUNK_SAMP_TYPE, CHUNK_SMP4_TYPE, CHUNK_REGN_TYPE};
    long key_count[NCHUNKS][NKEYS];
    char *keys_str[] = {KEY_TYPE_STR, KEY_VALTYPE_STR, KEY_GROUP_STR, KEY_OFFS_STR, KEY_SCALE_STR, KEY_COORD_STR, KEY_NAME_STR};
    long type_count[NCHUNKS][NTYPES];
    char *types_str[] = {TYPE_PROC_STR, TYPE_SLXI_STR, TYPE_SLXN_STR, TYPE_0FAM_STR, TYPE_1CY3_STR, TYPE_2TXR_STR, TYPE_3CY5_STR};
    int iread, ichunk, ikey, itype;

    while ((c = getopt(argc, argv, "l:")) != -1) {
        switch (c) {
        case 'l':
            if (1 != sscanf(optarg, "%d", &level_mode)) {
                fprintf(stderr,
                        "Otion -%c requires an operand\n", optopt);
                errflg++;
            }
	    break;
        case ':':       /* -? without operand */
            fprintf(stderr,
                    "Option -%c requires an operand\n", optopt);
            errflg++;
            break;
        case '?':
            fprintf(stderr,
                    "Unrecognised option: -%c\n", optopt);
            errflg++;
        }
    }

    if (errflg) {
	usage(1);
    }

    nfiles = (argc-optind);
    if( nfiles < 1 ){
        fprintf(stderr, "Please specify input archive name(s).\n");
        usage(1);
    }
    
    for (ifile=0; ifile<nfiles; ifile++) {
        HashTable *regn_hash;
        char bases[] = "ACGTN";
        uint64_t base_count[5];
        char type[5];

        input = argv[optind+ifile];
        printf("Reading archive %s.\n", input);

        for (iread=0; iread<NREADS; iread++)
	    read_count[iread] = 0;

        for (ichunk=0; ichunk<NCHUNKS; ichunk++) {
	    chunk_count[ichunk] = 0;
	    chunk_size[ichunk] = 0;
	}

        for (ichunk=0; ichunk<NCHUNKS; ichunk++)
            for (ikey=0; ikey<NKEYS; ikey++)
                key_count[ichunk][ikey] = 0;

        for (ichunk=0; ichunk<NCHUNKS; ichunk++)
            for (itype=0; itype<NTYPES; itype++)
                type_count[ichunk][itype] = 0;

        if (NULL == (regn_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE|HASH_FUNC_JENKINS3))) {
	    return 1;
        }
    
        memset(base_count, 0, 5 * sizeof(uint64_t));

        if( 0 == srf_info(input, level_mode, read_count,
			  chunk_count, chunk_size,
			  key_count, type_count, regn_hash, base_count) ){

            /* read counts */
            if( level_mode & LEVEL_READ ) {
                for (iread=0; iread<NREADS; iread++) {
                    if( read_count[iread] )
			printf("Reads: %s : %ld\n", read_str[iread], read_count[iread]);
                }
            }

            /* chunk, key and type counts */
            if( level_mode & LEVEL_CHUNK ) {
                for (ichunk=0; ichunk<NCHUNKS; ichunk++) {
                    if( chunk_count[ichunk] ) {
                        printf("Chunk: %s : %ld %"PRId64"\n",
			       ZTR_BE2STR(chunk_type[ichunk], type),
			       chunk_count[ichunk], chunk_size[ichunk]);
                        for (ikey=0; ikey<NKEYS; ikey++) {
                            if(key_count[ichunk][ikey]) {
                                printf("  Mdata key: %s : %ld\n", keys_str[ikey], key_count[ichunk][ikey]);
                                if (ikey == KEY_TYPE && (ichunk == CHUNK_SAMP || ichunk == CHUNK_SMP4)) {
                                    for (itype=0; itype<NTYPES; itype++)
                                        if(type_count[ichunk][itype])
                                            printf("    types: %s : %ld\n", types_str[itype], type_count[ichunk][itype]);
                                }
                                if (ikey == KEY_NAME && (ichunk == CHUNK_REGN)) {
                                    int ibucket;
                                    for (ibucket=0; ibucket<regn_hash->nbuckets; ibucket++) {
                                        HashItem *hi;
                                        for (hi = regn_hash->bucket[ibucket]; hi; hi = hi->next) {
                                            regn_t *regn = (regn_t *)hi->data.p;
                                            printf("    %s x%d\n", hi->key, regn->count);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /* base counts */
            if( level_mode & LEVEL_BASE ) {
                uint64_t total = 0;
                int i;
                for (i=0; i<5; i++) {
                    if( base_count[i] ){
                        printf("Bases: %c: %"PRId64"\n",
			       bases[i], base_count[i]);
                        total += base_count[i];
                    }
                }
                printf("Bases: TOTAL: %"PRId64"\n", total);
            }
        }
    }

    return 0;
}
