/*
 * Copyright (c) 2013 Genome Research Ltd.
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
 * A debugging program to dump out information on the layout of a CRAM file.
 * It's an abomination frankly, but isn't intended for production use.
 */

#include "io_lib_config.h"

#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>

#include <io_lib/cram.h>

void DumpMap2(cram_map **ma, FILE *fp, char *prefix, char *data,
	      HashTable *ds_h) {
    int i, j;
    uintptr_t k;
    for (i = 0; i < CRAM_MAP_HASH; i++) {
	cram_map *m;
	for (m = ma[i]; m; m = m->next) {
	    fprintf(fp, "%s%c%c%c => %16s {",
		    prefix ? prefix : "",
		    (m->key>>16) & 0xff ? (m->key>>16) & 0xff : ' ',
		    (m->key>> 8) & 0xff,
		    (m->key>> 0) & 0xff,
		    cram_encoding2str(m->encoding));
	    for (k = m->offset, j = 0; j < m->size; j++, k++) {
		printf(j ? ", %d" : "%d", (unsigned char)data[k]);
	    }
	    printf("}\n");

	    // Crude, only works with single byte ITF8 values
	    if (m->encoding == E_EXTERNAL ||
		m->encoding == E_BYTE_ARRAY_STOP ||
		m->encoding == E_BYTE_ARRAY_LEN) {
		HashData hd;
		hd.i = (unsigned char)data[m->offset + m->size-1];

		k = (m->key << 8) | hd.i;
		HashTableAdd(ds_h, (char *)k, 4, hd, NULL);
	    }
	}
    }
}

static cram_map *map_find(cram_map **map, unsigned char *key, int id) {
    cram_map *m;

    m = map[CRAM_MAP(key[0],key[1])];
    while (m && m->key != id)
	m= m->next;

    assert(m);

    return m;
}

void dump_core_block(cram_block *b, int verbose) {
    int i;
    int binary = 0, len;

    len = verbose ? b->uncomp_size : MIN(100, b->uncomp_size);
    for (i = 0; i < len; i++) {
	if (isprint(b->data[i]))
	    continue;
	binary++;
    }

    if (binary * 10 > len) {
	printf("Data = {");
	for (i = 0; i < len; i++) {
	    printf(i ? ", %02x" : "%02x", (unsigned char)b->data[i]);
	}

	if (i < b->uncomp_size)
	    puts(", ...}");
	else
	    puts("}");

    } else {
	for (i = 0; i < len; i++) {
	    if (isprint(b->data[i]))
		putchar(b->data[i]);
	    else
		printf("\\%03o", b->data[i]);
	}

	if (i < b->uncomp_size)
	    puts("...");
    }
}

void dump_seq_block(cram_block *b, int verbose) {
    int i;
    for (i = 0; (verbose || i < 100) && i < b->uncomp_size; i++) {
	if (isprint(b->data[i]))
	    putchar(b->data[i]);
	else
	    printf("\\%03o", b->data[i]);
    }
}

void dump_quality_block(cram_block *b, int verbose) {
    int i;
    for (i = 0; (verbose || i < 100) && i < b->uncomp_size; i++) {
	if (isprint(b->data[i] + '!'))
	    putchar(b->data[i] + '!');
	else
	    printf("\\%03o", b->data[i] + '!');
    }
    putchar('\n');
}

void dump_name_block(cram_block *b, int verbose) {
    int i;
    for (i = 0; (verbose || i < 100) && i < b->uncomp_size; i++) {
	if (isprint(b->data[i]))
	    putchar(b->data[i]);
	else
	    printf("\\%03o", b->data[i]);
    }
}

void dump_tag_block(cram_block *b, int verbose) {
    return dump_core_block(b, verbose);
}

int main(int argc, char **argv) {
    cram_fd *fd;
    cram_container *c;
    off_t pos, pos2, hpos;
    int verbose = 0;
    HashTable *bsize_h;
    HashTable *ds_h; // content_id to data-series lookup.
    HashTable *dc_h; // content_id to data-compression lookup

    static int bmax = 0;
    bsize_h = HashTableCreate(128, HASH_DYNAMIC_SIZE|
			    HASH_NONVOLATILE_KEYS |
			    HASH_INT_KEYS);
    ds_h = HashTableCreate(128, HASH_DYNAMIC_SIZE|
			   HASH_NONVOLATILE_KEYS |
			   HASH_INT_KEYS);
    dc_h = HashTableCreate(128, HASH_DYNAMIC_SIZE);

    if (argc >= 2 && strcmp(argv[1], "-v") == 0) {
	argc--;
	argv++;
	verbose = 1;
    }

    if (argc < 2) {
	fprintf(stderr, "Usage: cram_dump [-v] filename.cram\n");
	return 1;
    }

    if (NULL == (fd = cram_open(argv[1], "rb"))) {
	fprintf(stderr, "Error opening CRAM file '%s'.\n", argv[1]);
	return 1;
    }

    printf("File definition structure\n");
    printf("    Major vers: %d\n", fd->file_def->major_version);
    printf("    Minor vers: %d\n", fd->file_def->minor_version);
    printf("    file_id:    %.20s\n", fd->file_def->file_id);

    printf("\nBAM header:\n%.*s\n",
	   sam_hdr_length(fd->header),
	   sam_hdr_str(fd->header));

    pos = CRAM_IO_TELLO(fd);
    while ((c = cram_read_container(fd))) {
	int i, j;

	if (fd->err) {
	    perror("Cram container read");
	    return 1;
	}

	printf("\nContainer pos %"PRId64" size %d\n", (int64_t)pos, c->length);
	printf("    Ref id:            %d\n", c->ref_seq_id);
	printf("    Ref pos:           %d + %d\n", c->ref_seq_start, c->ref_seq_span);
	printf("    Rec counter:       %"PRId64"\n", c->record_counter);
       	printf("    No. recs:          %d\n", c->num_records);
	printf("    No. bases          %"PRId64"\n", c->num_bases);
	printf("    No. blocks:        %d\n", c->num_blocks);
	printf("    No. landmarks:     %d\n", c->num_landmarks);

	printf("    {");
	for (i = 0; i < c->num_landmarks; i++) {
	    printf(i ? ", %d" : "%d", c->landmark[i]);
	}
	printf("}\n");

	hpos = CRAM_IO_TELLO(fd);

	if (!c->length) {
	    //printf("\n    EMPTY BLOCK\n");
	    pos = CRAM_IO_TELLO(fd);
	    continue;
	}
	printf("\n    Container_header block pos %"PRId64"\n", (int64_t)hpos);
	if (!(c->comp_hdr_block = cram_read_block(fd)))
	    return 1;
	assert(c->comp_hdr_block->content_type == COMPRESSION_HEADER);

	c->comp_hdr = cram_decode_compression_header(fd, c->comp_hdr_block);
	if (!c->comp_hdr)
	    return 1;

	printf("      Preservation map:\n");
	HashTableDump(c->comp_hdr->preservation_map, stdout, "\t");
	printf("      Substitution map:\n");
	printf("        A: %.4s\n", c->comp_hdr->substitution_matrix[0]);
	printf("        C: %.4s\n", c->comp_hdr->substitution_matrix[1]);
	printf("        G: %.4s\n", c->comp_hdr->substitution_matrix[2]);
	printf("        T: %.4s\n", c->comp_hdr->substitution_matrix[3]);
	printf("        N: %.4s\n", c->comp_hdr->substitution_matrix[4]);

	printf("      TD map:\n");
	for (i = 0; i < c->comp_hdr->nTL; i++)
	    printf("        %3d: %s\n", i, c->comp_hdr->TL[i]);

	printf("\n      Record encoding map:\n");
	DumpMap2(c->comp_hdr->rec_encoding_map, stdout, "\t", 
		 (char *)c->comp_hdr_block->data, ds_h);

	printf("\n      Tag encoding map:\n");
	DumpMap2(c->comp_hdr->tag_encoding_map, stdout, "\t",
		 (char *)c->comp_hdr_block->data, ds_h);

	for (j = 0; j < c->num_landmarks; j++) {
	    cram_slice *s;
	    int id;
	    
	    pos2 = CRAM_IO_TELLO(fd);
	    assert(pos2 - pos - c->offset == c->landmark[j]);

	    s = cram_read_slice(fd);
	    printf("\n    Slice %d/%d, container offset %d, file offset %d\n", j+1, c->num_landmarks, (int)(pos2 - pos - c->offset), (int)pos2);
	    printf("\tSlice content type %s\n",
		   cram_content_type2str(s->hdr->content_type));

	    if (s->hdr->content_type == MAPPED_SLICE) {
		int i;
		printf("\tSlice ref seq    %d\n", s->hdr->ref_seq_id);
		printf("\tSlice ref start  %d\n", s->hdr->ref_seq_start);
		printf("\tSlice ref span   %d\n", s->hdr->ref_seq_span);
		printf("\tSlice MD5        ");
		for (i = 0; i < 16; i++)
		    printf("%02x", s->hdr->md5[i]);
		putchar('\n');
	    }
	    printf("\tRec counter      %"PRId64"\n", s->hdr->record_counter);
	    printf("\tNo. records      %d\n", s->hdr->num_records);
	    printf("\tNo. blocks       %d\n", s->hdr->num_blocks);
	    printf("\tBlk IDS:         {");
	    for (id = 0; id < s->hdr->num_content_ids; id++) {
		printf(id ? ", %d" : "%d", s->hdr->block_content_ids[id]);
	    }
	    printf("}\n");
	    if (s->hdr->content_type == MAPPED_SLICE) {
		printf("\tRef base id:     %d\n", s->hdr->ref_base_id);
	    }

	    if (s->hdr->tags) {
		HashIter *iter;
		HashItem *hi;

		iter = HashTableIterCreate();
		while ((hi = HashTableIterNext(s->hdr->tags, iter))) {
		    printf("\tOptional tag %c%c:%c:",
			   hi->key[0], hi->key[1], hi->key[2]);

		    switch(hi->key[2]) {
			uint32_t len;
			unsigned char *dat;
		    case 'i':
			printf("%"PRId64"\n", hi->data.i);
			break;
		    case 'f':
			printf("%f\n", hi->data.f);
			break;
		    case 'Z': case 'H':
			printf("%s\n", (char *)hi->data.p);
			break;
		    case 'A':
			printf("<%d>\n", (unsigned char)hi->data.i);
			break;
		    case 'B':
			dat = hi->data.p;
			len = dat[1] | (dat[2]<<8) | (dat[3]<<16)| (dat[4]<<24);
			switch(dat[0]) {
			case 's': case 'S':
			    len *= 2;
			    break;
			case 'i': case 'I': case 'f':
			    len *= 4;
			    break;
			default:
			    break;
			}
			putchar(dat[0]);
			dat += 5;
			while (len--) {
			    printf(",%02x", *dat++);
			}
			putchar('\n');
		    }
		}
	    }
	
	    for (id = 0; id < s->hdr->num_blocks; id++) {
		HashItem *hi;
		intptr_t k = s->block[id]->content_type == CORE
		    ? -1 : s->block[id]->content_id;
		hi = HashTableSearch(bsize_h, (char *)k, sizeof(k));
		if (hi) {
		    hi->data.i += s->block[id]->comp_size;
		} else {
		    HashData hd;
		    hd.i = s->block[id]->comp_size;
		    HashTableAdd(bsize_h, (char *)k, sizeof(k), hd, NULL);
		}

		// WARNING: scuppered by having high content_id values.
		if (bmax < s->block[id]->content_id)
		    bmax = s->block[id]->content_id;
	    }

	    for (id = 0; id < s->hdr->num_blocks; id++)
		cram_uncompress_block(s->block[id]);

	    /* Test decoding of 1st seq */
	    if (verbose) {
		cram_block *b = s->block[0];
		int32_t i32, bf, fn, prev_pos, rl;
		unsigned char cf;
		int out_sz, r, f;
		int rec;

		assert(b->content_type == CORE);

		for (rec = 0; rec < s->hdr->num_records; rec++) {
		    unsigned char ntags;

		    printf("Rec %d/%d at %d,%d\n", rec+1, s->hdr->num_records,
			   (int)b->byte, b->bit);

		    out_sz = 1; /* decode 1 item */
		    r = c->comp_hdr->codecs[DS_BF]->decode(s,c->comp_hdr->codecs[DS_BF], b, (char *)&bf, &out_sz);
		    printf("BF = %d => SAM 0x%x (ret %d, out_sz %d)\n", bf, fd->bam_flag_swap[bf], r, out_sz);
		    bf = fd->bam_flag_swap[bf];

		    if (!IS_CRAM_1_VERS(fd)) {
			r = c->comp_hdr->codecs[DS_CF]->decode(s,c->comp_hdr->codecs[DS_CF], b, (char *)&i32, &out_sz);
			printf("CF = %d (ret %d, out_sz %d)\n", i32, r, out_sz);
			cf = i32;
		    } else {
			r = c->comp_hdr->codecs[DS_CF]->decode(s,c->comp_hdr->codecs[DS_CF], b, (char *)&cf, &out_sz);
			printf("CF = %d (ret %d, out_sz %d)\n", cf, r, out_sz);
		    }

		    if (!IS_CRAM_1_VERS(fd) && s->hdr->ref_seq_id == -2) {
			int32_t ri;
			r |= c->comp_hdr->codecs[DS_RI]->decode(s, c->comp_hdr->codecs[DS_RI], b, (char *)&ri, &out_sz);
			printf("RI = %d (ret %d, out_sz %d)\n", ri, r, out_sz);
		    }

		    r = c->comp_hdr->codecs[DS_RL]->decode(s,c->comp_hdr->codecs[DS_RL], b, (char *)&rl, &out_sz);
		    printf("RL = %d (ret %d, out_sz %d)\n", rl, r, out_sz);

		    r = c->comp_hdr->codecs[DS_AP]->decode(s,c->comp_hdr->codecs[DS_AP], b, (char *)&i32, &out_sz);
		    printf("AP = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

		    r = c->comp_hdr->codecs[DS_RG]->decode(s,c->comp_hdr->codecs[DS_RG], b, (char *)&i32, &out_sz);
		    printf("RG = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

		    if (c->comp_hdr->read_names_included) {
			int32_t out_sz2 = 1;
			cram_block *dat = cram_new_block(EXTERNAL, 0);

			r = c->comp_hdr->codecs[DS_RN]->decode(s,c->comp_hdr->codecs[DS_RN], b, (char *)dat, &out_sz2);
			printf("RN = %.*s (ret %d, out_sz %d)\n", out_sz2, BLOCK_DATA(dat), r, out_sz2);
			cram_free_block(dat);
		    }

		    if (cf & CRAM_FLAG_DETACHED) {
			char mf;
			puts("Detached");
			/* MF, RN if !captureReadNames, NS, NP, IS */
			if (!IS_CRAM_1_VERS(fd)) {
			    r = c->comp_hdr->codecs[DS_MF]->decode(s,c->comp_hdr->codecs[DS_MF], b, (char *)&i32, &out_sz);
			    printf("MF = %d (ret %d, out_sz %d)\n", i32, r, out_sz);
			} else {
			    r = c->comp_hdr->codecs[DS_MF]->decode(s,c->comp_hdr->codecs[DS_MF], b, &mf, &out_sz);
			    printf("MF = %d (ret %d, out_sz %d)\n", mf, r, out_sz);
			}
			if (!c->comp_hdr->read_names_included) {
			    cram_block *dat = cram_new_block(EXTERNAL, 0);
			    int32_t out_sz2 = 1;

			    r = c->comp_hdr->codecs[DS_RN]->decode(s,c->comp_hdr->codecs[DS_RN], b, (char *)dat, &out_sz2);
			    printf("RN = %.*s (ret %d, out_sz %d)\n", out_sz2, BLOCK_DATA(dat), r, out_sz2);
			    cram_free_block(dat);
			}
		    
			r = c->comp_hdr->codecs[DS_NS]->decode(s,c->comp_hdr->codecs[DS_NS], b, (char *)&i32, &out_sz);
			printf("NS = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

			r = c->comp_hdr->codecs[DS_NP]->decode(s,c->comp_hdr->codecs[DS_NP], b, (char *)&i32, &out_sz);
			printf("NP = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

			r = c->comp_hdr->codecs[DS_TS]->decode(s,c->comp_hdr->codecs[DS_TS], b, (char *)&i32, &out_sz);
			printf("TS = %d (ret %d, out_sz %d)\n", i32, r, out_sz);
		    } else if (cf & CRAM_FLAG_MATE_DOWNSTREAM) {
			puts("Not detached, and mate is downstream");
			r = c->comp_hdr->codecs[DS_NF]->decode(s,c->comp_hdr->codecs[DS_NF], b, (char *)&i32, &out_sz);
			printf("NF = %d+%d = %d (ret %d, out_sz %d)\n", i32, rec+1, i32+rec+1, r, out_sz);
		    }

		    if (IS_CRAM_1_VERS(fd)) {
			r = c->comp_hdr->codecs[DS_TC]->decode(s,c->comp_hdr->codecs[DS_TC], b, (char *)&ntags, &out_sz);
			printf("TC = %d (ret %d, out_sz %d)\n", ntags, r, out_sz);

			for (f = 0; f < ntags; f++) {
			    int32_t id;
			    char key[3];
			    cram_map *m;
			    cram_block *tag = cram_new_block(EXTERNAL, 0);

			    r = c->comp_hdr->codecs[DS_TN]->decode(s, c->comp_hdr->codecs[DS_TN],
							      b, (char *)&id, &out_sz);
			    key[0] = (id>>16)&0xff;
			    key[1] = (id>>8)&0xff;
			    key[2] = id&0xff;
			    printf("%3d: TN= %.3s\n", f, key);

			    printf("id=%d\n", id);
			    if ((m = map_find(c->comp_hdr->tag_encoding_map,
					      (unsigned char *)key, id))) {
				int i, out_sz;

				BLOCK_SIZE(tag) = 0;
				r = m->codec->decode(s, m->codec, b, (char *)tag, &out_sz);
				printf("%3d: Val", f);
				for(i = 0; i < out_sz; i++) {
				    printf(" %02x", BLOCK_DATA(tag)[i]);
				}
				printf("\n");
			    } else {
				fprintf(stderr, "*** ERROR: unrecognised aux key ***\n");
			    }

			    cram_free_block(tag);
			    // skip decoding of tag data itself and hope it's
			    // in an external block.
			}
		    } else {
			int32_t tl, ntags;
			char *tn;
			r = c->comp_hdr->codecs[DS_TL]->decode(s,c->comp_hdr->codecs[DS_TL], b, (char *)&tl, &out_sz);
			printf("TL = %d (ret %d, out_sz %d)\n", tl, r, out_sz);

			tn = (char *)c->comp_hdr->TL[tl];
			ntags = strlen(tn)/3;

			for (f = 0; f < ntags; f++) {
			    int32_t id;
			    char key[3];
			    cram_map *m;
			    cram_block *tag = cram_new_block(EXTERNAL, 0);

			    key[0] = *tn++;
			    key[1] = *tn++;
			    key[2] = *tn++;
			    id = (key[0]<<16) | (key[1]<<8) | key[2];
			    printf("%3d: TN= %.3s\n", f, key);

			    printf("id=%d\n", id);
			    if ((m = map_find(c->comp_hdr->tag_encoding_map,
					      (unsigned char *)key, id))) {
				int i, out_sz;

				BLOCK_SIZE(tag) = 0;
				r = m->codec->decode(s, m->codec, b, (char *)tag, &out_sz);
				printf("%3d: Val", f);
				for(i = 0; i < out_sz; i++) {
				    printf(" %02x", BLOCK_DATA(tag)[i]);
				}
				printf("\n");
			    } else {
				fprintf(stderr, "*** ERROR: unrecognised aux key ***\n");
			    }

			    cram_free_block(tag);
			    // skip decoding of tag data itself and hope it's
			    // in an external block.
			}
		    }

		    if (!(bf & BAM_FUNMAP)) {
			r = c->comp_hdr->codecs[DS_FN]->decode(s,c->comp_hdr->codecs[DS_FN], b, (char *)&fn, &out_sz);
			printf("FN = %d (ret %d, out_sz %d)\n", fn, r, out_sz);

			prev_pos = 0;
			for (f = 0; f < fn; f++) {
			    char op;
			    int32_t pos;

			    r = c->comp_hdr->codecs[DS_FC]->decode(s,c->comp_hdr->codecs[DS_FC], b, &op, &out_sz);
			    printf("  %d: FC = %c (ret %d, out_sz %d)\n", f, op, r, out_sz);

			    r = c->comp_hdr->codecs[DS_FP]->decode(s,c->comp_hdr->codecs[DS_FP], b, (char *)&pos, &out_sz);
			    printf("  %d: FP = %d+%d = %d (ret %d, out_sz %d)\n", f, pos, prev_pos, pos + prev_pos, r, out_sz);

			    pos += prev_pos;
			    prev_pos = pos;

			    switch(op) {
			    case 'S': { // soft clip: IN
				char dat[100];
				int32_t out_sz2 = 1;

				dat[0]='?';dat[1]=0;
				if (IS_CRAM_1_VERS(fd)) {
				    if (c->comp_hdr->codecs[DS_IN]) {
					r = c->comp_hdr->codecs[DS_IN]->decode(s,c->comp_hdr->codecs[DS_IN], b, dat, &out_sz2);
					printf("  %d: IN(S) = %.*s (ret %d, out_sz %d)\n", f, out_sz2, dat, r, out_sz2);
				    }
				} else {
				    if (c->comp_hdr->codecs[DS_SC]) {
					r = c->comp_hdr->codecs[DS_SC]->decode(s,c->comp_hdr->codecs[DS_SC], b, dat, &out_sz2);
					printf("  %d: SC(S) = %.*s (ret %d, out_sz %d)\n", f, out_sz2, dat, r, out_sz2);
				    }
				}
				break;
			    }

			    case 'X': { // Substitution; BS
				char bs;
				r = c->comp_hdr->codecs[DS_BS]->decode(s,c->comp_hdr->codecs[DS_BS], b, &bs, &out_sz);
				printf("  %d: BS = %d (ret %d)\n", f, bs, r);
				break;
			    }

			    case 'D': { // Deletion; DL
				r = c->comp_hdr->codecs[DS_DL]->decode(s,c->comp_hdr->codecs[DS_DL], b, (char *)&i32, &out_sz);
				printf("  %d: DL = %d (ret %d)\n", f, i32, r);
				break;
			    }

			    case 'I': { // Insertion (several bases); IN
				static char *dat = NULL;
				static int dat_l = 0;
				int32_t out_sz2 = 1;

				if (dat_l < rl+1) {
				    dat = realloc(dat, rl+1);
				    dat_l = rl;
				}

				dat[0]='?';dat[1]=0;
				r = c->comp_hdr->codecs[DS_IN]->decode(s,c->comp_hdr->codecs[DS_IN], b, dat, &out_sz2);
				printf("  %d: IN(I) = %.*s (ret %d, out_sz %d)\n", f, out_sz2, dat, r, out_sz2);
				break;
			    }

			    case 'i': { // Insertion (single base); BA
				char cc;
				r = c->comp_hdr->codecs[DS_BA]->decode(s,c->comp_hdr->codecs[DS_BA], b, &cc, &out_sz);
				printf("  %d: BA = %c (ret %d)\n", f, cc, r);
				break;
			    }

			    case 'b': { // Read bases; BB
				static char *seq = NULL;
				static int seq_l = 0;
				int out_sz2;

				if (seq_l < rl) {
				    seq = realloc(seq, rl);
				    seq_l = rl;
				}

				r  = c->comp_hdr->codecs[DS_BB]->decode(s,c->comp_hdr->codecs[DS_BB], b, seq, &out_sz2);
				printf("  %d: BB(b) = %.*s (ret %d, out_sz %d)\n", f, out_sz2, seq, r, out_sz2);
				break;
			    }

			    case 'q': { // Read bases; QQ
				int out_sz2;
				static char *qual = NULL;
				static int qual_l = 0;

				if (qual_l < rl) {
				    qual = realloc(qual, rl);
				    qual_l = rl;
				}

				r  = c->comp_hdr->codecs[DS_QQ]->decode(s,c->comp_hdr->codecs[DS_QQ], b, qual, &out_sz2);
				printf("  %d: QQ(b) = %.*s (ret %d, out_sz %d)\n", f, out_sz2, qual, r, out_sz2);
				break;
			    }

			    case 'B': { // Read base; BA, QS
				char cc, qc;
				r  = c->comp_hdr->codecs[DS_BA]->decode(s,c->comp_hdr->codecs[DS_BA], b, &cc, &out_sz);
				r |= c->comp_hdr->codecs[DS_QS]->decode(s,c->comp_hdr->codecs[DS_QS], b, &qc, &out_sz);
				printf("  %d: BA/QS(B) = %c/%d (ret %d)\n", f, cc, qc, r);
				break;
			    }

			    case 'Q': { // Quality score; QS
				char qc;
				r = c->comp_hdr->codecs[DS_QS]->decode(s,c->comp_hdr->codecs[DS_QS], b, &qc, &out_sz);
				printf("  %d: QS = %d (ret %d)\n", f, qc, r);
				break;
			    }

			    case 'N':
				r = c->comp_hdr->codecs[DS_RS]->decode(s,c->comp_hdr->codecs[DS_RS], b, (char *)&i32, &out_sz);
				printf("  %d: RS = %d (ret %d)\n", f, i32, r);
				break;

			    case 'P':
				r = c->comp_hdr->codecs[DS_PD]->decode(s,c->comp_hdr->codecs[DS_PD], b, (char *)&i32, &out_sz);
				printf("  %d: PD = %d (ret %d)\n", f, i32, r);
				break;

			    case 'H':
				r = c->comp_hdr->codecs[DS_HC]->decode(s,c->comp_hdr->codecs[DS_HC], b, (char *)&i32, &out_sz);
				printf("  %d: HC = %d (ret %d)\n", f, i32, r);
				break;

			    default:
				abort();
			    }
			}

			r = c->comp_hdr->codecs[DS_MQ]->decode(s,c->comp_hdr->codecs[DS_MQ], b, (char *)&i32, &out_sz);
			printf("MQ = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

			if (cf & CRAM_FLAG_PRESERVE_QUAL_SCORES) {
			    char dat[1024];
			    int len = rl;

			    do {
				int32_t out_sz2 = len > 1024 ? 1024 : len, i;
				dat[0]='?';dat[1]=0;
				r = c->comp_hdr->codecs[DS_QS]->decode(s,c->comp_hdr->codecs[DS_QS], b, dat, &out_sz2);
				for (i = 0; i < out_sz2; i++)
				    dat[i] += '!';
				printf("QS = %.*s (ret %d, out_sz %d)\n", out_sz2, dat, r, out_sz2);
				len -= 1024;
			    } while (len > 0);
			}
		    } else {
			puts("Unmapped");
			char dat[1024];
			int len = rl;

			while (len > 0) {
			    int32_t out_sz2 = len > 1024 ? 1024 : len;
			    r = c->comp_hdr->codecs[DS_BA]->decode(s, c->comp_hdr->codecs[DS_BA], b, dat, &out_sz2);
			    printf("SQ = %.*s (out_sz %d)\n", out_sz2, dat, out_sz2);
			    len -= 1024;
			}

			if (cf & CRAM_FLAG_PRESERVE_QUAL_SCORES) {
			    int len = rl, i;

			    do {
				int32_t out_sz2 = len > 1024 ? 1024 : len;
				r = c->comp_hdr->codecs[DS_QS]->decode(s, c->comp_hdr->codecs[DS_QS], b, dat, &out_sz2);
				for (i = 0; i < out_sz2; i++)
				    dat[i] += '!';
				printf("QS = %.*s (out_sz %d)\n", out_sz2, dat, out_sz2);
				len -= 1024;
			    } while (len > 0);
			}
		    }
		}
	    }

	    for (id = 0; id < s->hdr->num_blocks; id++) {
		HashData hd = {0};
		cram_block *b = s->block[id];
		printf("\n\tBlock %d/%d\n", id+1, s->hdr->num_blocks);
		printf("\t    Size:         %d comp / %d uncomp\n",
		       b->comp_size, b->uncomp_size);
		printf("\t    Method:       %s\n",
		       cram_block_method2str(b->orig_method));
		struct {
		    int id;
		    enum cram_block_method method;
		} id_type = {b->content_id, b->orig_method};
		HashTableAdd(dc_h, (char *)&id_type, sizeof(id_type), hd, NULL);
		printf("\t    Content type: %s\n",
		       cram_content_type2str(b->content_type));
		printf("\t    Content id:   %d\n", b->content_id);

		if (b->method != RAW)
		    cram_uncompress_block(b);

		if (b->content_type == CORE) {
		    dump_core_block(b, verbose);
		} else {
		    int t, m;
		    enum cram_fields cf = 0;
		    char fields[1024], *fp = fields;
		    for (m = 0; m < 2; m++) {
			cram_map **ma = m
			    ? c->comp_hdr->tag_encoding_map
			    : c->comp_hdr->rec_encoding_map;
			
			for (t = 0; t < CRAM_MAP_HASH; t++) {
			    cram_map *m;
			    unsigned char *data = c->comp_hdr_block->data;
			    for (m = ma[t]; m; m = m->next) {
				if (m->encoding != E_EXTERNAL &&
				    m->encoding != E_BYTE_ARRAY_STOP &&
				    m->encoding != E_BYTE_ARRAY_LEN)
				    continue;
				if (data[m->offset + m->size-1] !=
				    b->content_id)
				    continue;
				if ((m->key>>16)&0xff)
				    *fp++ = (m->key>>16)&0xff;
				*fp++ = (m->key>> 8)&0xff;
				*fp++ = (m->key>> 0)&0xff;
				*fp++ = ' ';

				if (m->key>>16)             cf |= CRAM_aux;
				if (m->key == ('Q'<<8)+'S') cf |= CRAM_QS;
				if (m->key == ('R'<<8)+'N') cf |= CRAM_RN;
				if (m->key == ('S'<<8)+'C') cf |= CRAM_SC;
				if (m->key == ('I'<<8)+'N') cf |= CRAM_IN;
			    }
			}
		    }
		    *fp++ = 0;
		    printf("\t    Keys:         %s\n", fields);

		    if (cf == CRAM_aux)
			dump_tag_block(b, verbose);
		    else if (cf == CRAM_QS)
			dump_quality_block(b, verbose);
		    else if (cf == CRAM_SC)
			dump_seq_block(b, verbose);
		    else if (cf == CRAM_RN)
			dump_name_block(b, verbose);
		    else
			dump_core_block(b, verbose);
		}
	    }

	    cram_free_slice(s);
	}

	pos = CRAM_IO_TELLO(fd);
	assert(pos == hpos + c->length);

	cram_free_container(c);
    }

    cram_close(fd);

    {
	intptr_t id;

	puts("");
	for (id = -1; id <= bmax; id++) {
	    intptr_t k;
	    HashItem *hi;
	    HashIter *iter;

	    if (!(hi = HashTableSearch(bsize_h, (char *)id, sizeof(id))))
		continue;

	    k = (intptr_t) hi->key;
	    if (k == -1) {
		printf("Block CORE              , total size %11"PRId64"\n", hi->data.i);
		continue;
	    }

	    printf("Block content_id %7d, total size %11"PRId64" ",
		   (int) k, hi->data.i);

	    struct {
		int id;
		enum cram_block_method method;
	    } id_type = {k, 0};

	    id_type.method = GZIP;
	    putchar((HashTableSearch(dc_h, (char *)&id_type, sizeof(id_type)))?'g':' ');

	    id_type.method = BZIP2;
	    putchar((HashTableSearch(dc_h, (char *)&id_type, sizeof(id_type)))?'b':' ');

	    id_type.method = LZMA;
	    putchar((HashTableSearch(dc_h, (char *)&id_type, sizeof(id_type)))?'l':' ');

	    id_type.method = RANS0;
	    putchar((HashTableSearch(dc_h, (char *)&id_type, sizeof(id_type)))?'r':' ');

	    id_type.method = RANS1;
	    putchar((HashTableSearch(dc_h, (char *)&id_type, sizeof(id_type)))?'R':' ');

	    iter = HashTableIterCreate();
	    while ((hi = HashTableIterNext(ds_h, iter))) {
		int c;
		char buf[5];
		int x = 4;

		if (hi->data.i != k)
		    continue;
		
		c = ((uintptr_t) hi->key)>>8;
		
		buf[x--] = 0;
		while(c & 0xff) {
		    buf[x--] = c;
		    c >>= 8;
		}
		printf(" %s", &buf[x+1]);
	    }
	    putchar('\n');

	    HashTableIterDestroy(iter);
	}
    }

    HashTableDestroy(bsize_h, 0);
    HashTableDestroy(ds_h, 0);
    HashTableDestroy(dc_h, 0);

    return 0;
}
