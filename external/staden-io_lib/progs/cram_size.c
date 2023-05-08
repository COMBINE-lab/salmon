/*
 * Copyright (c) 2013-2015 Genome Research Ltd.
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
 * A cut down version of cram_dump.c to accumulate only size
 * information per data series.  This is much faster than cram_dump as
 * it does not require uncompression of data blocks.
 */

#include "io_lib_config.h"

#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>

#include <io_lib/cram.h>

// Accumulate per {data_series, content_id} combination.
void ParseMap(cram_map **ma, char *data, HashTable *ds_h) {
    int i;
    uintptr_t k;
    for (i = 0; i < CRAM_MAP_HASH; i++) {
	cram_map *m;
	for (m = ma[i]; m; m = m->next) {
	    // Crude, only works with single byte ITF8 values
	    if (m->encoding == E_EXTERNAL ||
		m->encoding == E_BYTE_ARRAY_STOP ||
		m->encoding == E_BYTE_ARRAY_LEN) {
		HashData hd;
		hd.i = (unsigned char)data[m->offset + m->size-1];

		k = (m->key << 8) | hd.i;

		cram_codec *c = cram_decoder_init(m->encoding,
						  data + m->offset,
						  m->size, E_BYTE_ARRAY, 0);
		int id1, id2;
		if (c) {
		    id1 = cram_codec_to_id(c, &id2);
		    c->free(c);
		    if (id1 >= 0) {
			hd.i = id1;
			HashTableAdd(ds_h, (char *)k, 4, hd, NULL);
		    }
		    if (id2 >= 0) {
			hd.i = id2;
			HashTableAdd(ds_h, (char *)k, 4, hd, NULL);
		    }
		} else {
		    HashTableAdd(ds_h, (char *)k, 4, hd, NULL);
		}
	    }
	}
    }
}

void print_sizes(HashTable *bsize_h, HashTable *ds_h, HashTable *dc_h, int bmax) {
    intptr_t id;

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

	enum cram_block_method methods[] = {GZIP, BZIP2, LZMA, RANS0, RANS1};
	int count[1+sizeof(methods)/sizeof(*methods)] = {0};

	int i;
	for (i = 0; i < sizeof(methods)/sizeof(*methods); i++) {
	    id_type.method = methods[i];
	    if ((hi = HashTableSearch(dc_h, (char *)&id_type, sizeof(id_type)))) {
		count[0] += hi->data.i;
		count[i+1] += hi->data.i;
	    }
	}

	for (i = 1; i <= sizeof(methods)/sizeof(*methods); i++) {
	    printf("%s%s%c%s",
		   (count[i]+0.01)/(count[0]+0.01) >= 0.50 ? "\033[7m":"",
		   (count[i]+0.01)/(count[0]+0.01) >= 0.10 ? "\033[4m":"",
		   " gblrR"[count[i]?i:0],
		   "\033[0m");
	}

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

int process_sizes(cram_fd *fd,
		  HashTable *bsize_h, HashTable *ds_h, HashTable *dc_h,
		  int *bmax) {
    cram_container *c;
    off_t pos, pos2, hpos;

    pos = CRAM_IO_TELLO(fd);
    while ((c = cram_read_container(fd))) {
	int j;

	if (fd->err) {
	    perror("Cram container read");
	    return 1;
	}

	hpos = CRAM_IO_TELLO(fd);

	if (!c->length) {
	    pos = CRAM_IO_TELLO(fd);
	    continue;
	}
	if (!(c->comp_hdr_block = cram_read_block(fd)))
	    return 1;
	assert(c->comp_hdr_block->content_type == COMPRESSION_HEADER);

	c->comp_hdr = cram_decode_compression_header(fd, c->comp_hdr_block);
	if (!c->comp_hdr)
	    return 1;

	ParseMap(c->comp_hdr->rec_encoding_map,
		 (char *)c->comp_hdr_block->data, ds_h);
	ParseMap(c->comp_hdr->tag_encoding_map,
		 (char *)c->comp_hdr_block->data, ds_h);

	for (j = 0; j < c->num_landmarks; j++) {
	    cram_slice *s;
	    int id;
	    
	    pos2 = CRAM_IO_TELLO(fd);
	    assert(pos2 - pos - c->offset == c->landmark[j]);

	    s = cram_read_slice(fd);
	    
	    // Accumulate per block content_id
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
		if (*bmax < s->block[id]->content_id)
		    *bmax = s->block[id]->content_id;
	    }

	    // Hack for rans1
	    for (id = 0; id < s->hdr->num_blocks; id++) {
		if (s->block[id]->comp_size >= 2 &&
		    s->block[id]->orig_method == RANS0 &&
		    s->block[id]->data[0] != 0)
		    s->block[id]->orig_method = RANS1;
	    }

	    // Accumulate per {id, compression_method} combo
	    for (id = 0; id < s->hdr->num_blocks; id++) {
		HashData hd;
		HashItem *hi;
		cram_block *b = s->block[id];
		struct {
		    int id;
		    enum cram_block_method method;
		} id_type = {b->content_id, b->orig_method};
		hd.i = 0;
		hi = HashTableAdd(dc_h, (char *)&id_type, sizeof(id_type), hd, NULL);
		hi->data.i++;

		int t, m;
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
			    if (data[m->offset + m->size-1] != b->content_id)
				continue;
			    if ((m->key>>16)&0xff)
				*fp++ = (m->key>>16)&0xff;
			    *fp++ = (m->key>> 8)&0xff;
			    *fp++ = (m->key>> 0)&0xff;
			    *fp++ = ' ';
			}
		    }
		}
		*fp++ = 0;
	    }

	    cram_free_slice(s);
	}

	pos = CRAM_IO_TELLO(fd);
	assert(pos == hpos + c->length);

	cram_free_container(c);
    }

    return 0;
}

int main(int argc, char **argv) {
    cram_fd *fd;
    HashTable *bsize_h;
    HashTable *ds_h; // content_id to data-series lookup.
    HashTable *dc_h; // content_id to data-compression lookup
    int bmax = 0;

    bsize_h = HashTableCreate(128, HASH_DYNAMIC_SIZE|
			    HASH_NONVOLATILE_KEYS |
			    HASH_INT_KEYS);
    ds_h = HashTableCreate(128, HASH_DYNAMIC_SIZE|
			   HASH_NONVOLATILE_KEYS |
			   HASH_INT_KEYS);
    dc_h = HashTableCreate(128, HASH_DYNAMIC_SIZE);

    if (argc < 2) {
	fprintf(stderr, "Usage: cram_size filename.cram\n");
	return 1;
    }

    if (NULL == (fd = cram_open(argv[1], "rb"))) {
	fprintf(stderr, "Error opening CRAM file '%s'.\n", argv[1]);
	return 1;
    }

    if (0 != process_sizes(fd, bsize_h, ds_h, dc_h, &bmax))
	return 1;
    cram_close(fd);

    print_sizes(bsize_h, ds_h, dc_h, bmax);

    HashTableDestroy(bsize_h, 0);
    HashTableDestroy(ds_h, 0);
    HashTableDestroy(dc_h, 0);

    return 0;
}
