#ifndef ARITH_STATIC_H
#define ARITH_STATIC_H

unsigned char *arith_compress(unsigned char *in, unsigned int in_size,
			      unsigned int *out_size, int order);
unsigned char *arith_uncompress(unsigned char *in, unsigned int in_size,
				unsigned int *out_size, int order);


#endif /* ARITH_STATIC_H */
