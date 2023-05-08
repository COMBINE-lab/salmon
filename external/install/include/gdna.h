#ifndef GDNA_H
#define GDNA_H
#include "GBase.h"

extern const char* IUPAC_DEFS;
extern const char* IUPAC_COMP;

char ntComplement(char c);

//in-place reverse complement of a nucleotide (sub)sequence
char* reverseComplement(char* seq, int slen=0);

bool gDnaInit();

byte gdna2bit(char* &nt, int n=4); //pack n bases into a byte (n can be 1..4)
char g2bit2base(byte v2bit); //convert the 2-bit value into 'A', 'C', 'G' or 'T'

#endif
