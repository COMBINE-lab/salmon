#include "gdna.h"
#include <string.h>

const char* IUPAC_2BIT  ="AACCTTGGTTAAAAAACCCCGGAAAAAACCAAAAAA";
const char* IUPAC_2BITN ="001133223300000011112200000011000000";
const char* IUPAC_DEFS ="AaCcTtGgUuMmRrWwSsYyKkVvHhDdBbNnXx-*";
const char* IUPAC_COMP ="TtGgAaCcAaKkYyWwSsRrMmBbDdHhVvNnXx-*";

#define A_2BIT 0 // 00
#define C_2BIT 1 // 01
#define G_2BIT 2 // 10
#define T_2BIT 3 // 11

static byte ntCompTable[256];
static byte nt2bit[256]; //maps any character to a 2bit base value (with N = A)
static char v_2bit2nt[4] = {'A','C','G','T'};

//----------------------

static bool gdna_Ready=gDnaInit();

//----------------------

byte gdna2bit(char* &nt, int n) {
// Pack n bases into a byte (n can be 1..4)
byte out = 0;
while (n && *nt) {
    n--;
    out <<= 2;
    out += nt2bit[(int)*nt];
    nt++;
    }
#ifdef GDEBUG
if (n) {
     GError("Error: attempt to read 6-mer beyond the end of the string!\n");
     }
#endif
return out;
}


char ntComplement(char c) {
 return ntCompTable[(int)c];
 }

char g2bit2base(byte v2bit) {
 return v_2bit2nt[v2bit & 0x03 ];
}

//in place reverse complement of nucleotide (sub)sequence
char* reverseComplement(char* seq, int slen) {
   if (slen==0) slen=strlen(seq);
   //reverseChars(seq,len);
   int l=0;
   int r=slen-1;
   char c;
   while (l<r) {
      c=seq[l];seq[l]=seq[r];
      seq[r]=c;   //this was: Gswap(str[l],str[r]);
      l++;r--;
      }
   for (int i=0;i<slen;i++) seq[i]=ntComplement(seq[i]);
   return seq;
 }

bool gDnaInit() {
       if (gdna_Ready) return true;
       int l=strlen(IUPAC_DEFS);
       ntCompTable[0]=0;
       nt2bit[0]=0;
       for (int ch=1;ch<256;ch++) {
          ntCompTable[ch]=0;
          nt2bit[ch]=0;
          for (int i=0;i<l;i++)
                if (ch==IUPAC_DEFS[i]) {
                  ntCompTable[ch]=IUPAC_COMP[i];
                  nt2bit[ch] = IUPAC_2BITN[i]-'0';
                  break;
                  }
          if (ntCompTable[ch]==0) {
              ntCompTable[ch]='N';
              }
          }
      gdna_Ready=true;
      return true;
     }




