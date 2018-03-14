#include <iostream>

extern "C" {
#include "kseq.h"
}

#include <zlib.h>

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

namespace v1Converter {

struct ReadSeq {
  std::string seq;
  std::string name;
  ~ReadSeq() {}
};

inline void copyRecord(kseq_t* seq, ReadSeq* s) {
  // Copy over the sequence and read name
  s->seq.assign(seq->seq.s, seq->seq.l);
  s->name.assign(seq->name.s, seq->name.l);
}
}

int main(int argc, char* argv[]) {
  using namespace v1Converter;
  // open the file and init the parser
  auto fp = gzopen(argv[1], "r");
  auto fp2 = gzopen(argv[2], "r");

  ReadSeq bc;
  ReadSeq seqFrag;
  ReadSeq umi;

  auto* seq = kseq_init(fp);
  auto* seq2 = kseq_init(fp2);

  int ksv = kseq_read(seq);
  int ksv2 = kseq_read(seq2);
  while (ksv >= 0 and ksv2 >= 0) {
    // The right read becomes the read
    copyRecord(seq2, &seqFrag);
    ksv2 = kseq_read(seq2);
    if (ksv2 >= 0) {
      // The left read becomes the barcode + umi
      bc.name.assign(seq->name.s, seq->name.l);
      bc.seq.assign(seq->seq.s, seq->seq.l);
      bc.seq.insert(seq->seq.l, seq2->seq.s, seq2->seq.l);
    } else {
      break;
    }


    std::cout << '>' << bc.name << "\n";
    std::cout << bc.seq<< "\n";

    std::cerr << '>' << seqFrag.name << "\n";
    std::cerr << seqFrag.seq << "\n";

    ksv = kseq_read(seq);
    ksv2 = kseq_read(seq2);
  }
  // destroy the parser and close the file
  kseq_destroy(seq);
  gzclose(fp);
  kseq_destroy(seq2);
  gzclose(fp2);

  return 0;
}
