#ifndef FASTA_PARSER
#define FASTA_PARSER

#include <vector>

class Transcript;
struct SalmonOpts;

class FASTAParser {
public:
  FASTAParser(const std::string& fname);
  void populateTargets(std::vector<Transcript>& transcripts, SalmonOpts& sopt);

private:
  std::string fname_;
};

#endif // FASTA_PARSER
