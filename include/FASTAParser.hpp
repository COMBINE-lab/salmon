#ifndef FASTA_PARSER
#define FASTA_PARSER

#include <vector>

class Transcript;
class SalmonOpts;

class FASTAParser {
public:
    FASTAParser(const std::string& fname);
    void populateTargets(std::vector<Transcript>& transcripts, SalmonOpts& sopt,
                        std::unordered_map<uint32_t,uint32_t> & alleleToSuperTxpMap);

private:
    std::string fname_;
};

#endif // FASTA_PARSER
