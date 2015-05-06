#ifndef FASTA_PARSER
#define FASTA_PARSER

#include <vector>

class Transcript;

class FASTAParser {
public:
    FASTAParser(const std::string& fname);
    void populateTargets(std::vector<Transcript>& transcripts);
private:
    std::string fname_;
};

#endif // FASTA_PARSER 
