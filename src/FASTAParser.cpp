#include <unistd.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <random>

#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

#include "FASTAParser.hpp"
#include "Transcript.hpp"
#include "SalmonStringUtils.hpp"
#include "SalmonOpts.hpp"

FASTAParser::FASTAParser(const std::string& fname): fname_(fname) {}

void FASTAParser::populateTargets(std::vector<Transcript>& refs,
                                  SalmonOpts& sopt,
                                  std::unordered_map<uint32_t,uint32_t> & alleleToSuperTxpMap) {
    using stream_manager = jellyfish::stream_manager<std::vector<std::string>::const_iterator>;
    using single_parser = jellyfish::whole_sequence_parser<stream_manager>;

    using std::string;
    using std::unordered_map;

    unordered_map<string, size_t> nameToID;
    for (auto& ref : refs) {
	nameToID[ref.RefName] = ref.id;
    }

    // Separators for the header (default ' ' and '\t')
    // If we have the gencode flag, then add '|'.
    std::string sepStr = " \t";
    if (sopt.gencodeRef) {
        sepStr += '|';
    }

    std::vector<std::string> readFiles{fname_};
    size_t maxReadGroup{1000}; // Number of files to read simultaneously
    size_t concurrentFile{1}; // Number of reads in each "job"
    stream_manager streams(readFiles.cbegin(), readFiles.cend(), concurrentFile);
    single_parser parser(4, maxReadGroup, concurrentFile, streams);

    constexpr char bases[] = {'A', 'C', 'G', 'T'};
    // Create a random uniform distribution
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_int_distribution<> dis(0, 3);
    uint64_t numNucleotidesReplaced{0};

    // All header names we encounter in the fasta file
    std::unordered_set<std::string> fastaNames;

    while(true) {
        typename single_parser::job j(parser); // Get a job from the parser: a bunch of read (at most max_read_group)
        if(j.is_empty()) break;           // If got nothing, quit

        for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read we got
            std::string& header = j->data[i].header;
            std::string name = header.substr(0, header.find_first_of(sepStr));
            fastaNames.insert(name);

            auto it = nameToID.find(name);
            if (it == nameToID.end()) {
              sopt.jointLog->warn("Transcript {} appears in the reference but did not appear in the BAM", name);
            } else {

	      std::string& seq = j->data[i].seq;
              size_t readLen = seq.length();

	      refs[it->second].setSAMSequenceOwned(salmon::stringtools::encodeSequenceInSAM(seq.c_str(), readLen));

	      // Replace non-ACGT bases
	      for (size_t b = 0; b < readLen; ++b) {
		seq[b] = ::toupper(seq[b]);
		int c = jellyfish::mer_dna::code(seq[b]);
		// Replace non-ACGT bases with pseudo-random bases
		if (jellyfish::mer_dna::not_dna(c)) {
		  char rbase = bases[dis(eng)];
		  c = jellyfish::mer_dna::code(rbase);
		  seq[b] = rbase;
		  ++numNucleotidesReplaced;
		}
	      }

	      // allocate space for the new copy
	      char* seqCopy = new char[seq.length()+1];
	      std::strcpy(seqCopy, seq.c_str());
	      refs[it->second].setSequenceOwned(seqCopy, sopt.gcBiasCorrect, sopt.gcSampFactor);
	      // seqCopy will only be freed when the transcript is destructed!
            }
        }
    }

    // Check that every sequence present in the BAM header was also present in the
    // transcriptome fasta.
    bool missingTxpError{false};
    for (auto& kv : nameToID) {
      auto& name = kv.first;
      if (fastaNames.find(name) == fastaNames.end()) {
        sopt.jointLog->critical("Transcript {} appeared in the BAM header, but was not in the provided FASTA file", name);
        missingTxpError = true;
      }
    }

    //////////////////////////////////////////////////////////
    //
    //Avi's Edits
    std::ifstream alleleFile(sopt.alleleFilePath.c_str());
    uint32_t avgLength, origNoOfTxps, matIndex, patIndex, catIndex, i=0;
    std::set<uint32_t> indexList;
    origNoOfTxps = refs.size();
    for (std::string token; alleleFile >> token; ++i) {
        matIndex = nameToID[token + "_mat"];
        patIndex = nameToID[token + "_pat"];
        indexList.insert(matIndex);
        indexList.insert(patIndex);
        catIndex = origNoOfTxps+i;
        avgLength = (refs[matIndex].RefLength + refs[patIndex].RefLength) / 2;
        refs.emplace_back(catIndex, (token+"_cat").c_str(), avgLength, 0.005);
        alleleToSuperTxpMap.insert(std::pair<uint32_t, uint32_t>(matIndex, catIndex));
        alleleToSuperTxpMap.insert(std::pair<uint32_t, uint32_t>(patIndex, catIndex));
    }

    std::set<uint32_t>::iterator it;
    for (uint32_t i=0; i<nameToID.size(); ++i){
        it=indexList.find(i);
        if (it == indexList.end()) {
            alleleToSuperTxpMap.insert(std::pair<uint32_t, uint32_t>(i,i));
        }
        else{
            indexList.erase(it);
        }
    }
    //////////////////////////////////////////////////////////

    if (missingTxpError) {
      sopt.jointLog->critical("Please provide a reference FASTA file that includes all targets present in the BAM header\n"
      "If you have access to the genome FASTA and GTF used for alignment \n"
      "consider generating a transcriptome fasta using a command like: \n"
      "gffread -w output.fa -g genome.fa genome.gtf\n"
      "you can find the gffread utility at (http://ccb.jhu.edu/software/stringtie/gff.shtml)");
      sopt.jointLog->flush();
      std::exit(1);
    }

    sopt.jointLog->info("replaced {} non-ACGT nucleotides with random nucleotides",  numNucleotidesReplaced);

}

