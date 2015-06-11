#include <unistd.h>
#include <cstdio>
#include <iostream>
#include <unordered_map>


#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

#include "FASTAParser.hpp"
#include "Transcript.hpp"
#include "SalmonStringUtils.hpp"

FASTAParser::FASTAParser(const std::string& fname): fname_(fname) {}

void FASTAParser::populateTargets(std::vector<Transcript>& refs) {
    using stream_manager = jellyfish::stream_manager<std::vector<std::string>::const_iterator>;
    using single_parser = jellyfish::whole_sequence_parser<stream_manager>;

    using std::string;
    using std::unordered_map;

    unordered_map<string, size_t> nameToID;
    for (auto& ref : refs) {
	nameToID[ref.RefName] = ref.id;
    }

    std::vector<std::string> readFiles{fname_};
    size_t maxReadGroup{1000}; // Number of files to read simultaneously
    size_t concurrentFile{1}; // Number of reads in each "job"
    stream_manager streams(readFiles.cbegin(), readFiles.cend(), concurrentFile);
    single_parser parser(4, maxReadGroup, concurrentFile, streams);

    while(true) {
        typename single_parser::job j(parser); // Get a job from the parser: a bunch of read (at most max_read_group)
        if(j.is_empty()) break;           // If got nothing, quit

        for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read we got
            std::string& header = j->data[i].header;
            std::string name = header.substr(0, header.find(' '));

            auto it = nameToID.find(name);
            if (it == nameToID.end()) {
                std::cerr << "WARNING: Transcript " << name << " appears in the reference but did not appear in the BAM\n";
            } else {
                std::string& seq = j->data[i].seq;
                refs[it->second].Sequence = salmon::stringtools::encodeSequenceInSAM(seq.c_str(), seq.size());
            }
        }
    }

}

