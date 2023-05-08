
#include <cstdlib>
#include <string>

#include "spdlog/spdlog.h"

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/BamHeader.h"
#include "SeqLib/BamRecord.h"

bool isTruth(SeqLib::BamRecord& rec, SeqLib::BamHeader& bh) {
    if (rec.MappedFlag() or rec.PairMappedFlag()) { 
        std::string gis = bh.IDtoName(rec.ChrID());
        std::string readgis = rec.Qname();
        size_t s = gis.find_first_of('|')+1;
        size_t e = gis.find_first_of('|', s);
        gis = gis.substr(s, e-s-1);
        s = readgis.find_first_of('|')+1;
        e = readgis.find_first_of('|', s);
        readgis = readgis.substr(s, e-s-1);
        return readgis == gis;
    }
    return false;
}
void filterFractionOffp(SeqLib::BamReader& br,
                 SeqLib::BamWriter& bw,
                 SeqLib::BamHeader& bh,
                 float prob,
                 std::shared_ptr<spdlog::logger>& logger) {
    uint64_t cntr = 0;
    SeqLib::BamRecord rec;
    SeqLib::BamRecord rec2;
    
    while (br.GetNextRecord(rec)) {        
        if ( (rec.MappedFlag() or rec.PairMappedFlag()) ) {
            if ( !isTruth(rec, bh) and ((double) rand() / (RAND_MAX)) <= prob ) {
                if (rec.PairedFlag())
                    br.GetNextRecord(rec2);
                continue;
            }
        }
        bw.WriteRecord(rec);
        cntr++;
        if (rec.PairedFlag()) {
            br.GetNextRecord(rec2);
            bw.WriteRecord(rec2);
            cntr++;
        }
        if (cntr % 10000000 == 0) {
            logger->info("{} records processed.", cntr);
        }
    }
}

void filterAllFp(SeqLib::BamReader& br,
                 SeqLib::BamWriter& bw,
                 SeqLib::BamHeader& bh,
                 std::shared_ptr<spdlog::logger>& logger) {
    uint64_t cntr = 0;
    SeqLib::BamRecord rec;
    SeqLib::BamRecord rec2;
    bool hasTruth = false;
    std::string curReadName = "dummy";
    SeqLib::BamRecordVector recvec;
    recvec.reserve(3000);
    while (br.GetNextRecord(rec)) {
        cntr++;
        if (rec.Qname() != curReadName) {
            curReadName = rec.Qname();
            if (hasTruth) {
                for (auto& rec_it : recvec) {
                    bw.WriteRecord(rec_it);
                }
            }
            recvec.clear();
            hasTruth = false;
        }
        if (isTruth(rec, bh)) {
            hasTruth = true;
        }        
        recvec.push_back(rec);  
        if (cntr % 10000000 == 0) {
            logger->info("{} records processed.", cntr);
        }           
    }
    // for last set of records
    if (hasTruth) {
        for (auto& rec_it : recvec) {
            bw.WriteRecord(rec_it);
        }
    }
    
}

int main(int argc, char* argv[]) {
    std::shared_ptr<spdlog::logger> logger = spdlog::stderr_color_mt("console");
    if (argc < 4) {
        logger->error("choose one of the filter types 1. all_fps 2. part_fps\n and at least two args, input and output BAM file names");
        exit(1);
    }

    std::string command = argv[1];
    std::string mfileName = argv[2];
    std::string moutfileName = argv[3];
    
    logger->info("command: {}", command);
    logger->info("input: {}", mfileName);
    logger->info("output: {}", moutfileName);

    SeqLib::BamReader br;
    SeqLib::BamWriter bw;
    br.Open(mfileName);   
    bw.Open(moutfileName); 
    SeqLib::BamHeader bh = br.Header();   
    bw.SetHeader(br.Header());  
    bw.WriteHeader();       
    logger->info("# of targets: {}", bh.NumSequences());            
    
    if (command == "all_fps") {
        filterAllFp(br, bw, bh, logger);    
    }
    else if (command == "part_fps") {
        if (argc < 5) {
            logger->error("require at least 3 input arguments: input bam (str), output bam (str), filtering threshold (number between 0 and 1)");
            exit(1);
        }
        float prob = std::stof(argv[4]);
        logger->info("filter thresh: {}", prob);    
        filterFractionOffp(br, bw, bh, prob, logger); 
    }
    
    br.Close();
    bw.Close();
}