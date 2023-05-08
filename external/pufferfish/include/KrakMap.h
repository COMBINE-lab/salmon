#ifndef __KRAKMAP_H__
#define __KRAKMAP_H__

#include <deque> // std::deque
//#include <queue> // std::priority_queue
#include "sparsepp/spp.h"
#include "Taxa.h"


class KrakMap {
    public:
        KrakMap(std::string& taxonomyTree_filename, 
                std::string& refId2TaxId_filename, 
                std::string pruneLevelIn, 
                double filteringThresholdIn);
        bool classify(std::string& mapperOutput_filename);
        void serialize(std::string& output_filename);
    private:
        bool readHeader(std::ifstream& mfile);
        void loadMappingInfo(std::ifstream& mfile);
        void walk2theRoot(TaxaNode* child);
        void propagateInfo();
        void assignRead(uint64_t readLen);
        void clearReadSubTree();
        
        spp::sparse_hash_map<uint32_t, TaxaNode> taxaNodeMap;
        spp::sparse_hash_map<std::string, uint32_t> refId2taxId;
        std::deque<TaxaNode*> hits;
        std::set<uint64_t> activeTaxa;
        Rank pruningLevel = Rank::SPECIES;
        uint64_t rootId = 1;
        double filteringThreshold = 0;
        spp::sparse_hash_map<uint64_t, TaxaInfo> mappedReadCntr;
        uint64_t readCntr = 0;
        bool isPaired = true;
};

#endif