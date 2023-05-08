#ifndef __bulkChunk__
#define __bulkChunk__

#include "BinWriter.hpp"
enum class TypeId : uint8_t {
    Bool,
    U8,
    U16,
    U32,
    U64,
    F32,
    F64,
    ARR,
    STR
};

struct TagDesc {
    std::string name;
    TypeId typeId;
    TagDesc(std::string nameIn, TypeId typeIdIn) : name(nameIn), typeId(typeIdIn) {}
};

struct TagSection {
    std::vector<TagDesc> tags;
    uint16_t count() {return static_cast<uint16_t>(tags.size());}
};

class FileTags {
    protected:
    BinWriter bw;
    TagSection fileLevel;
    TagSection readLevel;
    TagSection alignmentLevel;
    
    public:
    FileTags(): bw(10000) {}

    BinWriter &export2Buffer() {
        // file-level
        bw << fileLevel.count();
        for (auto &td : fileLevel.tags) {
            bw << td.name;
            bw << static_cast<uint8_t>(td.typeId);
        }
        // read-level
        bw << readLevel.count();
        for (auto &td : readLevel.tags) {
            bw << td.name;
            bw << static_cast<uint8_t>(td.typeId);
        }
        // alignment-level
        bw << alignmentLevel.count();
        for (auto &td : alignmentLevel.tags) {
            bw << td.name;
            bw << static_cast<uint8_t>(td.typeId);
        }
        return bw;
    }
};

class BulkTags : public FileTags {
    public:
    BulkTags(bool isPE = true) : FileTags() {
        

        // standard tags for bulk
        // fileLevel.tags.emplace_back("allowOrphs", TypeId::Bool);
        // fileLevel.tags.emplace_back("allowDiscords", TypeId::Bool);
        // fileLevel.tags.emplace_back("allowDoves", TypeId::Bool);
        // no file-level tags
        // read-level tags
        readLevel.tags.emplace_back("readName", TypeId::STR);
        readLevel.tags.emplace_back("alignmentCnt", TypeId::U32);
        readLevel.tags.emplace_back("readLen", TypeId::U16);
        if (isPE) {
            readLevel.tags.emplace_back("pairedReadLen", TypeId::U16);
        }

        // alignment-level
        alignmentLevel.tags.emplace_back("refId", TypeId::STR);
        // single-end or left-end
        alignmentLevel.tags.emplace_back("leftPos", TypeId::U64);
        alignmentLevel.tags.emplace_back("leftOri", TypeId::Bool);
        alignmentLevel.tags.emplace_back("leftCIGAR", TypeId::STR);
        alignmentLevel.tags.emplace_back("leftScore", TypeId::F32);

        if (isPE) {
            // right-end
            alignmentLevel.tags.emplace_back("righPos", TypeId::U64);
            alignmentLevel.tags.emplace_back("rightOri", TypeId::Bool);
            alignmentLevel.tags.emplace_back("rightCIGAR", TypeId::STR);
            alignmentLevel.tags.emplace_back("rightScore", TypeId::F32);
        }
    }
};

struct ChunkConfig {
    uint64_t num_chunks;
};

#endif