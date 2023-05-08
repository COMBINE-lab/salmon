#include <utility>

#ifndef __MAPPINGS_H__
#define __MAPPINGS_H__

#include "spdlog/spdlog.h"
#include <fstream>
#include <iostream>
#include <mutex>
#include "Taxa.h"

class Chunk {
public:
    Chunk() {
        chunk_.resize(allocatedSize_);
    }

    // assigning 100 bytes to each read alignment group
    Chunk(uint64_t readCnt) : allocatedSize_(100 * readCnt) {
        chunk_.resize(allocatedSize_);
    }

    bool allocate(uint64_t newSize) {
//        std::cerr << "new chunk size: " << newSize << "\n";
        chunkSize_ = newSize;
        currByte_ = 0;
        if (newSize > allocatedSize_) {
            allocatedSize_ = newSize;
            chunk_.clear();
            chunk_.resize(allocatedSize_);
            //delete[] chunk_;
            //chunk_ = new char[allocatedSize_];
            return true;
        }
        return false;
    }

    bool hasNext() { return currByte_ < chunkSize_; }

    template<class T>
    inline void fill(T &val) {
        memcpy(&val, chunk_.data() + currByte_, sizeof(val));
        currByte_ += sizeof(val);

        //std::cerr << "fill: "<< (uint64_t) val << "," << currByte_ << "\n";
    }

    inline void fill(std::string &val, uint64_t len) {
        char strCharBuff[len];
        memcpy(&strCharBuff, chunk_.data() + currByte_, len);
        val = std::string(strCharBuff, len);
        currByte_ += len;
        //std::cerr << "fill: " << len << "," << val << "," << currByte_ << "," << chunkSize_ <<"\n";
    }

    std::vector<char> chunk_;
    //char *chunk_ = new char[allocatedSize_];
private:
    uint64_t allocatedSize_{1000000};
    uint64_t chunkSize_{0};
    uint64_t currByte_{0};

};

class PAMReader {
public:
    PAMReader(std::string mfileName,
              std::shared_ptr<spdlog::logger> loggerIn) {
        logger = loggerIn;
        if (!inFile.is_open()) {
            inFile.open(mfileName, std::ios::binary);
        }
        if (!readHeader()) {
            logger->error("Invalid header for mapping output file.");
            std::exit(1);
        }
    }
    bool readHeader() {
        size_t refCount;
        inFile.read(reinterpret_cast<char *>(&isPaired), sizeof(bool));
        inFile.read(reinterpret_cast<char *>(&refCount), sizeof(size_t));
        logger->info("Total # of References: {}", refCount);
        refLengths.resize(refCount);
        refNames.resize(refCount);
        uint8_t refNameSize;
        refLenType refLen;
        //std::cout << "is paired: " << isPaired << "\n";
        for (size_t i = 0; i < refCount; i++) {
            inFile.read(reinterpret_cast<char *>(&refNameSize), sizeof(refNameSize));
            char *strChar = new char[refNameSize];
            inFile.read(strChar, refNameSize);
            std::string refName(strChar, refNameSize);
            inFile.read(reinterpret_cast<char *>(&refLen), sizeof(refLenType));
            refNames[i] = refName;
            refLengths[i] = refLen;
//            std::cerr << refName << " " << refLen << "\n";
        }
        return true;
    }

    bool readChunk(Chunk &chunk, std::mutex& iomutex) {
        std::lock_guard<std::mutex> l(iomutex);
        if (hasNext()) {
            uint64_t chunksize;
            inFile.read(reinterpret_cast<char *>(&chunksize), sizeof(chunksize));
//            std::cerr << "chunkSize: " << chunksize << "\n";
            if (!hasNext()) return false; // because you need to read last chunk from file first before the flag is set
            chunk.allocate(chunksize); // only allocates new space if chunksize > chunk.size()
            inFile.read(chunk.chunk_.data(), chunksize);
            return true;
        }
        return false;
    }

    bool hasNext() { return inFile.is_open() && inFile.good(); }


    std::ifstream inFile;
    bool isPaired = true;
    std::vector<refLenType> refLengths;
    std::vector<std::string> refNames;
    std::shared_ptr<spdlog::logger> logger;

};

class PuffMappingReader {
public:
    PuffMappingReader() {}

    PuffMappingReader(PAMReader* pr,
                      std::shared_ptr<spdlog::logger> loggerIn) {
        pamReader = pr;
        logger = loggerIn;
    }

    void load(std::string mfileName,
              std::shared_ptr<spdlog::logger> loggerIn) {
        logger = loggerIn;
        pamReader = new PAMReader(std::move(mfileName), loggerIn);
    }

    PAMReader* getReader() {return pamReader;}

    bool nextAlignmentGroup(std::vector<ReadInfo> &alignmentGrp,
                            std::mutex &iomutex,
                            uint32_t /*threadID*/,
                            bool needReadName = false) {
        uint32_t readsLeft{0};
        for (ReadInfo& rinfo: alignmentGrp) {
            if (nextRead(rinfo, iomutex, needReadName))
                readsLeft++;
            else
                break;
        }
        /*{
            std::lock_guard<std::mutex> l(iomutex);
            std::cerr << "Thread " << threadID << " got " << readsLeft << " in this cycle\n";
        }*/
        if (readsLeft < alignmentGrp.size())
            alignmentGrp.resize(readsLeft);
        return readsLeft > 0;
    }

    bool nextRead(ReadInfo &rinf, std::mutex& iomutex, bool needReadName = false) {
        if (!chunk.hasNext()) {// try to read a new chunk from file
            if (!pamReader->readChunk(chunk, iomutex)) // if nothing left to process return false
                return false;
        }
        rinf.mappings.clear();

        uint32_t puff_id, mcnt{0};
        rLenType rlen{0}, lcnt{0}, rcnt{0};//, ibeg{0}, ilen{0};
        refLenType refPos;
        std::string readName;
        uint8_t readNameLength;
        // read name
        chunk.fill(readNameLength);
        chunk.fill(readName, readNameLength);

        // NOTE: do we actually care about the read name?
        if (needReadName) {
            rinf.rid = readName;
        }

        // mapping count
        chunk.fill(mcnt);
        chunk.fill(rlen); // left read len
/*
        {
            std::lock_guard<std::mutex> l(iomutex);
            if (mcnt == 0)
                std::cerr << "\n\ncnt is zerooooo\n\n";
        }
*/

        rinf.cnt = mcnt;
        rinf.len = rlen;
        if (pamReader->isPaired) {
            chunk.fill(rlen); // right read len
            rinf.len += rlen;
        }
//            std::cerr << rinf.rid << " " << rinf.cnt << " " << rinf.len << "\n";
        rinf.mappings.reserve(rinf.cnt);
        TaxaNode dummy;
        TaxaNode *taxaPtr = &dummy;

        for (size_t mappingCntr = 0; mappingCntr < rinf.cnt; mappingCntr++) {
            double rscore{0.0}, lscore{0.0};
            chunk.fill(puff_id);
            // fetch the taxon from the map
//                std::cerr << "puff_id: " << puff_id << " ";
            rinf.mappings.emplace_back(puff_id);
            taxaPtr = &rinf.mappings.back();
            chunk.fill(lcnt);
//                std::cerr << "left intrvls:" << lcnt << " ";
            if (pamReader->isPaired) {
                chunk.fill(rcnt);
//                    std::cerr << "right intrvls:" << rcnt << " ";
            }
            if (lcnt) {
                chunk.fill(lscore);
                chunk.fill(refPos);
                taxaPtr->setFw(refPos & PuffMappingReader::HighBitMask, ReadEnd::LEFT);
                taxaPtr->setPos(refPos & PuffMappingReader::LowBitsMask, ReadEnd::LEFT);
//                std::cerr << "left score:" << lscore << " ";
//                std::cout << "left pos:" << refPos << " " << taxaPtr->getPos(ReadEnd::LEFT) << " " << taxaPtr->isFw(ReadEnd::LEFT) << "\n";
            }

            if (pamReader->isPaired) {
                if (rcnt) {
                    ReadEnd currRe = ReadEnd::RIGHT;
                    // NOTE, IMPORTANT! There is this implicit untold!!! assumption that
                    // if a mapping is orphan (one end is only mapped), we're gonna see that as the LEFT end!!
                    if (!lcnt) {
                        currRe = ReadEnd::LEFT;
                    }
                    chunk.fill(rscore);
                    chunk.fill(refPos);
                    taxaPtr->setFw(refPos & PuffMappingReader::HighBitMask, currRe);
                    taxaPtr->setPos(refPos & PuffMappingReader::LowBitsMask, currRe);
//                    std::cerr << "right score:" << rscore << " ";
//                    std::cout << "right pos: " << taxaPtr->getPos(currRe) << " " << taxaPtr->isFw(currRe) << "\n";
                }
            }
            //std::cout << "here\n";
            taxaPtr->cleanIntervals(ReadEnd::LEFT);
            taxaPtr->cleanIntervals(ReadEnd::RIGHT);
            taxaPtr->setScore(lscore + rscore);
            //taxaPtr->updateScore();
            //std::cout << "here\n";
        }
        return true;
    }

    const std::string &refName(size_t id) { return pamReader->refNames[id]; }

    size_t refLength(size_t id) { return pamReader->refLengths[id]; }

    size_t numRefs() const { return pamReader->refNames.size(); }

    bool isMappingPaired() { return pamReader->isPaired; }

private:

    Chunk chunk;
    PAMReader* pamReader;
    std::shared_ptr<spdlog::logger> logger;
    static constexpr const refLenType HighBitMask = 1u << (sizeof(refLenType) * 8 - 1);
    static constexpr const refLenType LowBitsMask = HighBitMask - 1;

    //ReadInfo readInfo;
};

#endif
