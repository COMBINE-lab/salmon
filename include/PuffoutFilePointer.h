//
// Created by Fatemeh Almodaresi on 5/17/18.
//

#ifndef PUFFOUTFILEPOINTER_H
#define PUFFOUTFILEPOINTER_H

#include <mutex>

typedef uint16_t rLenType;
typedef uint32_t refLenType;

class Chunk {
public:
    Chunk() {}

    // assigning 100 bytes to each read alignment group
    Chunk(uint64_t readCnt) : allocatedSize_(100 * readCnt) {}

    bool allocate(uint64_t newSize) {
        chunkSize_ = newSize;
        if (newSize > allocatedSize_) {
            allocatedSize_ = newSize;
            delete chunk_;
            chunk_ = new char(allocatedSize_);
            return true;
        }
        return false;
    }

    bool hasNext() { return currByte_ < chunkSize_; }

    template<class T>
    inline void fill(T &val) {
        memcpy(chunk_ + currByte_, &val, sizeof(val));
        currByte_ += sizeof(val);
    }

    inline void fill(std::string &val, uint64_t len) {
        char strCharBuff[len];
        memcpy(chunk_ + currByte_, &strCharBuff, len);
        val = std::string(strCharBuff, len);
        currByte_ += len;
    }

    char *chunk_ = new char(allocatedSize_);
private:
    uint64_t allocatedSize_{1000000};
    uint64_t chunkSize_{0};
    uint64_t currByte_{0};
};

class PuffoutFilePointer {
public:
    PuffoutFilePointer(std::string filename) {
        inFile.open(filename, std::ios::binary);
        if (!readHeader()) {
            logger->error("Invalid header for mapping output file.");
            std::exit(1);
        }
    }

    bool readChunk(Chunk &alnChunk) {
        if (hasNext()) {
            iomutex.lock();
            uint64_t chunksize;
            inFile.read(reinterpret_cast<char *>(&chunksize), sizeof(uint64_t));
            alnChunk.allocate(chunksize); // only allocates new space if chunksize > chunk.size()
            inFile.read(alnChunk.chunk_, chunksize);
            iomutex.unlock();
            return true;
        }
        return false;
    }

    const std::string &refName(size_t id) { return refNames[id]; }

    size_t refLength(size_t id) { return refLengths[id]; }

    size_t numRefs() const { return refNames.size(); }

    bool isMappingPaired() { return isPaired; }

private:
    bool hasNext() { return inFile.is_open() && inFile.good(); }

    bool readHeader() {
        iomutex.lock();
        size_t refCount;
        inFile.read(reinterpret_cast<char *>(&isPaired), sizeof(bool));
        inFile.read(reinterpret_cast<char *>(&refCount), sizeof(size_t));
        logger->info("Total # of References: {}", refCount);
        refLengths.reserve(refCount);
        refNames.reserve(refCount);
        uint8_t refNameSize;
        refLenType refLen;
        //std::cout << "is paired: " << isPaired << "\n";
        for (size_t i = 0; i < refCount; i++) {
            inFile.read(reinterpret_cast<char *>(&refNameSize), sizeof(refNameSize));
            char *strChar = new char[refNameSize];
            inFile.read(strChar, refNameSize);
            std::string refName(strChar, refNameSize);
            inFile.read(reinterpret_cast<char *>(&refLen), sizeof(refLenType));
            refNames.push_back(refName);
            refLengths.push_back(refLen);
            //std::cout << refName << " " << refLen << "\n";
        }
        iomutex.unlock();
        return true;
    }


    std::mutex iomutex;
    std::ifstream inFile;
    bool isPaired{true};
    std::vector <refLenType> refLengths;
    std::vector <std::string> refNames;
    std::shared_ptr <spdlog::logger> logger;

};

#endif //PUFFOUTFILEPOINTER_H
