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
    Chunk() {
        chunk_.resize(allocatedSize_);
    }

    // assigning 100 bytes to each read alignment group
    Chunk(uint64_t readCnt) : allocatedSize_(100 * readCnt) {
        chunk_.resize(allocatedSize_);
    }

    bool allocate(uint64_t newSize) {
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

class PuffoutFilePointer {
public:
    PuffoutFilePointer(std::string filename) {
        logger = spdlog::get("jointLog");
        logger->info("reading from pufferfish output: {}", filename);
        inFile.open(filename, std::ios::binary);
        if (!readHeader()) {
            logger->error("Invalid header for mapping output file.");
            std::exit(1);
        }
    }

    bool readChunk(Chunk &alnChunk, std::mutex &iomutex) {
        std::lock_guard<std::mutex> l(iomutex);
        if (hasNext()) {
            //iomutex.lock();
            uint64_t chunksize;
            inFile.read(reinterpret_cast<char *>(&chunksize), sizeof(chunksize));
            if (!hasNext()) return false; // because you need to read last chunk from file first before the flag is set
            //std::cerr << "chunk size, new: " << chunksize << " cur: " << alnChunk.chunk_.size() << "\n";
            //if (chunksize == 0) return false; // Fixme there is no guarantee that you'll read 0 from file
            alnChunk.allocate(chunksize); // only allocates new space if chunksize > chunk.size()
            //std::cerr << "alnChunk.chunk_.size(): " << alnChunk.chunk_.size() << "\n";
            inFile.read(alnChunk.chunk_.data(), chunksize);
            //iomutex.unlock();
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
        //iomutex.lock();
        size_t refCount;
        inFile.read(reinterpret_cast<char *>(&isPaired), sizeof(bool));
        inFile.read(reinterpret_cast<char *>(&refCount), sizeof(size_t));
        logger->info("Total # of References: {}", refCount);
        refLengths.reserve(refCount);
        refNames.reserve(refCount);
        uint8_t refNameSize;
        refLenType refLen;
        //std::cerr << "is paired: " << isPaired << "\n";
        for (size_t i = 0; i < refCount; i++) {
            inFile.read(reinterpret_cast<char *>(&refNameSize), sizeof(refNameSize));
            char *strChar = new char[refNameSize];
            inFile.read(strChar, refNameSize);
            std::string refName(strChar, refNameSize);
            inFile.read(reinterpret_cast<char *>(&refLen), sizeof(refLenType));
            refNames.push_back(refName);
            refLengths.push_back(refLen);
            //std::cerr << refName << " " << refLen << "\n";
        }
        //iomutex.unlock();
        return true;
    }


    //std::mutex iomutex;
    std::ifstream inFile;
    bool isPaired{false};
    std::vector <refLenType> refLengths;
    std::vector <std::string> refNames;
    std::shared_ptr <spdlog::logger> logger;

};

#endif //PUFFOUTFILEPOINTER_H
