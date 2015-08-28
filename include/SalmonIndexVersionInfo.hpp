#ifndef __SALMON_INDEX_VERSION_INFO_HPP__
#define __SALMON_INDEX_VERSION_INFO_HPP__

#include "spdlog/details/format.h"
#include "boost/filesystem.hpp"
#include "cereal/archives/json.hpp"

enum class IndexType : uint8_t { FMD, QUASI };


class SalmonIndexVersionInfo {
    public:
        /**
         * default constructor(s)
         */
       SalmonIndexVersionInfo() : indexVersion_(0), hasAuxKmerIndex_(false),
                                  auxKmerLength_(0), indexType_(IndexType::FMD) {}

       SalmonIndexVersionInfo(uint32_t indexVersionIn, bool hasAuxKmerIndexIn,
                              uint32_t auxKmerLengthIn, IndexType indexTypeIn) :
			      indexVersion_(indexVersionIn),
                              hasAuxKmerIndex_(hasAuxKmerIndexIn), auxKmerLength_(auxKmerLengthIn),
			      indexType_(indexTypeIn) {}

       /**
        * Read the index version info from file
        */
       bool load(boost::filesystem::path& versionFile) {
           namespace bfs = boost::filesystem;
            if(!bfs::exists(versionFile)) {
                fmt::MemoryWriter infostr;
                infostr << "Error: The index version file " << versionFile.string()
                        << " doesn't seem to exist.  Please try re-building the salmon "
                           "index.";
                throw std::invalid_argument(infostr.str());
            }
            std::ifstream ifs(versionFile.string());
            {
                cereal::JSONInputArchive iarchive(ifs); // Create an input archive
                iarchive(cereal::make_nvp("indexVersion", indexVersion_),
                        cereal::make_nvp("hasAuxIndex", hasAuxKmerIndex_),
                        cereal::make_nvp("auxKmerLength", auxKmerLength_),
            			cereal::make_nvp("indexType", indexType_));
            }
            ifs.close();
            return true;
       }

       bool save(boost::filesystem::path& versionFile) {
           std::ofstream ofs(versionFile.string());
           {
               cereal::JSONOutputArchive oarchive(ofs);
               oarchive(cereal::make_nvp("indexVersion", indexVersion_),
                       cereal::make_nvp("hasAuxIndex", hasAuxKmerIndex_),
                       cereal::make_nvp("auxKmerLength", auxKmerLength_),
		               cereal::make_nvp("indexType", indexType_));
           }
           ofs.close();
           return true;
       }

       bool hasAuxKmerIndex() { return hasAuxKmerIndex_; }
       void hasAuxKmerIndex(bool val) { hasAuxKmerIndex_ = val; }

       uint32_t indexVersion() { return indexVersion_; }
       void indexVersion(uint32_t version) { indexVersion_ = version; }

       uint32_t auxKmerLength() { return auxKmerLength_; }
       void auxKmerLength(uint32_t len) { auxKmerLength_ = len; };

       IndexType indexType() { return indexType_; }
       void indexType(IndexType indexTypeIn) { indexType_ = indexTypeIn; };
    private:
        uint32_t indexVersion_;
        bool hasAuxKmerIndex_;
        uint32_t auxKmerLength_;
        IndexType indexType_;
};

#endif // __SALMON_INDEX_VERSION_INFO_HPP__

