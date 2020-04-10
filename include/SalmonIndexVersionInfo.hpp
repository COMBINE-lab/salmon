#ifndef __SALMON_INDEX_VERSION_INFO_HPP__
#define __SALMON_INDEX_VERSION_INFO_HPP__

#include "SalmonConfig.hpp"
#include "boost/filesystem.hpp"
#include "cereal/archives/json.hpp"
#include "spdlog/fmt/fmt.h"
#include "json.hpp"

enum class SalmonIndexType : uint8_t { FMD=0, QUASI=1, PUFF=2 };

class SalmonIndexVersionInfo {
public:
  /**
   * default constructor(s)
   */
  SalmonIndexVersionInfo()
      : indexVersion_(0), hasAuxKmerIndex_(false), auxKmerLength_(0),
        indexType_(SalmonIndexType::PUFF) {}

  SalmonIndexVersionInfo(uint32_t indexVersionIn, bool hasAuxKmerIndexIn,
                         uint32_t auxKmerLengthIn, SalmonIndexType indexTypeIn, std::string sverIn)
      : indexVersion_(indexVersionIn), hasAuxKmerIndex_(hasAuxKmerIndexIn),
        auxKmerLength_(auxKmerLengthIn), indexType_(indexTypeIn), salmonVersion_(sverIn) {}

  /**
   * Read the index version info from file
   */
  bool load(boost::filesystem::path& versionFile) {
    namespace bfs = boost::filesystem;
    if (!bfs::exists(versionFile)) {
      fmt::MemoryWriter infostr;
      infostr << "Error: The index version file " << versionFile.string()
              << " doesn't seem to exist.  Please try re-building the salmon "
                 "index.";
      throw std::invalid_argument(infostr.str());
    }
    std::ifstream ifs(versionFile.string());
    {
      // NOTE: Let's keep this around for one version (until 1.3.0)
      // just for reference.
      // previous cereal-based implementation to load version information metadata
      /*
      cereal::JSONInputArchive iarchive(ifs); // Create an input archive
      iarchive(cereal::make_nvp("indexVersion", indexVersion_),
               cereal::make_nvp("hasAuxIndex", hasAuxKmerIndex_),
               cereal::make_nvp("auxKmerLength", auxKmerLength_),
               cereal::make_nvp("indexType", indexType_),
               cereal::make_nvp("salmonVersion", salmonVersion_));
      */

      nlohmann::json j;
      ifs >> j;

      // Right now, we will make the reading of this optional 
      // as we don't want to break the backward compatibility 
      // of salmon 1.2.0 to _read_ indices made with earlier 
      // versions.
      if (j.find("salmonVersion") != j.end()) {
        salmonVersion_ = j["salmonVersion"];
      } else {
        salmonVersion_ = "0.0.0";        
      }

      indexVersion_ = j["indexVersion"].get<uint32_t>();
      hasAuxKmerIndex_ = j["hasAuxIndex"].get<bool>();
      auxKmerLength_ = j["auxKmerLength"].get<uint32_t>();
      // This is one (small) place where cereal handles things
      // more nicely than nlohmann::json.
      switch(j["indexType"].get<uint8_t>()) {
        case 0:
          indexType_ = SalmonIndexType::FMD;
          break;
        case 1:
          indexType_ = SalmonIndexType::QUASI;
          break;
        case 2:
          indexType_ = SalmonIndexType::PUFF;
          break;
        default: {
          fmt::MemoryWriter infostr;
          infostr << "Unknown index type tag : " << j["indexType"].get<uint32_t>() << ".";
          throw std::invalid_argument(infostr.str());
        }
      }

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
               cereal::make_nvp("indexType", indexType_),
               cereal::make_nvp("salmonVersion", salmonVersion_));
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

  SalmonIndexType indexType() { return indexType_; }
  void indexType(SalmonIndexType indexTypeIn) { indexType_ = indexTypeIn; };

  std::string salmonVersion() const { return salmonVersion_; }
  void salmonVersion(const std::string& sv) { salmonVersion_ = sv; }

private:
  uint32_t indexVersion_;
  bool hasAuxKmerIndex_;
  uint32_t auxKmerLength_;
  SalmonIndexType indexType_;
  std::string salmonVersion_;
};

#endif // __SALMON_INDEX_VERSION_INFO_HPP__
