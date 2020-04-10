#ifndef __SALMON_INDEX_HPP__
#define __SALMON_INDEX_HPP__

#include <memory>

#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

#include "cereal/archives/json.hpp"
#include "cereal/types/vector.hpp"
#include "spdlog/spdlog.h"
#include "json.hpp"

#include "pufferfish/ProgOpts.hpp"
#include "pufferfish/PufferfishIndex.hpp"
#include "pufferfish/PufferfishSparseIndex.hpp"

#include "SalmonUtils.hpp"
#include "SalmonConfig.hpp"
#include "SalmonIndexVersionInfo.hpp"

// declaration of quasi index function
int pufferfishIndex(pufferfish::IndexOptions& indexOpts);

class SalmonIndex {
public:
  SalmonIndex(std::shared_ptr<spdlog::logger>& logger,
              SalmonIndexType indexType)
      : loaded_(false), versionInfo_(0, false, 0, indexType, ""), logger_(logger),
        seqHash256_(""), nameHash256_(""), seqHash512_(""), nameHash512_(""),
        decoySeqHash256_(""), decoyNameHash256_("") {}

  void load(const boost::filesystem::path& indexDir) {
    namespace bfs = boost::filesystem;

    // Check if version file exists and, if so, read it.
    boost::filesystem::path versionPath = indexDir / "versionInfo.json";
    versionInfo_.load(versionPath);
    if (versionInfo_.indexVersion() == 0) {
      fmt::MemoryWriter infostr;
      infostr << "Error: The index version file " << versionPath.string()
              << " doesn't seem to exist.  Please try re-building the salmon "
                 "index.";
      throw std::invalid_argument(infostr.str());
    }
    // Check index version compatibility here

    auto indexType = versionInfo_.indexType();
    // Load the appropriate index type
    if (indexType == SalmonIndexType::FMD) {
      fmt::MemoryWriter infostr;
      infostr << "Error: This version of salmon does not support FMD indexing.";
      throw std::invalid_argument(infostr.str());
    } else if (indexType == SalmonIndexType::QUASI) {
      fmt::MemoryWriter infostr;
      infostr << "Error: This version of salmon does not support indexing using the RapMap index.";
      throw std::invalid_argument(infostr.str());
    } else if (indexType == SalmonIndexType::PUFF) {
      loadPuffIndex_(indexDir);
    } else {
      fmt::MemoryWriter infostr;
      infostr << "Error: Unknown index type.";
      throw std::invalid_argument(infostr.str());
    }

    loaded_ = true;
  }

  bool build(boost::filesystem::path indexDir, pufferfish::IndexOptions& idxOpt) {
    namespace bfs = boost::filesystem;
    switch (versionInfo_.indexType()) {
    case SalmonIndexType::QUASI:
      logger_->error("This version of salmon does not support a RapMap-based index.");
      return false;
    case SalmonIndexType::PUFF:
      return buildPuffIndex_(indexDir, idxOpt);
    case SalmonIndexType::FMD:
      logger_->error("This version of salmon does not support FMD indexing.");
      return false;
    default:
      logger_->warn("Unexpected index type; cannot build");
      return false;
    }
  }

  bool loaded() const { return loaded_; }
  bool isSparse() const { return sparse_; }
  bool is64BitQuasi() const { return largeQuasi_; }
  bool isPerfectHashQuasi() const { return perfectHashQuasi_; }

  PufferfishIndex* puffIndex() { return pfi_.get(); }
  PufferfishSparseIndex* puffSparseIndex() { return pfi_sparse_.get(); }

  SalmonIndexType indexType() { return versionInfo_.indexType(); }

  std::string seqHash256() const { return seqHash256_; }
  std::string nameHash256() const { return nameHash256_; }
  std::string seqHash512() const { return seqHash512_; }
  std::string nameHash512() const { return nameHash512_; }
  std::string decoySeqHash256() const { return decoySeqHash256_; }
  std::string decoyNameHash256() const { return decoyNameHash256_; }

  salmon::utils::DuplicateTargetStatus index_retains_duplicates() const { 
    return keep_duplicates_; 
  }

private:
  bool buildPuffIndex_(boost::filesystem::path indexDir, pufferfish::IndexOptions& idxOpt) {
    namespace bfs = boost::filesystem;
    std::cerr << "out : " << idxOpt.outdir << "\n";
    int ret = pufferfishIndex(idxOpt);
    bfs::path versionFile = indexDir / "versionInfo.json";
    versionInfo_.indexVersion(salmon::indexVersion);
    versionInfo_.hasAuxKmerIndex(false);
    versionInfo_.auxKmerLength(idxOpt.k);
    versionInfo_.indexType(SalmonIndexType::PUFF);
    versionInfo_.salmonVersion(salmon::version);
    versionInfo_.save(versionFile);
    return ret;
  }

  bool loadPuffIndex_(const boost::filesystem::path& indexDir) {
    namespace bfs = boost::filesystem;
    logger_->info("Loading pufferfish index");
    // Read the actual Quasi index
    { // quasi-based
      boost::filesystem::path indexPath = indexDir;
      std::string indexStr = indexDir.string();
      if (indexStr.back() != '/') {
        indexStr.push_back('/');
      }

      std::string sampling_type_;
      int version_;
      std::ifstream indexStream(indexStr + "info.json");
      if (indexStream.is_open()) {
        nlohmann::json info_arch;
        indexStream >> info_arch;

        sampling_type_ = info_arch["sampling_type"];
        seqHash256_ = info_arch["SeqHash"];
        nameHash256_ = info_arch["NameHash"];
        seqHash512_ = info_arch["SeqHash512"];
        nameHash512_ = info_arch["NameHash512"];
        decoySeqHash256_ = info_arch["DecoySeqHash"];
        decoyNameHash256_ = info_arch["DecoyNameHash"];
        version_ = info_arch["index_version"];
        (void) version_;
        if (info_arch.find("keep_duplicates") != info_arch.end()) {
          bool kd = info_arch["keep_duplicates"];
          keep_duplicates_ = kd ? salmon::utils::DuplicateTargetStatus::RETAINED_DUPLICATES : 
            salmon::utils::DuplicateTargetStatus::REMOVED_DUPLICATES;
        } else {
          logger_->warn("The index did not record if the `--keepDuplicates` flag was used. "
          "Please consider re-indexing with a newer version of salmon that will propagate this information.");
        }
      } else {
        logger_->critical(
            "Could not properly open the info.json file from the index : {}.",
            indexStr + "info.json");
      }
      indexStream.close();
      /*
      cereal::JSONInputArchive ar(indexStream);
      ar(cereal::make_nvp("sampling_type", sampling_type_),
         cereal::make_nvp("SeqHash", seqHash256_),
         cereal::make_nvp("NameHash", nameHash256_),
         cereal::make_nvp("SeqHash512", seqHash512_),
         cereal::make_nvp("NameHash512", nameHash512_),
         cereal::make_nvp("DecoySeqHash", decoySeqHash256_),
         cereal::make_nvp("DecoyNameHash", decoyNameHash256_),
         cereal::make_nvp("index_version", version_));
      }
      indexStream.close();
      */

      /*
      if (h.version() != salmon::requiredQuasiIndexVersion) {
        logger_->critical(
            "I found an index with version {}, but I require {}. "
            "Please re-index the reference.",
            h.version(), salmon::requiredQuasiIndexVersion);
        std::exit(1);
      }
      if (h.indexType() != SalmonIndexType::PUFF) {
        logger_->critical("The index {} does not appear to be of the "
                          "appropriate type (pufferfish)",
                          indexStr);
        std::exit(1);
      }

      seqHash256_ = h.seqHash256();
      nameHash256_ = h.nameHash256();
      seqHash512_ = h.seqHash512();
      nameHash512_ = h.nameHash512();
      decoySeqHash256_ = h.decoySeqHash256();
      decoyNameHash256_ = h.decoyNameHash256();
      */

      // Is the quasi-index using a perfect hash
      // perfectHashQuasi_ = h.perfectHash();
      sparse_ = (sampling_type_ == "sparse");

      if (sparse_) {
        logger_->info("Loading sparse pufferfish index.");
        pfi_sparse_.reset(new PufferfishSparseIndex(indexStr));
      } else {
        logger_->info("Loading dense pufferfish index.");
        pfi_.reset(new PufferfishIndex(indexStr));
      }
    }
    logger_->info("done");
    return true;
  }

  /*
  bool isDecoy(uint64_t tid){
    bool decoy{false};
    if (perfectHashQuasi_) {
      decoy = quasiIndexPerfectHash32_->isDecoy(tid);
    } else {
      decoy = quasiIndex32_->isDecoy(tid);
    }
    return decoy;
  }
  */

  bool loaded_;
  SalmonIndexVersionInfo versionInfo_;
  // Can't think of a generally better way to do this now
  // without making the entire code-base look crazy
  bool largeQuasi_{false};
  bool perfectHashQuasi_{false};

  salmon::utils::DuplicateTargetStatus keep_duplicates_{salmon::utils::DuplicateTargetStatus::UNKNOWN};
  bool sparse_{false};
  std::unique_ptr<PufferfishIndex> pfi_{nullptr};
  std::unique_ptr<PufferfishSparseIndex> pfi_sparse_{nullptr};

  std::shared_ptr<spdlog::logger> logger_;
  std::string seqHash256_;
  std::string nameHash256_;
  std::string seqHash512_;
  std::string nameHash512_;
  std::string decoySeqHash256_;
  std::string decoyNameHash256_;
};

#endif //__SALMON_INDEX_HPP
