#ifndef __SALMON_INDEX_HPP__
#define __SALMON_INDEX_HPP__

#include <memory>

#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

#include "cereal/archives/json.hpp"
#include "cereal/types/vector.hpp"
#include "spdlog/spdlog.h"

#include "BooMap.hpp"
#include "FrugalBooMap.hpp"
#include "IndexHeader.hpp"
#include "RapMapSAIndex.hpp"
#include "SalmonConfig.hpp"
#include "SalmonIndexVersionInfo.hpp"

// declaration of quasi index function
int rapMapSAIndex(int argc, char* argv[]);

template <typename IndexT>
using DenseHash =
    spp::sparse_hash_map<uint64_t, rapmap::utils::SAInterval<IndexT>,
                         rapmap::utils::KmerKeyHasher>;
template <typename IndexT>
using PerfectHash = FrugalBooMap<uint64_t, rapmap::utils::SAInterval<IndexT>>;

class SalmonIndex {
public:
  SalmonIndex(std::shared_ptr<spdlog::logger>& logger,
              SalmonIndexType indexType)
      : loaded_(false), versionInfo_(0, false, 0, indexType), logger_(logger),
        seqHash256_(""), nameHash256_(""), seqHash512_(""), nameHash512_("") {}

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
    } else {
      loadQuasiIndex_(indexDir);
    }

    loaded_ = true;
  }

  bool build(boost::filesystem::path indexDir, std::vector<std::string>& argVec,
             uint32_t k) {
    namespace bfs = boost::filesystem;
    switch (versionInfo_.indexType()) {
    case SalmonIndexType::QUASI:
      return buildQuasiIndex_(indexDir, argVec, k);
    case SalmonIndexType::FMD:
      logger_->error("This version of salmon does not support FMD indexing.");
      return false;
    default:
      logger_->warn("Unexpected index type; cannot build");
      return false;
    }
  }

  bool loaded() { return loaded_; }

  bool is64BitQuasi() { return largeQuasi_; }
  bool isPerfectHashQuasi() { return perfectHashQuasi_; }

  RapMapSAIndex<int32_t, DenseHash<int32_t>>* quasiIndex32() {
    return quasiIndex32_.get();
  }
  RapMapSAIndex<int64_t, DenseHash<int64_t>>* quasiIndex64() {
    return quasiIndex64_.get();
  }

  RapMapSAIndex<int32_t, PerfectHash<int32_t>>* quasiIndexPerfectHash32() {
    return quasiIndexPerfectHash32_.get();
  }
  RapMapSAIndex<int64_t, PerfectHash<int64_t>>* quasiIndexPerfectHash64() {
    return quasiIndexPerfectHash64_.get();
  }

  SalmonIndexType indexType() { return versionInfo_.indexType(); }

  const char* transcriptomeSeq() {
    if (loaded_) {
      if (is64BitQuasi()) {
        return quasiIndex64_->seq.c_str();
      } else {
        return quasiIndex32_->seq.c_str();
      }
    } else {
      return nullptr;
    }
  }

  uint64_t transcriptOffset(uint64_t id) {
    if (loaded_) {
      if (is64BitQuasi()) {
        return quasiIndex64_->txpOffsets[id];
      } else {
        return quasiIndex32_->txpOffsets[id];
      }
    } else {
      return std::numeric_limits<uint64_t>::max();
    }
  }

  std::string seqHash256() const { return seqHash256_; }
  std::string nameHash256() const { return nameHash256_; }
  std::string seqHash512() const { return seqHash512_; }
  std::string nameHash512() const { return nameHash512_; }

private:
  bool buildQuasiIndex_(boost::filesystem::path indexDir,
                        std::vector<std::string>& quasiArgVec, uint32_t k) {
    namespace bfs = boost::filesystem;
    int32_t quasiArgc = static_cast<int32_t>(quasiArgVec.size());
    char** quasiArgv = new char*[quasiArgc];
    for (int32_t i = 0; i < quasiArgc; ++i) {
      auto& arg = quasiArgVec[i];
      quasiArgv[i] = new char[arg.size() + 1];
      std::strcpy(quasiArgv[i], arg.c_str());
    }

    int ret = rapMapSAIndex(quasiArgc, quasiArgv);

    bfs::path versionFile = indexDir / "versionInfo.json";
    versionInfo_.indexVersion(salmon::indexVersion);
    versionInfo_.hasAuxKmerIndex(false);
    versionInfo_.auxKmerLength(k);
    versionInfo_.indexType(SalmonIndexType::QUASI);
    versionInfo_.save(versionFile);

    // Free the memory used for the arg vector
    for (int32_t i = 0; i < quasiArgc; ++i) {
      // I hate manual memory management
      delete[] quasiArgv[i];
    }
    delete[] quasiArgv;

    return (ret == 0);
  }

  bool loadQuasiIndex_(const boost::filesystem::path& indexDir) {
    namespace bfs = boost::filesystem;
    logger_->info("Loading Quasi index");
    // Read the actual Quasi index
    { // quasi-based
      boost::filesystem::path indexPath = indexDir;
      std::string indexStr = indexDir.string();
      if (indexStr.back() != '/') {
        indexStr.push_back('/');
      }

      IndexHeader h;
      std::ifstream indexStream(indexStr + "header.json");
      {
        cereal::JSONInputArchive ar(indexStream);
        ar(h);
      }
      indexStream.close();

      if (h.version() != salmon::requiredQuasiIndexVersion) {
        logger_->critical(
            "I found a quasi-index with version {}, but I require {}. "
            "Please re-index the reference.",
            h.version(), salmon::requiredQuasiIndexVersion);
        std::exit(1);
      }
      if (h.indexType() != IndexType::QUASI) {
        logger_->critical("The index {} does not appear to be of the "
                          "appropriate type (quasi)",
                          indexStr);
        std::exit(1);
      }

      seqHash256_ = h.seqHash256();
      nameHash256_ = h.nameHash256();
      seqHash512_ = h.seqHash512();
      nameHash512_ = h.nameHash512();

      // Is the quasi-index using a perfect hash
      perfectHashQuasi_ = h.perfectHash();

      if (h.bigSA()) {
        largeQuasi_ = true;
        logger_->info("Loading 64-bit quasi index");
        if (perfectHashQuasi_) {
          quasiIndexPerfectHash64_.reset(
              new RapMapSAIndex<int64_t, PerfectHash<int64_t>>);
          if (!quasiIndexPerfectHash64_->load(indexStr)) {
            fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
            fmt::print(stderr, "Please make sure that 'salmon index' has been "
                               "run successfully\n");
            std::exit(1);
          }
        } else {
          quasiIndex64_.reset(new RapMapSAIndex<int64_t, DenseHash<int64_t>>);
          if (!quasiIndex64_->load(indexStr)) {
            fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
            fmt::print(stderr, "Please make sure that 'salmon index' has been "
                               "run successfully\n");
            std::exit(1);
          }
        }
      } else { // 32-bit index
        logger_->info("Loading 32-bit quasi index");

        if (perfectHashQuasi_) {
          quasiIndexPerfectHash32_.reset(
              new RapMapSAIndex<int32_t, PerfectHash<int32_t>>);
          if (!quasiIndexPerfectHash32_->load(indexStr)) {
            fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
            fmt::print(stderr, "Please make sure that 'salmon index' has been "
                               "run successfully\n");
            std::exit(1);
          }
        } else {
          quasiIndex32_.reset(new RapMapSAIndex<int32_t, DenseHash<int32_t>>);
          if (!quasiIndex32_->load(indexStr)) {
            fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
            fmt::print(stderr, "Please make sure that 'salmon index' has been "
                               "run successfully\n");
            std::exit(1);
          }
        }
      }
    }
    logger_->info("done");
    return true;
  }

  bool isDecoy(uint64_t tid){
    bool decoy{false};
    if (largeQuasi_) {
      if (perfectHashQuasi_) {
        decoy = quasiIndexPerfectHash64_->isDecoy(tid);
      } else {
        decoy = quasiIndex64_->isDecoy(tid);
      }
    } else {
      if (perfectHashQuasi_) {
        decoy = quasiIndexPerfectHash32_->isDecoy(tid);
      } else {
        decoy = quasiIndex32_->isDecoy(tid);
      }
    }
    // should not get here
    return decoy;
  }

  bool loaded_;
  SalmonIndexVersionInfo versionInfo_;
  // Can't think of a generally better way to do this now
  // without making the entire code-base look crazy
  bool largeQuasi_{false};
  bool perfectHashQuasi_{false};

  std::unique_ptr<RapMapSAIndex<int32_t, DenseHash<int32_t>>> quasiIndex32_{
      nullptr};
  std::unique_ptr<RapMapSAIndex<int64_t, DenseHash<int64_t>>> quasiIndex64_{
      nullptr};

  std::unique_ptr<RapMapSAIndex<int32_t, PerfectHash<int32_t>>>
      quasiIndexPerfectHash32_{nullptr};
  std::unique_ptr<RapMapSAIndex<int64_t, PerfectHash<int64_t>>>
      quasiIndexPerfectHash64_{nullptr};

  std::shared_ptr<spdlog::logger> logger_;
  std::string seqHash256_;
  std::string nameHash256_;
  std::string seqHash512_;
  std::string nameHash512_;
};

#endif //__SALMON_INDEX_HPP
