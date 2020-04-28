#ifndef __GZIP_WRITER_HPP__
#define __GZIP_WRITER_HPP__

#include <memory>
#include <mutex>

#include "spdlog/spdlog.h"

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "ReadExperiment.hpp"
#include "SalmonOpts.hpp"
#include "SalmonSpinLock.hpp"
#include "AlevinOpts.hpp"
#include "Graph.hpp"
#include "MappingStatistics.hpp"

class GZipWriter {
public:
  GZipWriter(const boost::filesystem::path path,
             std::shared_ptr<spdlog::logger> logger);

  ~GZipWriter();

  void close_all_streams();

  template <typename ExpT>
  bool writeEquivCounts(const SalmonOpts& opts, ExpT& experiment);

  template <typename ExpT, typename ProtocolT>
  bool writeEquivCounts(const AlevinOpts<ProtocolT>& aopts,
                        ExpT& experiment);

  template <typename ExpT>
  bool writeBFH(boost::filesystem::path& outDir,
                ExpT& experiment, size_t umiLength,
                std::vector<std::string>& bcSeqVec);

  template <typename ExpT>
  bool writeMeta(const SalmonOpts& opts, const ExpT& experiment, const MappingStatistics& mstats);

  template <typename ProtocolT>
  bool writeMetaAlevin(const AlevinOpts<ProtocolT>& opts,
                       boost::filesystem::path aux_dir);

  template <typename ExpT>
  bool writeEmptyMeta(const SalmonOpts& opts, const ExpT& experiment,
                      std::vector<std::string>& errors);

  template <typename ExpT>
  bool writeAbundances(const SalmonOpts& sopt, ExpT& readExp, bool explicitSum);

  bool writeAbundances(std::vector<double>& alphas,
                       std::vector<Transcript>& transcripts);

  bool writeAbundances(std::string& bcName,
                       std::string& features,
                       uint8_t featureCode,
                       std::vector<double>& alphas,
                       std::vector<uint8_t>& tiers,
                       bool dumpUmiGraph);

  bool writeSparseAbundances(std::string& bcName,
                             std::string& features,
                             std::string& arboData,
                             uint8_t featureCode,
                             std::vector<double>& alphas,
                             std::vector<uint8_t>& tiers,
                             bool dumpArboFragCounts,
                             bool dumpUmiGraph);

  bool writeBootstraps(std::string& bcName,
                       std::vector<double>& alphas,
                       std::vector<double>& variance,
                       bool useAllBootstraps,
                       std::vector<std::vector<double>>& sampleEstimates);

  bool writeSparseBootstraps(std::string& bcName,
                             std::vector<double>& alphas,
                             std::vector<double>& variance,
                             bool useAllBootstraps,
                             std::vector<std::vector<double>>& sampleEstimates);

  template <typename ExpT>
  bool writeEmptyAbundances(const SalmonOpts& sopt, ExpT& readExp);

  template <typename T>
  bool writeBootstrap(const std::vector<T>& abund, bool quiet = false);

  bool writeCellEQVec(size_t barcode, const std::vector<uint32_t>& offsets,
                      const std::vector<uint32_t>& counts, bool quiet = true);

  bool writeDedupCellEQVec(size_t barcode, 
                           const std::vector<std::vector<uint32_t>>& labels,
                           const std::vector<uint32_t>& counts,
                           bool quiet = true);

  bool writeUmiGraph(alevin::graph::Graph& g, std::string& trueBarcodeStr);

  bool setSamplingPath(const SalmonOpts& sopt);

  void writeMtx(std::shared_ptr<spdlog::logger>& jointLog, 
                boost::filesystem::path& outputDirectory,
                size_t numGenes, size_t numCells, size_t totalExpGeneCounts);

private:
  boost::filesystem::path path_;
  boost::filesystem::path bsPath_;
  std::shared_ptr<spdlog::logger> logger_;
  std::unique_ptr<boost::iostreams::filtering_ostream> bsStream_{nullptr};
  std::unique_ptr<boost::iostreams::filtering_ostream> countMatrixStream_{nullptr};
  std::unique_ptr<boost::iostreams::filtering_ostream> meanMatrixStream_{nullptr};
  std::unique_ptr<boost::iostreams::filtering_ostream> varMatrixStream_{nullptr};
  std::unique_ptr<boost::iostreams::filtering_ostream> fullBootstrapMatrixStream_{nullptr};
  std::unique_ptr<boost::iostreams::filtering_ostream> tierMatrixStream_{nullptr};
  std::unique_ptr<boost::iostreams::filtering_ostream> umiGraphStream_{nullptr};
  std::unique_ptr<boost::iostreams::filtering_ostream> cellEQStream_{nullptr};
  std::unique_ptr<boost::iostreams::filtering_ostream> cellDedupEQStream_{nullptr};
  std::unique_ptr<boost::iostreams::filtering_ostream> arboMatrixStream_{nullptr};
  std::unique_ptr<std::ofstream> bcNameStream_{nullptr};
  std::unique_ptr<std::ofstream> bcFeaturesStream_{nullptr};
  std::unique_ptr<std::ofstream> bcBootNameStream_{nullptr};
// only one writer thread at a time
#if defined __APPLE__
  spin_lock writeMutex_;
#else
  std::mutex writeMutex_;
#endif
  std::atomic<uint32_t> numBootstrapsWritten_{0};
};

#endif //__GZIP_WRITER_HPP__
