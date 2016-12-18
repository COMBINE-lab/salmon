#ifndef __GZIP_WRITER_HPP__
#define __GZIP_WRITER_HPP__

#include <memory>
#include <mutex>

#include "spdlog/spdlog.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "SalmonSpinLock.hpp"
#include "SalmonOpts.hpp"
#include "ReadExperiment.hpp"

class GZipWriter {
  public:
    GZipWriter(const boost::filesystem::path path, std::shared_ptr<spdlog::logger> logger);

    ~GZipWriter();

    template <typename ExpT>
    bool writeEquivCounts(
	const SalmonOpts& opts,
	ExpT& experiment);

    template <typename ExpT>
    bool writeMeta(
	const SalmonOpts& opts,
	const ExpT& experiment);

    template <typename ExpT>
    bool writeAbundances(
      const SalmonOpts& sopt,
      ExpT& readExp);

    template <typename T>
    bool writeBootstrap(const std::vector<T>& abund, bool quiet=false);

  bool setSamplingPath(const SalmonOpts& sopt);
   private:
     boost::filesystem::path path_;
     boost::filesystem::path bsPath_;
     std::shared_ptr<spdlog::logger> logger_;
     std::unique_ptr<boost::iostreams::filtering_ostream> bsStream_{nullptr};
// only one writer thread at a time
#if defined __APPLE__
        spin_lock writeMutex_;
#else
        std::mutex writeMutex_;
#endif
        std::atomic<uint32_t> numBootstrapsWritten_{0};
};

#endif //__GZIP_WRITER_HPP__
