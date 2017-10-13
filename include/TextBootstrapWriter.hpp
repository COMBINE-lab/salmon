#ifndef __TEXT_BOOTSTRAP_WRITER_HPP__
#define __TEXT_BOOTSTRAP_WRITER_HPP__

#include <fstream>
#include <memory>
#include <mutex>
#include <vector>

#include "BootstrapWriter.hpp"
#include "SalmonSpinLock.hpp"
#include "Transcript.hpp"
#include "spdlog/spdlog.h"

class TextBootstrapWriter : public BootstrapWriter {
public:
  TextBootstrapWriter(boost::filesystem::path& outputPath,
                      std::shared_ptr<spdlog::logger> logger)
      : outputPath_(outputPath), logger_(logger) {
    // Create the directory if it doesn't exist
    if (!boost::filesystem::exists(outputPath_.parent_path())) {
      boost::filesystem::create_directories(outputPath_);
    }
    // open the file
    ofile_.open(outputPath.string());
  }

  ~TextBootstrapWriter() {
#if defined __APPLE__
    spin_lock::scoped_lock sl(writeMutex_);
#else
    std::lock_guard<std::mutex> lock(writeMutex_);
#endif
    ofile_.close();
  }

  bool writeHeader(std::string& comments,
                   std::vector<Transcript>& transcripts) override {
#if defined __APPLE__
    spin_lock::scoped_lock sl(writeMutex_);
#else
    std::lock_guard<std::mutex> lock(writeMutex_);
#endif
    ofile_ << comments;
    size_t numTxps = transcripts.size();
    if (numTxps == 0) {
      return false;
    }
    for (size_t tn = 0; tn < numTxps; ++tn) {
      auto& t = transcripts[tn];
      ofile_ << t.RefName;
      if (tn < numTxps - 1) {
        ofile_ << '\t';
      }
    }
    ofile_ << '\n';
    /*
    for (size_t tn = 0; tn < numTxps; ++tn) {
        auto& t  = transcripts[tn];
        ofile_ << t.EffectiveLength;
        if (tn < numTxps - 1) {
            ofile_ << '\t';
        }
    }
    ofile_ << '\n';
    */
    return true;
  }

  bool writeBootstrap(std::vector<double>& abund) override {
#if defined __APPLE__
    spin_lock::scoped_lock sl(writeMutex_);
#else
    std::lock_guard<std::mutex> lock(writeMutex_);
#endif
    size_t numTxps = abund.size();
    for (size_t tn = 0; tn < numTxps; ++tn) {
      auto& a = abund[tn];
      ofile_ << a;
      if (tn < numTxps - 1) {
        ofile_ << '\t';
      }
    }
    ofile_ << '\n';
    logger_->info("wrote {} bootstraps", numWritten_.load() + 1);
    ++numWritten_;
    return true;
  }

private:
  boost::filesystem::path outputPath_;
  std::ofstream ofile_;
  std::shared_ptr<spdlog::logger> logger_;
// only one writer thread at a time
#if defined __APPLE__
  spin_lock writeMutex_;
#else
  std::mutex writeMutex_;
#endif
  std::atomic<uint32_t> numWritten_{0};
};

#endif // __TEXT_BOOTSTRAP_WRITER_HPP__
