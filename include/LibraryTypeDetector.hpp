#ifndef __LIBRARY_TYPE_DETECTOR__
#define __LIBRARY_TYPE_DETECTOR__

#include "spdlog/fmt/fmt.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/spdlog.h"

#include "LibraryFormat.hpp"

#include <atomic>
#include <mutex>

class LibraryTypeDetector {
public:
  LibraryTypeDetector(ReadType type)
      : type_(type), libTypeCounts_(std::vector<std::atomic<uint64_t>>(
                         LibraryFormat::maxLibTypeID() + 1)) {}

  LibraryTypeDetector(const LibraryTypeDetector& other) {
    active_.store(other.active_.load());
    type_ = other.type_;
    numSamplesNeeded_.store(other.numSamplesNeeded_.load());
    libTypeCounts_ =
        std::vector<std::atomic<uint64_t>>(LibraryFormat::maxLibTypeID() + 1);
    for (size_t i = 0; i < other.libTypeCounts_.size(); ++i) {
      libTypeCounts_[i] = other.libTypeCounts_[i].load();
    }
  }

  bool isActive() { return active_.load(); }
  bool canGuess() { return numSamplesNeeded_ <= 0; }

  bool mostLikelyType(LibraryFormat& ifmt) {
    bool ret{false};
    if (mut_.try_lock()) {
      if (active_.load()) {
        if (type_ == ReadType::SINGLE_END) {
          uint64_t nf{0};
          uint64_t nr{0};
          for (size_t i = 0; i <= LibraryFormat::maxLibTypeID(); ++i) {
            auto f = LibraryFormat::formatFromID(i);
            auto c = libTypeCounts_[i].load();
            nf += (f.strandedness == ReadStrandedness::S) ? c : 0;
            nr += (f.strandedness == ReadStrandedness::A) ? c : 0;
          }
          double ratio = -1.0;
          if (nf + nr > 0) {
            ratio = static_cast<double>(nf) / (nf + nr);
          }

          ifmt.type = type_;
          ifmt.orientation = ReadOrientation::NONE;

          // If we have some issue computing this, leave as unstranded
          if (ratio < 0.0) {
            ifmt.strandedness = ReadStrandedness::U;
          } else if (ratio < 0.3) {
            // If we map to the forward strand < 30% of the time, we are
            // antisesnse
            ifmt.strandedness = ReadStrandedness::A;
          } else if (ratio < 0.7) {
            // Between 30% and 70% of the time, we are unstranded
            ifmt.strandedness = ReadStrandedness::U;
          } else {
            // Greater than 70% of the time, we are sense
            ifmt.strandedness = ReadStrandedness::S;
          }
        } else { // paired end
          uint64_t nsf{0};
          uint64_t nsr{0};

          uint64_t ninward{0};
          uint64_t noutward{0};
          uint64_t nsame{0};
          for (size_t i = 0; i <= LibraryFormat::maxLibTypeID(); ++i) {
            auto f = LibraryFormat::formatFromID(i);
            auto c = libTypeCounts_[i].load();
            nsf += (f.strandedness == ReadStrandedness::S or
                    f.strandedness == ReadStrandedness::SA)
                       ? c
                       : 0;
            nsr += (f.strandedness == ReadStrandedness::A or
                    f.strandedness == ReadStrandedness::AS)
                       ? c
                       : 0;
            ninward += (f.orientation == ReadOrientation::TOWARD) ? c : 0;
            noutward += (f.orientation == ReadOrientation::AWAY) ? c : 0;
            nsame += (f.orientation == ReadOrientation::SAME) ? c : 0;
          }

          ifmt.type = type_;
          if ((ninward + noutward + nsame > 0) and (nsf + nsr > 0)) {
            auto numOrient = ninward + noutward + nsame;
            double ratioIn = static_cast<double>(ninward) / numOrient;
            double ratioOut = static_cast<double>(noutward) / numOrient;
            double ratioSame = static_cast<double>(nsame) / numOrient;

            ifmt.orientation = ReadOrientation::NONE;
            bool same{false};

            if (ratioIn >= ratioOut and ratioIn >= ratioSame) {
              ifmt.orientation = ReadOrientation::TOWARD;
            } else if (ratioOut >= ratioIn and ratioOut >= ratioSame) {
              ifmt.orientation = ReadOrientation::AWAY;
            } else {
              ifmt.orientation = ReadOrientation::SAME;
              same = true;
            }

            auto numStrand = nsf + nsr;
            double ratioFW = static_cast<double>(nsf) / numStrand;
            if (ratioFW < 0.3) {
              ifmt.strandedness =
                  (same) ? ReadStrandedness::A : ReadStrandedness::AS;
            } else if (ratioFW < 0.7) {
              ifmt.strandedness = ReadStrandedness::U;
            } else {
              ifmt.strandedness =
                  (same) ? ReadStrandedness::S : ReadStrandedness::SA;
            }

          } else {
            ifmt.orientation = ReadOrientation::TOWARD;
            ifmt.strandedness = ReadStrandedness::U;
          }
        } // end paired-end

        auto log = spdlog::get("jointLog");
        log->info("Automatically detected most likely library type as {}\n",
                  ifmt.toString());

        active_.store(false);
        ret = true;
      } // end if active_

      mut_.unlock(); // release the lock
    }                // end try_lock()
    return ret;
  }

  void addSample(LibraryFormat f) {
    if (f.type == type_ and numSamplesNeeded_ >= 0) {
      ++libTypeCounts_[f.formatID()];
      --numSamplesNeeded_;
    }
  }

private:
  // set to false once we have guessed the type
  std::atomic_bool active_{true};
  std::mutex mut_;

  // single or paired-end
  ReadType type_;
  // number of samples needed before we can guess a type
  std::atomic<int64_t> numSamplesNeeded_{50000};

  // the counts for each library type
  std::vector<std::atomic<uint64_t>> libTypeCounts_;
};

#endif //__LIBRARY_TYPE_DETECTOR__
