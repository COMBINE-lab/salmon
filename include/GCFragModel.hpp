#ifndef __GC_FRAG_MODEL__
#define __GC_FRAG_MODEL__

#include "DistributionUtils.hpp"
#include "Eigen/Dense"
#include "SalmonMath.hpp"

#include <boost/iostreams/filtering_stream.hpp>

#include <iostream>
#include <vector>

struct GCDesc {
  int32_t fragFrac;
  int32_t contextFrac;

  // assumes 101 bins
  int32_t fragBin() { return fragFrac; }
  int32_t contextBin() { return contextFrac; }

  int32_t fragBin(int32_t n) {
    double w = (100.0 / n);
    return std::min(n - 1, static_cast<int32_t>(fragFrac / w));
  }
  int32_t contextBin(int32_t n) {
    double w = (100.0 / n);
    return std::min(n - 1, static_cast<int32_t>(contextFrac / w));
  }

  bool operator==(const GCDesc& other) const {
    return fragFrac == other.fragFrac and contextFrac == other.contextFrac;
  }
  bool operator==(const GCDesc&& other) const {
    return fragFrac == other.fragFrac and contextFrac == other.contextFrac;
  }

  friend std::ostream& operator<<(std::ostream& os, const GCDesc& c);
};

inline std::ostream& operator<<(std::ostream& os, const GCDesc& c) {
    os << "{ fragFrac : " << c.fragFrac << ", contextFrac : " << c.contextFrac << "}\n";
    return os;
}



class GCFragModel {
public:
  GCFragModel(size_t condBins = 3, size_t numGCBins = 101,
              distribution_utils::DistributionSpace dspace =
                  distribution_utils::DistributionSpace::LOG)
      : condBins_(condBins), numGCBins_(numGCBins), dspace_(dspace),
        normalized_(false) {
    counts_ = Eigen::MatrixXd(condBins_, numGCBins_);
    if (dspace_ == distribution_utils::DistributionSpace::LOG) {
      counts_.setOnes();
      counts_ *= salmon::math::LOG_0;
    } else {
      counts_.setZero();
    }
    // set the total vector to be the right size and full of 0's.
    modelTotals_.resize(condBins_, 0.0);
  }

  bool writeBinary(boost::iostreams::filtering_ostream& out) const {
    auto* mutThis = const_cast<GCFragModel*>(this);
    int32_t dtype =
        (dspace_ == distribution_utils::DistributionSpace::LINEAR) ? 0 : 1;
    out.write(reinterpret_cast<char*>(&dtype), sizeof(dtype));
    typename Eigen::MatrixXd::Index rows = counts_.rows(),
                                    cols = counts_.cols();
    out.write(reinterpret_cast<char*>(&rows),
              sizeof(typename Eigen::MatrixXd::Index));
    out.write(reinterpret_cast<char*>(&cols),
              sizeof(typename Eigen::MatrixXd::Index));
    out.write(reinterpret_cast<char*>(const_cast<double*>(modelTotals_.data())),
              sizeof(double) * rows);
    out.write(reinterpret_cast<char*>(mutThis->counts_.data()),
              rows * cols * sizeof(typename Eigen::MatrixXd::Scalar));
    return true;
  }

  GCFragModel(const GCFragModel&) = default;
  GCFragModel(GCFragModel&&) = default;
  GCFragModel& operator=(const GCFragModel&) = default;
  GCFragModel& operator=(GCFragModel&&) = default;

  /*
  double likelihood_(uint32_t numBins) {
  }
  uint32_t optNumBins() {
    if (numGCBins_ != 101) {
      std::cerr << "Selecting the optimal number of bins is currently "
                << "only supported when the initial histograms are generated "
                << "using 101 bins.\n";
      std::exit(1);
    }
    std::vector<uint32_t> nbins{5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 101};
    // Using the "BR" penalty
    // see:
    // Davies, Laurie, et al. "A comparison of automatic histogram
  constructions."
    // ESAIM: Probability and Statistics 13 (2009): 181-196.
    std::vector<double> scores;
    for (auto nb : nbins) {
      likelihood_(nb)
    }
  }
  */

  void reset(distribution_utils::DistributionSpace dspace =
                 distribution_utils::DistributionSpace::LOG) {
    normalized_ = false;
    dspace_ = dspace;
    if (dspace_ == distribution_utils::DistributionSpace::LOG) {
      counts_.setOnes();
      counts_ *= salmon::math::LOG_0;
    } else {
      counts_.setZero();
    }
  }

  GCFragModel ratio(GCFragModel& other, double maxRatio) {
    if (!normalized_) {
      normalize();
    }
    if (!other.normalized_) {
      other.normalize();
    }
    double minRatio = 1.0 / maxRatio;

    GCFragModel ratioModel(condBins_, numGCBins_, dspace_);
    for (size_t r = 0; r < condBins_; ++r) {
      for (size_t c = 0; c < numGCBins_; ++c) {
        double rat = (counts_(r, c) / other.counts_(r, c));
        if (rat > maxRatio) {
          rat = maxRatio;
        }
        if (rat < minRatio) {
          rat = minRatio;
        }
        ratioModel.counts_(r, c) = rat;
      }
    }
    return ratioModel;
  }

  void inc(GCDesc desc,
           double fragWeight //< the weight associated with this fragment
  ) {
    auto ctx = (condBins_ > 1) ? desc.contextBin(condBins_) : 0;
    auto frag = (numGCBins_ != 101) ? desc.fragBin(numGCBins_) : desc.fragBin();

    if (dspace_ == distribution_utils::DistributionSpace::LOG) {
      counts_(ctx, frag) = salmon::math::logAdd(counts_(ctx, frag), fragWeight);
    } else {
      counts_(ctx, frag) += fragWeight;
    }
  }

  double get(GCDesc desc) {
    auto ctx = (condBins_ > 1) ? desc.contextBin(condBins_) : 0;
    auto frag = (numGCBins_ != 101) ? desc.fragBin(numGCBins_) : desc.fragBin();
    return counts_(ctx, frag);
  }

  distribution_utils::DistributionSpace distributionSpace() const {
    return dspace_;
  }

  void combineCounts(const GCFragModel& other) {
    if (dspace_ != other.dspace_) {
      std::cerr
          << "Cannot combine distributions that live in a different space!\n";
      std::exit(1);
    }
    if (dspace_ == distribution_utils::DistributionSpace::LOG) {
      for (size_t r = 0; r < condBins_; ++r) {
        for (size_t c = 0; c < numGCBins_; ++c) {
          counts_(r, c) =
              salmon::math::logAdd(counts_(r, c), other.counts_(r, c));
        }
      }
    } else {
      for (size_t r = 0; r < condBins_; ++r) {
        for (size_t c = 0; c < numGCBins_; ++c) {
          counts_(r, c) += other.counts_(r, c);
        }
      }
    }
  }

  /**
   * NOTE: Improve interface --- also converts out of log space
   */
  void normalize(double prior = 0.1) {
    if (!normalized_) {
      if (dspace_ == distribution_utils::DistributionSpace::LOG) {
        prior = std::log(prior);
        for (size_t r = 0; r < condBins_; ++r) {
          double rowMass{salmon::math::LOG_0};
          for (size_t c = 0; c < numGCBins_; ++c) {
            rowMass = salmon::math::logAdd(
                prior, salmon::math::logAdd(rowMass, counts_(r, c)));
          }
          if (!salmon::math::isLog0(rowMass)) {
            for (size_t c = 0; c < numGCBins_; ++c) {
              counts_(r, c) = std::exp(
                  salmon::math::logAdd(prior, counts_(r, c)) - rowMass);
            }
            modelTotals_[r] = std::exp(rowMass);
          }
          // if rowMass is LOG_0, then leave modelTotals_[r] as 0.0
        }
      } else {
        for (size_t r = 0; r < condBins_; ++r) {
          double rowMass = 0.0;
          for (size_t c = 0; c < numGCBins_; ++c) {
            rowMass += (prior + counts_(r, c));
          }
          if (rowMass > 0.0) {
            double norm = 1.0 / rowMass;
            for (size_t c = 0; c < numGCBins_; ++c) {
              counts_(r, c) = (prior + counts_(r, c)) * norm;
            }
            modelTotals_[r] = rowMass;
          }
          // if rowMass is 0.0, just leave modelTotals_[r] as 0.0;
        }
      }
      normalized_ = true;
      dspace_ = distribution_utils::DistributionSpace::LINEAR;
    }
  }

private:
  size_t condBins_;
  size_t numGCBins_;
  distribution_utils::DistributionSpace dspace_;
  bool normalized_;
  Eigen::MatrixXd counts_;
  std::vector<double> modelTotals_;
};

#endif //__GC_FRAG_MODEL__
