#ifndef __PCA__HPP__
#define __PCA__HPP__

#include "Eigen/Dense"
#include <iostream>
#include <unordered_set>

class PCA {
public:
  PCA(Eigen::MatrixXd& m) { // Eigen::Map<Eigen::MatrixXd>& m) {
    dat_ = m;
    /*
    dat_ = Eigen::MatrixXd(m.rows(), m.cols());
    for(size_t i = 0; i < m.rows(); ++i) {
        for (size_t j = 0; j < m.cols(); ++j) {
            dat_(i,j) = m(i,j);
        }
    }
    */
  }

  void performDecomposition(bool standardize = true) {
    if (standardize) {
      std::unordered_set<uint32_t> droppedCols;

      size_t nrows = dat_.rows();
      Eigen::VectorXd means(dat_.cols());
      means = dat_.colwise().mean();
      Eigen::VectorXd stdDev(dat_.cols());
      decltype(dat_.cols()) i, j;
      for (i = 0; i < dat_.cols(); ++i) {
        for (j = 0; j < dat_.rows(); ++j) {
          dat_(j, i) -= means(i);
        }
        auto x = dat_.col(i).dot(dat_.col(i));
        auto sd = std::sqrt(x / nrows);
        if (sd < 1e-20) {
          droppedCols.insert(i);
        } else {
          dat_.col(i) /= sd;
        }
      }

      // If we dropped any columns for having a "0" standard deviation,
      // then create a new normalized data matrix without them.
      if (droppedCols.size() > 0) {
        Eigen::MatrixXd tmpDat(dat_.rows(), dat_.cols() - droppedCols.size());
        size_t curCol{0};
        decltype(dat_.cols()) i;
        for (i = 0; i < dat_.cols(); ++i) {
          if (droppedCols.find(i) == droppedCols.end()) {
            tmpDat.col(curCol) = dat_.col(i);
            ++curCol;
          }
        }
        std::swap(dat_, tmpDat);
      }
    }

    // std::cerr << "DATA = " << dat_ << "\n";
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(dat_, Eigen::ComputeThinU |
                                                    Eigen::ComputeThinV);
    proj_ = svd.matrixV();
    variances_ = svd.singularValues();
    // std::cerr << "V = " << proj_ << "\n";
    // std::cerr << "S = " << variances_;
  }

  uint32_t varFrac(double frac, bool verbose = false) {
    uint32_t cutoff{0};
    bool haveCutoff{false};

    if (verbose) {
      std::cerr << "variances: " << variances_ << "\n";
    }
    Eigen::VectorXd fracs(variances_.rows());
    double totalVar{0.0};
    for (uint32_t i = 0; i < variances_.rows(); ++i) {
      variances_(i) = variances_(i) * variances_(i);
      totalVar += variances_(i);
    }
    double partialSum{0.0};
    for (uint32_t i = 0; i < variances_.rows(); ++i) {
      partialSum += variances_(i);
      fracs(i) = partialSum / totalVar;
      if (!haveCutoff and fracs(i) >= frac) {
        cutoff = i;
        haveCutoff = true;
      }
    }
    if (verbose) {
      std::cerr << "fracs = " << fracs << "\n";
      std::cerr << "include the first " << cutoff + 1
                << " components to reach cutoff of " << frac << "\n";
    }
    return cutoff + 1;
  }

  Eigen::MatrixXd projectedData(double frac, bool verbose = false) {
    uint32_t numVec = varFrac(frac, verbose);
    return dat_ * proj_.block(0, 0, proj_.rows(), numVec);
  }

private:
  Eigen::MatrixXd dat_;
  Eigen::MatrixXd proj_;
  Eigen::VectorXd variances_;
};

#endif //__PCA__HPP__
