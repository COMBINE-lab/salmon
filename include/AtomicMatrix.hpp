#ifndef ATOMIC_MATRIX
#define ATOMIC_MATRIX

#include "tbb/concurrent_vector.h"

#include "SalmonMath.hpp"
#include "SalmonUtils.hpp"

#include <vector>

template <typename T> class AtomicMatrix {
public:
  AtomicMatrix() {
    nRow_ = 0;
    nCol_ = 0;
    alpha_ = salmon::math::LOG_0;
    logSpace_ = true;
  }

  AtomicMatrix(size_t nRow, size_t nCol, T alpha, bool logSpace = true)
      : nRow_(nRow), nCol_(nCol), alpha_(alpha), logSpace_(logSpace) {

          decltype(storage_) storage_tmp(nRow * nCol);
          std::swap(storage_, storage_tmp);
          T e = logSpace ? std::log(alpha) : alpha;
          std::fill(storage_.begin(), storage_.end(), e);

          decltype(rowsums_) rowsums_tmp(nRow);
          std::swap(rowsums_, rowsums_tmp);
          T ers = logSpace ? std::log(nCol * alpha) : nCol * alpha;
          std::fill(rowsums_.begin(), rowsums_.end(), ers);
        }

  AtomicMatrix& operator=(AtomicMatrix&& o) {
    std::swap(storage_, o.storage_);
    std::swap(rowsums_, o.rowsums_);
    nRow_ = o.nRow_;
    nCol_ = o.nCol_;
    alpha_ = o.alpha_;
    logSpace_ = o.logSpace_;
    return *this;
  }

  void incrementUnnormalized(size_t rowInd, size_t colInd, T amt) {
    using salmon::math::logAdd;
    size_t k = rowInd * nCol_ + colInd;
    if (logSpace_) {
      salmon::utils::incLoopLog(storage_[k], amt);
    } else {
      salmon::utils::incLoop(storage_[k], amt);
    }
  }

  void computeRowSums() {
    for (size_t rowInd = 0; rowInd < nRow_; ++rowInd) {
      T rowSum = salmon::math::LOG_0;
      for (size_t colInd = 0; colInd < nCol_; ++colInd) {
        size_t k = rowInd * nCol_ + colInd;
        rowSum = logAdd(rowSum, storage_[k]);
      }
      rowsums_[rowInd] = rowSum;
    }
  }

  void increment(size_t rowInd, size_t colInd, T amt) {
    using salmon::math::logAdd;
    size_t k = rowInd * nCol_ + colInd;
    if (logSpace_) {
      salmon::utils::incLoopLog(storage_[k], amt);
      salmon::utils::incLoopLog(rowsums_[rowInd], amt);
    } else {
      salmon::utils::incLoop(storage_[k], amt);
      salmon::utils::incLoop(rowsums_[rowInd], amt);
    }
  }

  T operator()(size_t rowInd, size_t colInd /*, bool normalized = true*/) {
    size_t k = rowInd * nCol_ + colInd;
    if (logSpace_) {
      return storage_[k] - rowsums_[rowInd];
    } else {
      return storage_[k] / rowsums_[rowInd];
    }
  }

  size_t nRow() const { return nRow_; }
  size_t nCol() const { return nCol_; }

private:
  std::vector<std::atomic<T>> storage_;
  std::vector<std::atomic<T>> rowsums_;
  size_t nRow_, nCol_;
  T alpha_;
  bool logSpace_;
};

#endif // ATOMIC_MATRIX
