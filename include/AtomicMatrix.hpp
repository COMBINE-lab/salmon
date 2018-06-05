#ifndef ATOMIC_MATRIX
#define ATOMIC_MATRIX

#include "tbb/atomic.h"
#include "tbb/concurrent_vector.h"

#include "SalmonMath.hpp"

#include <vector>

template <typename T> class AtomicMatrix {
public:
  AtomicMatrix(size_t nRow, size_t nCol, T alpha, bool logSpace = true)
      : storage_(nRow * nCol, logSpace ? std::log(alpha) : alpha),
        rowsums_(nRow, logSpace ? std::log(nCol * alpha) : nCol * alpha),
        nRow_(nRow), nCol_(nCol), alpha_(alpha), logSpace_(logSpace) {}

  void incrementUnnormalized(size_t rowInd, size_t colInd, T amt) {
    using salmon::math::logAdd;
    size_t k = rowInd * nCol_ + colInd;
    if (logSpace_) {
      T oldVal = storage_[k];
      T retVal = oldVal;
      T newVal = logAdd(oldVal, amt);
      do {
        oldVal = retVal;
        newVal = logAdd(oldVal, amt);
        retVal = storage_[k].compare_and_swap(newVal, oldVal);
      } while (retVal != oldVal);

    } else {
      T oldVal = storage_[k];
      T retVal = oldVal;
      T newVal = oldVal + amt;
      do {
        oldVal = retVal;
        newVal = oldVal + amt;
        retVal = storage_[k].compare_and_swap(newVal, oldVal);
      } while (retVal != oldVal);
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
      T oldVal = storage_[k];
      T retVal = oldVal;
      T newVal = logAdd(oldVal, amt);
      do {
        oldVal = retVal;
        newVal = logAdd(oldVal, amt);
        retVal = storage_[k].compare_and_swap(newVal, oldVal);
      } while (retVal != oldVal);

      oldVal = rowsums_[rowInd];
      retVal = oldVal;
      newVal = logAdd(oldVal, amt);
      do {
        oldVal = retVal;
        newVal = logAdd(oldVal, amt);
        retVal = rowsums_[rowInd].compare_and_swap(newVal, oldVal);
      } while (retVal != oldVal);
    } else {
      T oldVal = storage_[k];
      T retVal = oldVal;
      T newVal = oldVal + amt;
      do {
        oldVal = retVal;
        newVal = oldVal + amt;
        retVal = storage_[k].compare_and_swap(newVal, oldVal);
      } while (retVal != oldVal);

      oldVal = rowsums_[rowInd];
      retVal = oldVal;
      newVal = oldVal + amt;
      do {
        oldVal = retVal;
        newVal = oldVal + amt;
        retVal = rowsums_[rowInd].compare_and_swap(newVal, oldVal);
      } while (retVal != oldVal);
    }
  }

  T operator()(size_t rowInd, size_t colInd, bool normalized = true) {
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
  std::vector<tbb::atomic<T>> storage_;
  std::vector<tbb::atomic<T>> rowsums_;
  size_t nRow_, nCol_;
  T alpha_;
  bool logSpace_;
};

#endif // ATOMIC_MATRIX
