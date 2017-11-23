#ifndef __TRYABLE_SPIN_LOCK_HPP__
#define __TRYABLE_SPIN_LOCK_HPP__

#include <atomic>

// Taken from
// http://stackoverflow.com/questions/26583433/c11-implementation-of-spinlock-using-atomic
class TryableSpinLock {
  std::atomic_flag locked = ATOMIC_FLAG_INIT;

public:
  void lock() {
    while (locked.test_and_set(std::memory_order_acquire)) {
      ;
    }
  }

  // from
  // http://stackoverflow.com/questions/19742993/implementing-a-spinlock-in-boost-example-neededhttp://stackoverflow.com/questions/19742993/implementing-a-spinlock-in-boost-example-needed
  // is this legit?
  bool try_lock() { return !locked.test_and_set(std::memory_order_acquire); }

  void unlock() { locked.clear(std::memory_order_release); }
};

#endif //__TRYABLE_SPIN_LOCK_HPP__
