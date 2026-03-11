#ifndef __SALMON_SPIN_LOCK_HPP__
#define __SALMON_SPIN_LOCK_HPP__

#include <atomic>

/**
 * since std::mutex is *VERY SLOW* on OSX, we use
 * this there instead.
 * Taken from
 * (http://stackoverflow.com/questions/22899053/why-is-stdmutex-so-slow-on-osx)
 */
class spin_lock {
  std::atomic<bool> _lock;

public:
  spin_lock(const spin_lock&) = delete;
  spin_lock& operator=(const spin_lock&) = delete;

  spin_lock() : _lock(false) {}

  class scoped_lock {
    spin_lock& _lock;

  public:
    scoped_lock(const scoped_lock&) = delete;
    scoped_lock& operator=(const scoped_lock&) = delete;

    scoped_lock(spin_lock& lock) : _lock(lock) {
      bool expect = false;
      while (!_lock._lock.compare_exchange_weak(expect, true)) {
        expect = false;
      }
    }
    ~scoped_lock() { _lock._lock = false; }
  };
};

#endif // __SALMON_SPIN_LOCK_HPP__
