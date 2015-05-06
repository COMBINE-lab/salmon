/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef __JELLYFISH_LOCKS_PTHREAD_HPP__
#define __JELLYFISH_LOCKS_PTHREAD_HPP__

#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace jellyfish { namespace locks{ namespace pthread {
class cond
{
  pthread_mutex_t   _mutex;
  pthread_cond_t	_cond;

public:
  cond() {
    pthread_mutex_init(&_mutex, NULL);
    pthread_cond_init(&_cond, NULL);
  }

  ~cond() {
    pthread_cond_destroy(&_cond);
    pthread_mutex_destroy(&_mutex);
  }

  inline void lock() { pthread_mutex_lock(&_mutex); }
  inline void unlock() { pthread_mutex_unlock(&_mutex); }
  inline void wait() { pthread_cond_wait(&_cond, &_mutex); }
  inline void signal() { pthread_cond_signal(&_cond); }
  inline void broadcast() { pthread_cond_broadcast(&_cond); }
  inline int timedwait(struct timespec *abstime) {
    return pthread_cond_timedwait(&_cond, &_mutex, abstime);
  }
  inline int timedwait(time_t seconds) {
    struct timespec curtime;
#ifdef HAVE_CLOCK_GETTIME
    clock_gettime(CLOCK_REALTIME, &curtime);
#else
    struct timeval timeofday;
    gettimeofday(&timeofday, 0);
    curtime.tv_sec  = timeofday.tv_sec;
    curtime.tv_nsec = timeofday.tv_usec * 1000;
#endif
    curtime.tv_sec += seconds;
    return timedwait(&curtime);
  }
};

class mutex {
  pthread_mutex_t     _mutex;

public:
  mutex(int type = PTHREAD_MUTEX_DEFAULT) {
    pthread_mutexattr_t attr;
    pthread_mutexattr_init(&attr);
    pthread_mutexattr_settype(&attr, type);
    pthread_mutex_init(&_mutex, &attr);
  }

  ~mutex() {
    pthread_mutex_destroy(&_mutex);
  }

  inline void lock() { pthread_mutex_lock(&_mutex); }
  inline void unlock() { pthread_mutex_unlock(&_mutex); }
  inline bool try_lock() { return !pthread_mutex_trylock(&_mutex); }
};

class mutex_recursive : public mutex {
public:
  mutex_recursive() : mutex(PTHREAD_MUTEX_RECURSIVE) { }
};

class mutex_lock {
  mutex& m_;
public:
  explicit mutex_lock(mutex& m) : m_(m) { m_.lock(); }
  ~mutex_lock() { m_.unlock(); }
};

class Semaphore {
  int _value, _wakeups;
  cond _cv;
public:
        explicit Semaphore(int value) :
  _value(value),
  _wakeups(0)
  {
    // nothing to do
  }

  ~Semaphore() {}

  inline void wait() {
    _cv.lock();
    _value--;
    if (_value < 0) {
      do {
        _cv.wait();
      } while(_wakeups < 1);
      _wakeups--;
    }
    _cv.unlock();
  }

  inline void signal() {
    _cv.lock();
    _value++;
    if(_value <= 0) {
      _wakeups++;
      _cv.signal();
    }
    _cv.unlock();
  }
};

#if defined(_POSIX_BARRIERS) && (_POSIX_BARRIERS - 20012L) >= 0
class barrier
{
  pthread_barrier_t _barrier;

public:
      explicit barrier(unsigned count) {

    pthread_barrier_init(&_barrier, NULL, count);
  }

  ~barrier() {
    pthread_barrier_destroy(&_barrier);
  }

  /// Return true if serial thread.
  inline bool wait() {
    return pthread_barrier_wait(&_barrier) == PTHREAD_BARRIER_SERIAL_THREAD;
  }
};

#else
// #  ifndef PTHREAD_BARRIER_SERIAL_THREAD
// #    define  PTHREAD_BARRIER_SERIAL_THREAD 1
// #  endif

class barrier
{
  int count; // required # of threads
  int current;    // current # of threads that have passed thru
  mutex barlock;  // protect current
  Semaphore barrier1; // implement the barrier
  Semaphore barrier2;

public:
      explicit barrier(unsigned cnt)
  : count(cnt), current(0), barrier1(0), barrier2(0) {
  }

  ~barrier() {}

  inline bool wait() {
    bool ret = false;
    barlock.lock();
    current += 1;
    if(current == count) {
      ret = true;
      for(int i=0; i<count;i++) {
        barrier1.signal();
      }
    }
    barlock.unlock();
    barrier1.wait(); // wait for n threads to arrive

    barlock.lock();
    current -= 1;
    if(current == 0) {
      for(int i=0;i<count;i++) {
        barrier2.signal();
      }
    }
    barlock.unlock();
    barrier2.wait();
    return ret;
  }
};

#endif
} //namespace pthread {

typedef pthread::cond cond;
typedef pthread::mutex mutex;
typedef pthread::barrier barrier;
} } // namespace jellyfish { namespace { locks

#endif /* __JELLYFISH_LOCKS_PTHREAD_HPP__ */
