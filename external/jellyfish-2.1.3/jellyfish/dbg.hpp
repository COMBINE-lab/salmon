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

#ifndef __DBG_HPP__
#define __DBG_HPP__

#include <iostream>
#include <iomanip>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <signal.h>

#include <jellyfish/time.hpp>

namespace dbg {
  pid_t gettid();

  class stringbuf : public std::stringbuf {
  public:
    stringbuf() : std::stringbuf(std::ios_base::out) { }
    explicit stringbuf(const std::string &str) : 
      std::stringbuf(str, std::ios_base::out) { }

    bool end_is_space() {
      if(pptr() == pbase())
        return true;
      return isspace(*(pptr() - 1));
    }
    friend class print_t;
  };

  class str {
    const char  *_s;
    const size_t _l;
  public:
    str(const char *s, size_t len) : _s(s), _l(len) {}
    friend class print_t;
  };

  class xspace { };
  class no_flush { };

  class print_t {
    static pthread_mutex_t _lock;
    static volatile pid_t  _print_tid;

    stringbuf              _strbuf;
    std::ostream           _buf;
    bool                   _flush;
  public:
    print_t(const char *file, const char *function, int line) :
      _buf(&_strbuf), _flush(true)
    {
      const char *file_basename = strrchr(file, '/');
      if(!file_basename)
        file_basename = file;
      _buf << pthread_self() << "/" << gettid() << ":"
           << file_basename << ":" << function << ":" << line << ": ";
    }

    ~print_t() {
      if(_print_tid == 0 || gettid() == _print_tid) {
        pthread_mutex_lock(&_lock);
        std::cerr.write(_strbuf.pbase(), _strbuf.pptr() - _strbuf.pbase());
        if(_flush)
          std::cerr << std::endl;
        else
          std::cerr << "\n";
        pthread_mutex_unlock(&_lock);
      }
    }

    static int set_signal(int signum = SIGUSR1);
    static void signal_handler(int signum, siginfo_t *info, void *context);
    static pid_t print_tid() { return _print_tid; }
    static void print_tid(pid_t new_tid) { _print_tid = new_tid; }

    print_t & operator<<(const char *a[]) {
      for(int i = 0; a[i]; i++)
        _buf << (i ? "\n" : "") << a[i];
      return *this;
    }
    print_t & operator<<(const std::exception &e) {
      _buf << e.what();
      return *this;
    }
    print_t & operator<<(const str &ss) {
      _buf.write(ss._s, ss._l);
      return *this;
    }
    print_t & operator<<(const xspace &xs) {
      if(!_strbuf.end_is_space())
        _buf << " ";
      return *this;
    }
    print_t &operator<<(const no_flush &nf) {
      _flush = false;
      return *this;
    }
    print_t & operator<<(const Time &t) {
      _buf << t.str();
      return *this;
    }
    template<typename T>
    print_t & operator<<(const T &x) {
      _buf << x;
      return *this;
    }
  };

  class no_print_t {
  public:
    no_print_t() {}
    
    template<typename T>
    no_print_t & operator<<(const T &x) { return *this; }
  };

  void tic();
  Time toc();
}

#ifdef DEBUG
#define DBG if(1) dbg::print_t(__FILE__, __FUNCTION__, __LINE__)
#define NFDBG if(1) dbg::print_t(__FILE__, __FUNCTION__, __LINE__) << dbg::no_flush()
#define V(v) dbg::xspace() << #v ":" << v
#else
#define DBG if(1) dbg::no_print_t()
#define NFDBG if(1) dbg::no_print_t()
#define V(v) v
#endif

#endif /* __DBG_HPP__ */
