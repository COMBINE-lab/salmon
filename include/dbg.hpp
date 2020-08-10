#ifndef __DEBOUTPUT_H__
#define __DEBOUTPUT_H__

#include <iostream>
#include <sstream>

struct dbg {
  std::ostream&      os_;
  std::ostringstream ss_;
  const bool         flush_;

  dbg(std::ostream& os = std::cerr, bool flush = false)
    : os_(os), flush_(flush)
  { }
  ~dbg() {
    os_ << ss_.str();
    if(flush_)
      os_ << std::flush;
  }
};

template<typename T>
dbg& operator<<(dbg&& d, const T& x) {
  return d << x;
}

// template<typename T>
// dbg& operator<<(dbg d, const T& x) {
//   d.ss_ << x;
//   return d;
// }

template<typename T>
dbg& operator<<(dbg& d, const T& x) {
  d.ss_ << x;
  return d;
}

#endif /* __DEBOUTPUT_H__ */
