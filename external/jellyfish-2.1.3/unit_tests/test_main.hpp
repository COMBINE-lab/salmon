#ifndef __TEST_MAIN_HPP__
#define __TEST_MAIN_HPP__

#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>

#include <string>

uint64_t random_bits(int length);
inline uint64_t random_bits() { return random_bits(64); }

struct file_unlink {
  std::string path;
  bool do_unlink;
  explicit file_unlink(const char* s, bool d = true) : path(s), do_unlink(d) { }
  explicit file_unlink(const std::string& s, bool d = true) : path(s), do_unlink(d) { }

  ~file_unlink() {
    if(do_unlink)
      unlink(path.c_str());
  }
};

#endif /* __TEST_MAIN_HPP__ */

