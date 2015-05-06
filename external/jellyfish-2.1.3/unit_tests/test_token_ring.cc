#include <iostream>
#include <sstream>
#include <gtest/gtest.h>
#include <jellyfish/token_ring.hpp>
#include <jellyfish/thread_exec.hpp>

namespace {
using jellyfish::thread_exec;
struct write_number : public thread_exec {
  typedef jellyfish::token_ring<jellyfish::locks::pthread::cond> token_ring;

  std::ostream *os_;
  int           nb_threads_;
  token_ring    ring;

public:
  write_number(std::ostream* os, int nb_threads) :
    os_(os), nb_threads_(nb_threads), ring(nb_threads)
  { }

  virtual void start(int i) {
    token_ring::token& token = ring[i];
    for(int j = 0; j < 10; ++j) {
      token.wait();
      *os_ << (nb_threads_ * j + i) << "\n";
      token.pass();
    }
  }
};

TEST(TokenRing, Count) {
  static const int nb_threads = 5;
  std::ostringstream out;

  write_number write(&out, nb_threads);
  write.exec_join(nb_threads);

  std::ostringstream res;
  for(int i = 0; i < nb_threads * 10; ++i)
    res << i << "\n";

  EXPECT_EQ(res.str(), out.str());
}
} // namespace {
