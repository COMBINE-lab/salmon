#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstdio>
#include <fstream>
#include <istream>
#include <memory>

#include <gtest/gtest.h>

#include <unit_tests/test_main.hpp>
#include <jellyfish/stdio_filebuf.hpp>

namespace {
typedef jellyfish::stdio_filebuf<char> stdio_filebuf;
TEST(StdioFileBuf, Read) {
  const char* file_name = "test_stdio_file_buf_read";
  file_unlink fu(file_name);
  std::vector<char> buffer(8192);
  for(size_t i = 0; i < buffer.size(); ++i)
    buffer[i] = '@' + (i % 64);
  {
    std::ofstream out(file_name);
    out.write(buffer.data(), buffer.size());
  }

  int fd = open(file_name, O_RDONLY);
  EXPECT_GT(fd, -1);
  FILE* file = fopen(file_name, "r");
  EXPECT_NE((FILE*)0, file);

  std::unique_ptr<stdio_filebuf> fd_buf(new stdio_filebuf(fd, std::ios::in));
  std::istream fd_stream(fd_buf.get());
  std::unique_ptr<stdio_filebuf> file_buf(new stdio_filebuf(file, std::ios::in));
  std::istream file_stream(file_buf.get());

  size_t have_read = 0;
  char read_buffer[256];
  while(have_read < buffer.size()) {
    const size_t to_read     = random_bits(8);
    const size_t expect_read = std::min(buffer.size() - have_read, to_read);

    fd_stream.read(read_buffer, to_read);
    EXPECT_EQ(expect_read, fd_stream.gcount());
    EXPECT_EQ(0, std::strncmp(&buffer[have_read], read_buffer, expect_read));

    file_stream.read(read_buffer, to_read);
    EXPECT_EQ(expect_read, fd_stream.gcount());
    EXPECT_EQ(0, std::strncmp(&buffer[have_read], read_buffer, expect_read));

    have_read += expect_read;
  }
  EXPECT_TRUE(fd_stream.eof());
}
}
