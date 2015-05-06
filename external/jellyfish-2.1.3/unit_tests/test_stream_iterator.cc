#include <stdlib.h>
#include <unistd.h>

#include <gtest/gtest.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <jellyfish/err.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/stream_iterator.hpp>

namespace {
typedef std::vector<const char*>                  path_vector;
typedef path_vector::const_iterator               path_iterator;
typedef jellyfish::stream_iterator<path_iterator> stream_iterator;

class StreamIterator : public ::testing::Test {
protected:
  static path_vector paths;
  static int         total_lines;
  static char*       tmpdir;

  static void SetUpTestCase() {
    tmpdir = strdup("/tmp/stream_iterator_XXXXXX");
    if(!mkdtemp(tmpdir))
      eraise(std::runtime_error) << "Failed to create tmp directory" << jellyfish::err::no;

    int nb_files = jellyfish::random_bits(5);
    for(int i = 0; i < nb_files; ++i) {
      std::ostringstream path;
      path << tmpdir << "/" << i;
      paths.push_back(strdup(path.str().c_str()));
      std::ofstream tmp_file(path.str().c_str());
      if(tmp_file.fail())
        eraise(std::runtime_error) << "Failed to open file '" << path.str() << "' for writing";
      int nb_lines = jellyfish::random_bits(6);
      total_lines += nb_lines;
      for(int j = 0; j < nb_lines; ++j)
        tmp_file << "line " << j << "\n";
    }
  }

  static void TearDownTestCase() {
    for(path_iterator it = paths.begin(); it != paths.end(); ++it) {
      if(unlink(*it) == -1)
        eraise(std::runtime_error) << "Failed to unlink file '" << *it << jellyfish::err::no;
      free((void*)*it);
    }
    if(rmdir(tmpdir) == -1)
      eraise(std::runtime_error) << "Failed to rmdir '" << tmpdir << jellyfish::err::no;
    free(tmpdir);
  }
};
path_vector StreamIterator::paths;
int         StreamIterator::total_lines = 0;
char*       StreamIterator::tmpdir;

TEST_F(StreamIterator, Empty) {
  stream_iterator sit(paths.begin(), paths.begin());
  stream_iterator send;

  int nb_files = 0;
  int nb_lines = 0;
  for( ; sit != send; ++sit, ++nb_files)
    for(std::string line; std::getline(*sit, line); ++nb_lines) ;

  EXPECT_EQ(0, nb_files);
  EXPECT_EQ(0, nb_lines);
}


TEST_F(StreamIterator, RandomFiles) {
  SCOPED_TRACE(::testing::Message() << "nb_files:" << paths.size() << " nb_lines:" << total_lines);

  stream_iterator sit(paths.begin(), paths.end());
  stream_iterator send;

  int nb_files = 0;
  int nb_lines = 0;
  for( ; sit != send; ++sit, ++nb_files) {
    for(std::string line; std::getline(*sit, line); ++nb_lines) ;
    EXPECT_TRUE(sit->eof());
  }

  EXPECT_EQ(paths.size(), (size_t)nb_files);
  EXPECT_EQ(total_lines, nb_lines);
}
}
