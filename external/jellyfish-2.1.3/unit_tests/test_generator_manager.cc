#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <fstream>

#include <gtest/gtest.h>
#include <unit_tests/test_main.hpp>
#include <jellyfish/generator_manager.hpp>

namespace {
static const int nb_pipes = 5;
TEST(TmpPipes, CreateDestroy) {
  std::vector<std::string> pipe_paths;

  {
    jellyfish::tmp_pipes pipes(nb_pipes);
    EXPECT_EQ(nb_pipes, pipes.size());
    for(int i = 0; i < nb_pipes; ++i) {
      struct stat buf;
      ASSERT_EQ(0, stat(pipes[i], &buf));
      EXPECT_TRUE(S_ISFIFO(buf.st_mode));
      pipe_paths.push_back(pipes[i]);
    }
   }

  for(auto it = pipe_paths.begin(); it != pipe_paths.end(); ++it) {
    struct stat buf;
    EXPECT_EQ(-1, stat(it->c_str(), &buf));
  }
}

TEST(TmpPipes, Discard) {
  jellyfish::tmp_pipes pipes(nb_pipes);
  for(int i = 0; i < nb_pipes; ++i) {
    ASSERT_STRNE("", pipes[i]);
    std::string path(pipes[i]);
    int fd = open(path.c_str(), O_RDONLY|O_NONBLOCK);
    ASSERT_NE(-1, fd);
    pipes.discard(i);

    EXPECT_STREQ("", pipes[i]);
    struct stat buf;
    EXPECT_EQ(-1, stat(path.c_str(), &buf));
    char rbuf[1];
    EXPECT_EQ((ssize_t)0, read(fd, rbuf, 1));
    close(fd);
  }
}

TEST(GeneratorManager, OneLiners) {
  static const char* cmds_file = "./cmds_file";
  file_unlink unlink_cmds(cmds_file);

  {
    std::ofstream cmds(cmds_file, std::ios::out|std::ios::trunc);
    ASSERT_TRUE(cmds.good()) << "Failed to open cmd file '" << cmds_file << "'";
    cmds << "echo hello\n"
         << "date\n"
         << "uptime\n"
         << "uname\n";
    ASSERT_TRUE(cmds.good()) << "Failed to write to cmd file";
  }

  std::ifstream cmds(cmds_file);
  ASSERT_TRUE(cmds.good()) << "Failed top open cmd file '" << cmds_file << "'";
  static const int             nb_pipes     = 2;
  int                          active_pipes = nb_pipes;
  bool                         active[nb_pipes];
  jellyfish::generator_manager manager(cmds_file, nb_pipes);
  manager.start();
  for(int i = 0; i < nb_pipes; ++i)
    active[i] = true;
  std::vector<std::string> lines;
  while(active_pipes > 0) {
    for(int i = 0; i < nb_pipes; ++i) {
      if(!active[i])
        continue;
      std::ifstream p(manager.pipes()[i]);
      if(!p.good()) {
        active[i] = false;
        --active_pipes;
        continue;
      }
      std::string line;
      while(std::getline(p, line))
        lines.push_back(line);
    }
  }
  EXPECT_TRUE(manager.wait());

  EXPECT_EQ((size_t)4, lines.size());
}
}
