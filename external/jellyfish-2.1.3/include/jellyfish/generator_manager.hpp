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


#ifndef __JELLYFISH_SPAWN_EXTERNAL_HPP_
#define __JELLYFISH_SPAWN_EXTERNAL_HPP_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <iostream>
#include <vector>
#include <map>
#include <stdexcept>

#ifdef HAVE_EXT_STDIO_FILEBUF_H
#include <ext/stdio_filebuf.h>
#define STDIO_FILEBUF_TYPE __gnu_cxx::stdio_filebuf<std::istream::char_type>
#else
#include <jellyfish/stdio_filebuf.hpp>
#define STDIO_FILEBUF_TYPE jellyfish::stdio_filebuf<std::istream::char_type>
#endif

#include <jellyfish/err.hpp>

namespace jellyfish {
// Open a path and set CLOEXEC flags
int open_cloexec(const char* path, int flags);

// Input stream (inherit from std::istream, behaves mostly like an
// ifstream), with flag O_CLOEXEC (close-on-exec) turned on.
class cloexec_istream : public std::istream
{
  static std::streambuf* open_file(const char* path) {
    int fd = open_cloexec(path, O_RDONLY);
    return new STDIO_FILEBUF_TYPE(fd, std::ios::in);
  }

public:
  cloexec_istream(const cloexec_istream&) = delete;
  cloexec_istream(const char* path) :
  std::istream((open_file(path)))
  { }
  cloexec_istream(const std::string& path) :
    std::istream(open_file(path.c_str()))
  { }
  virtual ~cloexec_istream() { close(); }
  void close() { delete std::istream::rdbuf(0); }
};


// This class is responsible for creating a tmp directory and
// populating it with fifos.
class tmp_pipes {
  static std::string create_tmp_dir();
  std::vector<std::string> create_pipes(const std::string& tmpdir, int nb_pipes);

  std::string              tmpdir_;
  std::vector<std::string> pipes_;
  std::vector<const char*> pipes_paths_;

public:
  tmp_pipes(int nb_pipes):
  tmpdir_(create_tmp_dir()),
  pipes_(create_pipes(tmpdir_, nb_pipes))
  {
    for(auto it = pipes_.cbegin(); it != pipes_.cend(); ++it)
      pipes_paths_.push_back(it->c_str());
  }
  ~tmp_pipes() { cleanup(); }

  size_t size() const { return pipes_.size(); }
  const char* operator[](int i) const { return pipes_[i].c_str(); }
  std::vector<const char*>::const_iterator begin() const { return pipes_paths_.cbegin(); }
  std::vector<const char*>::const_iterator end() const { return pipes_paths_.cend(); }

  // Discard a pipe: unlink it while it is open for writing. The
  // reading process will get no data and won't be able to reopen the
  // file, marking the end of this pipe.
  void discard(int i);
  // Discard all pipes
  void cleanup();
};

// This class creates a new process which manages a bunch of
// "generators", sub-processes that writes into a fifo (named pipe)
// and generate sequence.
class generator_manager {
  cloexec_istream cmds_;
  tmp_pipes       pipes_;
  pid_t           manager_pid_;
  const char*     shell_;
  int             kill_signal_; // if >0, process has received that signal

  struct cmd_info_type {
    std::string command;
    int         pipe;
  };
  typedef std::map<pid_t, cmd_info_type> pid2pipe_type;
  pid2pipe_type pid2pipe_;

public:
  generator_manager(const char* cmds, int nb_pipes, const char* shell = 0) :
    cmds_(cmds),
    pipes_(nb_pipes),
    manager_pid_(-1),
    shell_(shell),
    kill_signal_(0)
  {
    if(!cmds_.good())
      eraise(std::runtime_error) << "Failed to open cmds file '" << cmds << "'";
    if(!shell_)
      shell_ = getenv("SHELL");
    if(!shell_)
      shell_ = "/bin/sh";
  }
  ~generator_manager() { wait(); }

  const tmp_pipes& pipes() const { return pipes_; }
  pid_t pid() const { return manager_pid_; }

  // Start the manager process
  void start();
  // Wait for manager process to finish. Return true if it finishes
  // with no error, false otherwise.
  bool wait();

private:
  /// Read commands from the cmds stream. There is one command per
  /// line. Empty lines or lines whose first non-white space character
  /// is a # are ignored. Return an empty string when no more commands
  /// are available.
  std::string get_cmd();
  void start_commands();
  void start_one_command(const std::string& command, int pipe);
  bool display_status(int status, const std::string& command);
  void setup_signal_handlers();
  void unset_signal_handlers();
  static void signal_handler(int signal);
  void cleanup();
};
}

#endif /* __JELLYFISH_SPAWN_EXTERNAL_HPP_ */
