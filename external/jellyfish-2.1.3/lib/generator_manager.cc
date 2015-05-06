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

#include <config.h>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <signal.h>
#include <assert.h>

#include <cstdlib>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <sstream>

#include <jellyfish/generator_manager.hpp>
#include <jellyfish/err.hpp>

namespace jellyfish {
int open_cloexec(const char* path, int flags) {
#ifdef O_CLOEXEC
    int fd = open(path, flags|O_CLOEXEC);
#else
    int fd = open(path, flags);
    if(fd != -1)
      fcntl(fd, F_SETFD, FD_CLOEXEC);
#endif
    return fd;
}

std::string tmp_pipes::create_tmp_dir() {
  std::vector<const char*> prefixes;
  const char* tmpdir = getenv("TMPDIR");
  if(tmpdir)
    prefixes.push_back(tmpdir);
#ifdef P_tmpdir
  prefixes.push_back(P_tmpdir);
#endif
  prefixes.push_back(".");

  for(auto it = prefixes.begin(); it != prefixes.end(); ++it) {
    size_t len = strlen(*it) + 6 + 1;
    std::unique_ptr<char[]> tmppath(new char[len]);
    sprintf(tmppath.get(), "%sXXXXXX", *it);
    const char* res = mkdtemp(tmppath.get());
    if(res)
      return std::string(res);
  }
  eraise(std::runtime_error) << "Failed to create a temporary directory for the pipes. Set the variable TMPDIR properly";
  return "";
}

std::vector<std::string> tmp_pipes::create_pipes(const std::string& tmpdir, int nb_pipes)
{
  std::vector<std::string> pipes;
  for(int i = 0; i < nb_pipes; ++i) {
    std::ostringstream path;
    path << tmpdir << "/fifo" << i;
    if(mkfifo(path.str().c_str(), S_IRUSR|S_IWUSR) == -1)
      eraise(std::runtime_error) << "Failed to create named fifos";
    pipes.push_back(path.str());
  }
  return pipes;
}

void tmp_pipes::discard(int i) {
  if(pipes_[i].empty())
    return;
  // First rename the fifo so no new reader will open it, then open
  // the fifo (with its new name for writing, in non-blocking mode. If
  // we get a valid file descriptor, some readers are blocked reading
  // on the fifo: we close the fifo and free the readers. Otherwise,
  // no readers is blocked and no action is required. Finally we
  // unlink the fifo for good.
  std::string discarded_name(pipes_[i]);
  discarded_name += "_discarded";
  if(rename(pipes_[i].c_str(), discarded_name.c_str()) == -1)
    return;
  pipes_[i].clear();
  pipes_paths_[i] = 0;
  int fd = open(discarded_name.c_str(), O_WRONLY|O_NONBLOCK);
  if(fd != -1)
    close(fd);
  unlink(discarded_name.c_str());
}

void tmp_pipes::cleanup() {
  for(size_t i = 0; i < pipes_.size(); ++i) {
    discard(i);
  }
  rmdir(tmpdir_.c_str());
}

void generator_manager::start()  {
  if(manager_pid_ != -1)
    return;
  manager_pid_ = fork();
  switch(manager_pid_) {
  case -1:
    eraise(std::runtime_error) << "Failed to start manager process";
    break;
  case 0:
    manager_pid_ = -1;
    break;
  default:
    cmds_.close();
    return;
  }


  // In child
  setup_signal_handlers();
  start_commands(); // child start commands
  int signal = kill_signal_;
  if(signal == 0)
    exit(EXIT_SUCCESS);

  // Got killed. Kill all children, cleanup and kill myself with the
  // same signal (and die from it this time :). We do not wait on the
  // dead children as we are going to die soon as well, and we don't
  // care about the return value at that point. Let init take care of
  // that for us...
  cleanup();
  unset_signal_handlers();
  kill(getpid(), signal); // kill myself
  exit(EXIT_FAILURE); // Should not be reached
}

static generator_manager* manager = 0;
void generator_manager::signal_handler(int signal) {
  manager->kill_signal_ = signal;
}
void generator_manager::setup_signal_handlers() {
  int res;
  struct sigaction act;
  memset(&act, '\0', sizeof(act));
  act.sa_handler = signal_handler;
  res = sigaction(SIGTERM, &act, 0);
  assert(res == 0);
  // Should we redefine other signals as well? Like SIGINT, SIGQUIT?
}

void generator_manager::unset_signal_handlers() {
  struct sigaction act;
  memset(&act, '\0', sizeof(act));
  act.sa_handler = SIG_DFL;
  assert(sigaction(SIGTERM, &act, 0) == 0);
}

bool generator_manager::wait() {
  if(manager_pid_ == -1) return false;
  pid_t pid = manager_pid_;
  manager_pid_ = -1;
  int status;
  if(pid != waitpid(pid, &status, 0))
    return false;
  return WIFEXITED(status) && (WEXITSTATUS(status) == 0);
}

void generator_manager::cleanup() {
  for(auto it = pid2pipe_.begin(); it != pid2pipe_.end(); ++it) {
    kill(it->first, SIGTERM);
    pipes_.discard(it->second.pipe);
  }
  pipes_.cleanup();
}

void generator_manager::start_one_command(const std::string& command, int pipe)
{
  cmd_info_type info = { command, pipe };
  pid_t child = fork();
  switch(child) {
  case -1:
    std::cerr << "Failed to fork. Command '" << command << "' not run" << std::endl;
    return;
  case 0:
    break;
  default:
    pid2pipe_[child] = info;
    return;
  }

  // In child
  int dev_null = open_cloexec("/dev/null", O_RDONLY);
  if(dev_null != -1)
    dup2(dev_null, 0);

  int pipe_fd = open_cloexec(pipes_[pipe], O_WRONLY);
  if(pipe_fd == -1) {
    std::cerr << "Failed to open output pipe. Command '" << command << "' not run" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(dup2(pipe_fd, 1) == -1) {
    std::cerr << "Failed to dup pipe to stdout. Command '" << command << "' not run" << std::endl;
    exit(EXIT_FAILURE);
  }
  execl(shell_, shell_, "-c", command.c_str(), (char*)0);
  std::cerr << "Failed to exec. Command '" << command << "' not run" << std::endl;
  exit(EXIT_FAILURE);
}

std::string generator_manager::get_cmd() {
  std::string command;

  while(std::getline(cmds_, command)) {
    size_t pos = command.find_first_not_of(" \t\n\v\f\r");
    if(pos != std::string::npos && command[pos] != '#')
      break;
    command.clear();
  }
  return command;
}

void generator_manager::start_commands()
{
  std::string command;
  size_t i;
  for(i = 0; i < pipes_.size(); ++i) {
    command = get_cmd();
    if(command.empty())
      break;
    start_one_command(command, i);
  }
  for( ; i < pipes_.size(); ++i)
    pipes_.discard(i);

  while(!pid2pipe_.empty()) {
    int status;
    int res = ::wait(&status);
    if(res == -1) {
      if(errno == EINTR) continue;
      break;
    }
    cmd_info_type info = pid2pipe_[res];
    pid2pipe_.erase(info.pipe);
    command = get_cmd();
    if(!command.empty()) {
      start_one_command(command, info.pipe);
    } else {
      pipes_.discard(info.pipe);
    }
    if(!display_status(status, info.command)) {
      cleanup();
      exit(EXIT_FAILURE);
    }
  }
}

bool generator_manager::display_status(int status, const std::string& command)
{
  if(WIFEXITED(status) && WEXITSTATUS(status) != 0) {
    std::cerr << "Command '" << command
              << "' exited with error status " << WEXITSTATUS(status) << std::endl;
    return false;
  } else if(WIFSIGNALED(status)) {
    std::cerr << "Command '" << command
              << "' killed by signal " << WTERMSIG(status) << std::endl;
    return false;
  }
  return true;
}

} // namespace jellyfish
