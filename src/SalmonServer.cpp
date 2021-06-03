#include "SalmonServer.hpp"
#include <asm-generic/socket.h>
#include <bits/types/struct_iovec.h>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <cerrno>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <map>
#include <poll.h>
#include <signal.h>
#include <stack>
#include <stdexcept>
#include <sys/socket.h>
#include <sys/types.h> /* See NOTES */
#include <sys/types.h>
#include <sys/un.h>
#include <sys/wait.h>
#include <unistd.h>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <utility>
namespace po = boost::program_options;

// Same error status as salmoneServer
static int serveIndex(const std::string& socketPath, const std::string& indexPath, int notifyFD,
                int& argc, argvtype& argv,
                std::unique_ptr<SalmonIndex>& salmonIndex);
static int serverMainLoop(int&argc, argvtype& argv, int unix_socket);
static int contactServer(const std::string& socketPath, const std::string& subcommand, const std::vector<std::string>& opts);

int salmonServerClient(const std::string& subcommand, int argc, argvtype argv) {
  using std::string;
  try {
    po::options_description cliento("client options");
    cliento.add_options()
      ("help,h", "this message")
      ("server", po::value<string>(), "server socket path")
      ("args", po::value<std::vector<std::string>>(), "Args");
    po::positional_options_description pd;
    pd.add("args", -1);
    po::command_line_parser parser{argc, argv};
    parser.options(cliento).positional(pd).allow_unregistered();
    auto parsed = parser.run();
    po::variables_map vm;
    po::store(parsed, vm);
    if(vm.count("help") && vm.count("server")) {
      const auto hstring = R"(
Client
==========
Contact the server listening on the UNIX socket given by --server to do the
processing. The subcommand and every switch is sent to the server. Rerun without
the --server switch to get the help message of the actual subcommand.
)";
      std::cout << hstring << '\n'
                << cliento << std::endl;
      return 0;
    }
    po::notify(vm);

    if (!vm.count("server")) // No --server, do nothing
      return -1;

    auto opts = po::collect_unrecognized(parsed.options, po::include_positional);
    return contactServer(vm["server"].as<string>(), subcommand, opts);
  } catch (std::exception& e) {
    std::cerr << "Client error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Client unknown error" << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_FAILURE; // Should never be reached
}

int salmonServerServer(int& argc, argvtype& argv,
                       std::unique_ptr<SalmonIndex>& salmonIndex) {
  using std::string;
  try {
    po::options_description servero("Server options");
    servero.add_options()
      ("help,h", "print this message")
      ("index,i", po::value<string>()->required(), "path to salmon index")
      ("server-notify", po::value<int>()->default_value(-1), "file descriptor to send ready notification")
      ("socket", po::value<std::string>()->required(), "server socket path");
    po::positional_options_description pd;
    pd.add("socket", 1);
    po::command_line_parser parser{argc, argv};
    parser.positional(pd).options(servero);
    auto parsed = parser.run();
    po::variables_map vm;
    po::store(parsed, vm);

    if(vm.count("help")) {
      const auto hstring = R"(
Server
==========
Start the salmon server holding the index given by --index. Wait for requests
from clients from the UNIX socket at path given by --socket.
)";
      std::cout << hstring << '\n'
                << servero << std::endl;
      return 0;
    }
    po::notify(vm);

    return serveIndex(vm["socket"].as<string>(), vm["index"].as<string>(),
                      vm["server-notify"].as<int>(), argc, argv, salmonIndex);
  } catch (std::exception& e) {
    std::cerr << "Server error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Server unknown error" << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_FAILURE; // Should never be reached
}

int salmonServer(int /* argc */, argvtype /* argv */, std::unique_ptr<SalmonIndex>& /* salmonIndex */) {
  throw std::runtime_error("BUG: this function should never be called!");
}

void cpperror(const std::string& msg) {
  throw std::runtime_error(msg + ": " + strerror(errno));
}

// utility functions
//

template <typename C> size_t dataSize(const C& container) {
  return container.size() * sizeof(typename C::value_type);
}

// Send full buffer to fd, restarting after interruptions or partial send
ssize_t sendBuffer(int fd, const char* buf, size_t size) {
  size_t offset = 0;
  while (offset < size) {
    ssize_t sent = send(fd, buf + offset, size - offset, 0);
    if (sent == -1) {
      if (errno == EINTR)
        continue;
      return -1;
    }
    offset += sent;
  }
  return offset;
}
// Same. Convenience function for std::vector<type> for any type
template <typename C> ssize_t sendBuffer(int fd, const C& container) {
  return sendBuffer(fd, (const char*)container.data(), dataSize(container));
}

// Receive exactly size bytes from fd, restarting after interruptions and
// partial receive
ssize_t receiveBuffer(int fd, char* buffer, size_t size) {
  size_t offset = 0;
  while (offset < size) {
    ssize_t received = recv(fd, buffer + offset, size - offset, 0);
    if (received == -1) {
      if (errno == EINTR)
        continue;
      return -1;
    }
    if (received == 0)
      return 0;
    offset += received;
  }
  return offset;
}
template <typename C> ssize_t receiveBuffer(int fd, C& container) {
  return receiveBuffer(fd, (char*)container.data(), dataSize(container));
}

// Defer utility objects
struct deferClose {
  int fd_;
  deferClose(int fd) : fd_(fd) {}
  ~deferClose() { close(fd_); }
};
struct deferCloseDir {
  DIR* d;
  deferCloseDir(DIR* dir) : d(dir) {}
  ~deferCloseDir() { closedir(d); }
};
// Unlink file, but only if we are the parent process
struct deferUnlink {
  static bool parent;
  const char* path_;
  deferUnlink(const char* path) : path_(path) {}
  ~deferUnlink() {
    if (parent)
      unlink(path_);
  }
};
bool deferUnlink::parent = true;

int serveIndex(const std::string& socketPath, const std::string& indexPath,
               int notifyFd, int& argc, argvtype& argv,
               std::unique_ptr<SalmonIndex>& salmonIndex) {
  // Create unix socket and bind it
  struct sockaddr_un addr;
  if (socketPath.size() >= sizeof(addr.sun_path) - 1)
    cpperror("Path for server socket is too long");
  int unix_socket = socket(AF_UNIX, SOCK_STREAM, 0);
  if (unix_socket == -1)
    cpperror("Failed to open the server socket");
  deferClose closeSocket(unix_socket);

  addr.sun_family = AF_UNIX;
  strcpy(addr.sun_path, socketPath.c_str());
  if (bind(unix_socket, (struct sockaddr*)&addr, sizeof(addr)) == -1)
    cpperror("Failed to bind server socket");
  deferUnlink unlinkSocket(addr.sun_path);

  // Load index
  boost::filesystem::path indexDirectory(indexPath);
  auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
  //  consoleSink->set_color(spdlog::level::warn, consoleSink->magenta);
  auto consoleLog = spdlog::create("servertderrLog", {consoleSink});
  salmonIndex = checkLoadIndex(indexDirectory, consoleLog);

  // Ready to serve
  if (notifyFd != -1) {
    sendBuffer(notifyFd, "ready", 6);
    close(notifyFd);
  }

  // This in a child process with -1 to continue execution of Salmon, or in the
  // parent process with an error code
  return serverMainLoop(argc, argv, unix_socket);
}

static volatile bool done = false;
void term_handler(int) { done = true; }
void chld_handler(int) { /* Do nothing */
}

struct DefineSignals {
  struct sigaction actionInt, actionTerm, actionChld, actionPipe;

  void setup() {
    struct sigaction act;
    memset(&act, '\0', sizeof(act));
    act.sa_handler = term_handler;
    if (sigaction(SIGINT, &act, &actionInt) == -1 ||
        sigaction(SIGTERM, &act, &actionTerm) == -1)
      cpperror("Error redirecting termination signals");
    act.sa_handler = chld_handler;
    if (sigaction(SIGCHLD, &act, &actionChld) == -1)
      cpperror("Error ignoring SIGCHLD");
    if (sigaction(SIGPIPE, &act, &actionPipe) == -1)
      cpperror("Error ignore SIGPIPE");
  }

  void reset() {
    sigaction(SIGINT, &actionInt, nullptr);
    sigaction(SIGTERM, &actionTerm, nullptr);
    sigaction(SIGCHLD, &actionChld, nullptr);
    sigaction(SIGPIPE, &actionPipe, nullptr);
  }
};

//
// Server side, parent process
//
// The server accepts connections on the UNIX socket. For each new connection,
// it forks a new child process that will do the actual work. Keep track of
// the connection "child pid" -> "socket fd". When a child dies, write the
// exit status on the socket fd (for the client to report), then close the
// socket fd.
//
// If a socket fd hangs up (poll returns HUP), which means the client closed
// its side (e.g., the client got Ctrl-C and died), then we send SIGTERM to
// the corresponding child process that is doing the work on behalf of this
// client.
//

// When a child is done (process waited for), send it it's status and close
// socket
void handleDoneChildren(std::map<pid_t, int> & childrenSocket) {
  while (true) {
    int status;
    pid_t pid = waitpid(-1, &status, WNOHANG);
    if (pid == 0)
      break; // No child has stopped
    if (pid == -1) {
      if (errno == EINTR)
        continue;
      if (errno == ECHILD)
        break;
      std::cerr << "Warning: error while waiting for a child: "
                << strerror(errno) << std::endl;
      break;
    }
    auto it = childrenSocket.find(pid);
    if (it == childrenSocket.end()) {
      std::cerr << "Warning: caught unknown child process " << pid
                << std::endl;
      continue;
    }

    if (it->second !=
        -1) { // If -1, already closed by client, so can't send status
      while (true) {
        auto sent = send(it->second, &status, sizeof(status), 0);
        if (sent == -1) {
          if (errno == EINTR)
            continue;
          std::cerr << "Warning: failed to send status to process " << pid
                    << ' ' << strerror(errno) << std::endl;
        }
        break;
      }
      close(it->second);
    }
    childrenSocket.erase(it);
  }
}

int handleChild(int fd, int& argc, argvtype& argv);
int serverMainLoop(int& argc, argvtype& argv, int unix_socket) {
  DefineSignals signals;
  signals.setup();

  std::cerr << "Server waiting for requests. Ctrl-C to stop." << std::endl;
  if (listen(unix_socket, 5) == -1)
    cpperror("Error listening on unix socket");

  std::map<pid_t, int> childrenSocket;
  std::vector<struct pollfd> pollfds;
  std::vector<pid_t> pollpids;

  while (!done) {
    handleDoneChildren(childrenSocket);

    pollfds.resize(childrenSocket.size() + 1);
    pollpids.resize(childrenSocket.size() + 1);
    pollfds[0].fd = unix_socket;
    pollfds[0].events = POLLIN;
    pollfds[0].revents = 0;
    int i = 1;
    for (const auto& child : childrenSocket) {
      pollfds[i].fd = child.second;
      pollfds[i].events = 0;
      pollfds[i].revents = 0;
      pollpids[i] = child.first;
      ++i;
    }

    // Can't use ppoll (not POSIX, not on mac) to avoid race conditions, use a
    // timeout
    int res = poll(pollfds.data(), pollfds.size(), 10000);
    if (res == 0)
      continue;
    if (res == -1) {
      if (errno == EINTR)
        continue;
      cpperror("Error polling file descriptors");
    }

    if (pollfds[0].revents & POLLIN) {
      int fd = accept(pollfds[0].fd, nullptr, nullptr);
      if (fd == -1) {
        if (errno == EINTR)
          continue;
        cpperror("Error accepting on unix socket");
      }

      pid_t pid = fork();
      switch (pid) {
      case -1:
        std::cerr << "Warning: failed to create child process: "
                  << strerror(errno) << std::endl;
        close(fd); // Summary termination error sent to client
        break;

      case 0:
        deferUnlink::parent = false;
        signals.reset();
        for (const auto poll : pollfds)
          close(poll.fd);
        return handleChild(fd, argc, argv);
        break;

      default:
        childrenSocket[pid] = fd;
        break;
      }
    }

    for (size_t i = 1; i < pollfds.size(); ++i) {
      if (pollfds[i].revents &
          POLLHUP) { // Client closed it's socket. Forget about it
        close(pollfds[i].fd);
        kill(pollpids[i], SIGTERM);
        childrenSocket[pollpids[i]] = -1;
      }
    }
  }

  close(unix_socket);
  std::cerr << "Waiting for remaining children" << std::endl;
  handleDoneChildren(childrenSocket);
  return 0;
}

//
// Server side, Child process.
//
// Receive arguments and file descriptor.
// After resetting the environment, return with -1 to keep on processing
// with the quantification
///

// dup file descriptor from[i] to to[i] and close from[i] so as to reproduce
// the fd environment of the client. The difficulty is that dup may enter in
// collision (close) with another from[j]. First we detect chain of collisions
// to do the dups in the right order. In case of a collision cycle, break the
// cycle by duping a fd from the cycle to a very high number first.
void dupClose(int from, int to) {
  if (dup2(from, to) == -1)
    cpperror("Dup file descriptors");
  close(from);
}
void redirectFds(const std::vector<int>& from, const std::vector<int>& to) {
  std::map<int, int> fromTo;
  int max = 0;
  for (size_t i = 0; i < from.size(); ++i) {
    max = std::max(max, from[i]);
    max = std::max(max, to[i]);
    fromTo[from[i]] = to[i];
  }

  std::stack<std::pair<int, int>> stack;
  while (!fromTo.empty()) {
    // Find the dependency chain
    stack.push(*fromTo.begin());
    int start = stack.top().first;
    bool cycle = false;
    while (true) {
      auto it = fromTo.find(stack.top().second);
      if (it == fromTo.end())
        break;                   // End of dependency chain
      if (it->second == start) { // Found a cycle. Break it
        dupClose(it->first, max + 1);
        cycle = true;
        fromTo.erase(it);
        break;
      }
      stack.push(*it);
    }

    // Unwind stack / dependency chain
    while (!stack.empty()) {
      dupClose(stack.top().first, stack.top().second);
      fromTo.erase(stack.top().first);
      stack.pop();
    }

    // Last one, the cycle breaker if any. It was moved to max+1 and is
    // supposed to go to start.first.
    if (cycle)
      dupClose(max + 1, start);
  }
}

constexpr size_t maxFds = 100; // Maximum number of file descriptors to send
static constexpr size_t maxArgvLen = 1024 * 1024; // Maximum argv size
static std::vector<char> rawArgv;
static std::vector<const char*> childArgv;
int handleChild(int socket, int& argc, argvtype& argv) {
  size_t offset = 0;

  size_t lens[2];
  struct iovec io = {.iov_base = lens, .iov_len = sizeof(lens)};
  std::vector<int> fds;
  std::vector<char> msgbuf(
      CMSG_SPACE(maxFds * sizeof(decltype(fds)::value_type)));
  struct msghdr msg;
  memset(&msg, '\0', sizeof(msg));
  msg.msg_iov = &io;
  msg.msg_iovlen = 1;
  msg.msg_control = msgbuf.data();
  msg.msg_controllen = msgbuf.size();
  while (true) {
    ssize_t received = recvmsg(socket, &msg, 0);
    if (received == -1) {
      if (errno == EINTR)
        continue;
      cpperror("Failed to receive client argument length");
    }
    if (received == 0)
      throw std::runtime_error("Premature closure of client socket");
    break;
  }
  const size_t fdNum = lens[0];
  const size_t argvLen = lens[1];
  if (fdNum > maxFds)
    throw std::runtime_error("Client sent too many file descriptors");
  if (fdNum == 0)
    throw std::runtime_error(
        "Client sent no file descriptors"); // Should always have at least
                                            // current directory
  if (argvLen >= maxArgvLen)
    throw std::runtime_error("Client argument length too long");

  // Get transferred file descriptors from ancillary data
  for (auto cmsg = CMSG_FIRSTHDR(&msg); cmsg != nullptr;
       cmsg = CMSG_NXTHDR(&msg, cmsg)) {
    if (cmsg->cmsg_level == SOL_SOCKET && cmsg->cmsg_type == SCM_RIGHTS) {
      fds.resize(fdNum);
      memcpy(fds.data(), CMSG_DATA(cmsg), dataSize(fds));
      break;
    }
  }
  if (fds.empty())
    cpperror("Failed to receive client file descriptors");

  // Get the original position of the file descriptors
  std::vector<int> ofds(fdNum);
  switch (receiveBuffer(socket, ofds)) {
  case -1:
    cpperror("Failed to get the old file descriptor positions");
    break;
  case 0:
    cpperror("Premature closure from client");
    break;
  default:
    break;
  }

  // Last file descriptor is current directory
  if (fchdir(fds.back()))
    cpperror("Failed to change current working directories");
  close(fds.back());
  fds.pop_back();
  ofds.pop_back();
  redirectFds(fds, ofds);

  // First copy the existing argv[0]. Then copy the arguments sent over by the
  // client. Then append -i /. Add 2 '\0' to make sure there is an extra null
  // character after the end of the last string (as there is in the original
  // argv).
  //
  // XXX: This type of arithmetic is error prone!
  size_t arg0Len = strlen(argv[0]);
  rawArgv.resize(argvLen + arg0Len + 5 + 2, '\0');
  char* cur = rawArgv.data();
  strcpy(cur, argv[0]);
  cur += arg0Len + 1;
  ssize_t off = receiveBuffer(socket, cur, argvLen);
  switch (off) {
  case -1:
    cpperror("Failed to received argument from client");
    break;
  case 0:
    throw std::runtime_error("Premature closure of client socket");
    break;
  default:
    break;
  }
  cur += off;
  *cur++ = '-';
  *cur++ = 'i';
  ++cur;
  *cur++ = '/';
  ++cur;

  cur = rawArgv.data();
  while (cur < rawArgv.data() + argvLen + arg0Len + 5) {
    childArgv.push_back(cur);
    cur += strlen(cur) + 1;
  }
  argc = childArgv.size();
  argv = childArgv.data();

  return -1;
}

//
// Client side
//
// Pack the command line argument and my file descriptor to send to the
// server (child process) via the UNIX socket. Then wait to receive status
// of quantification by reading the UNIX socket.
//

// Returns all the open file descriptor (except those in the exception list).
std::vector<int> collectOpenFds(int exception) {
  DIR* dir = opendir("/dev/fd");
  if (!dir)
    cpperror("Failed to open /dev/fd");
  deferCloseDir closeDevFd(dir);

  std::vector<int> res;

  char* endptr;
  for (auto entry = readdir(dir); entry; entry = readdir(dir)) {
    if (entry->d_name[0] == '.')
      continue;
    long fd = strtol(entry->d_name, &endptr, 10);
    if (*endptr != '\0')
      continue; // Conversion error
    if (exception == fd)
      continue;
    if (dirfd(dir) == fd)
      continue;
    res.push_back(fd);
  }
  return res;
}

// Send command line arguments and stdint, stdout, stderr and current
// directory
void sendArguments(int unixSocket, const std::string& subcommand, const std::vector<std::string>& opts) {
  // Data for file descriptors
  auto fds = collectOpenFds(unixSocket);
  if (fds.size() > maxFds)
    cpperror("Too many file descriptors open");
  int fdDir = open(".", O_RDONLY);
  if (fdDir == -1)
    cpperror("Failed top open current directory");
  deferClose closeFdDir(fdDir);
  fds.push_back(fdDir);

  // Linear copies of the arguments
  size_t argvLen = subcommand.size() + 1;
  for (const auto& opt : opts)
    argvLen += opt.size() + 1;

  std::vector<char> rawArgv(argvLen, '\0');
  char* cur = rawArgv.data();
  strcpy(cur, subcommand.c_str());
  cur += subcommand.size() + 1;
  for (const auto& opt : opts) {
    strcpy(cur, opt.c_str());
    cur += opt.size() + 1;
  }

  // Collect number of fds, argvLen and fds (as auxiliary) in a msg for
  // sendmsg

  size_t lens[2] = {fds.size(), argvLen};
  struct iovec io = {.iov_base = lens, .iov_len = sizeof(lens)};
  std::vector<char> msgBuf(CMSG_SPACE(dataSize(fds)));
  struct msghdr msg;
  memset(&msg, '\0', sizeof(msg));
  msg.msg_iov = &io;
  msg.msg_iovlen = 1;
  msg.msg_control = msgBuf.data();
  msg.msg_controllen = msgBuf.size();
  struct cmsghdr* cmsg = CMSG_FIRSTHDR(&msg);
  cmsg->cmsg_level = SOL_SOCKET;
  cmsg->cmsg_type = SCM_RIGHTS;
  cmsg->cmsg_len = CMSG_LEN(dataSize(fds));
  memcpy(CMSG_DATA(cmsg), fds.data(), dataSize(fds));
  while (true) {
    ssize_t sent = sendmsg(unixSocket, &msg, 0);
    if (sent == -1) {
      if (errno == EINTR)
        continue;
      cpperror("Failed to send sizes of file descriptors and arguments");
    }
    break;
  }

  // Sends the original list of positions of the file descriptors
  if (sendBuffer(unixSocket, fds) == -1)
    cpperror("Failed to send list of file descriptors");

  // Now send the content of rawArgv
  if (sendBuffer(unixSocket, rawArgv) == -1)
    cpperror("Failed to send arguments");
}

int contactServer(const std::string& socketPath,
                  const std::string& subcommand,
                  const std::vector<std::string>& opts) {
  // Create unix socket and connect it
  struct sockaddr_un addr;
  if (socketPath.size() >= sizeof(addr.sun_path) - 1)
    cpperror("Path for server socket is too long");
  int unixSocket = socket(AF_UNIX, SOCK_STREAM, 0);
  if (unixSocket == -1)
    cpperror("Failed to open the server socket");
  deferClose closeSocket(unixSocket);
  addr.sun_family = AF_UNIX;
  strcpy(addr.sun_path, socketPath.c_str());
  if (connect(unixSocket, (struct sockaddr*)&addr, sizeof(addr)) == -1)
    cpperror("Failed to connect to server");

  sendArguments(unixSocket, subcommand, opts);

  // Wait for return status
  int status;
  while (true) {
    ssize_t received = recv(unixSocket, &status, sizeof(status), 0);
    if (received == -1) {
      if (errno == EINTR)
        continue;
      cpperror("Failed to get return status from server");
    }
    if (received == 0) {
      cpperror("Premature close from server");
    }
    break;
  }

  if (WIFEXITED(status))
    return WEXITSTATUS(status);
  if (WIFSIGNALED(status)) {
    struct sigaction act;
    memset(&act, '\0', sizeof(act));
    sigaction(WTERMSIG(status), &act, nullptr);
    kill(getpid(), WTERMSIG(status));
  }
  return EXIT_FAILURE;
}
