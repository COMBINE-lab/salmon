#ifndef __SALMON_SERVER_HPP__
#define __SALMON_SERVER_HPP__

#include <memory>

#include "SalmonIndex.hpp"

typedef const char** argvtype;

// If no --server switch is in the arguments, this function is a NOOP.
//
// Otherwise, if both --server and --index,-i are given, it behaves like a
// server: it loads the index, opens a server unix socket, and waits for
// requests. Upon a requests, it forks, redirect stdout and stderr to what the
// client provides, override argc and argv to what the client provides and
// returns to continue processing that request normally (in a sub-process).
//
// If only --server is passed, it behaves like a client. It opens a client unix
// socket, sends it the stdout and stderr, and the content of argc and argv
// (minus --server), and then wait until the socket is closed and exits. The
// client function never returns as the processing is done by a process forked
// by the server.
//
// Returns -1 if work should continue normally (either no server, or inside a
// child process). Otherwise, it returns an error status.
// int salmonServer(int& argc, argvtype& argv,
//                  std::unique_ptr<SalmonIndex>& salmonIndex);

#endif // __SALMON_SERVER_HPP__
