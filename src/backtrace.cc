/*  This file is part of Jellyfish.

    This work is dual-licensed under 3-Clause BSD License or GPL 3.0.
    You can choose between one of them if you use this work.

`SPDX-License-Identifier: BSD-3-Clause OR  GPL-3.0`
*/

#include <iostream>
#include <stdexcept>
#include <execinfo.h>
#include <stdlib.h>
#include <typeinfo>
#include <cxxabi.h>

void print_backtrace() {
  void  *trace_elems[20];
  int    trace_elem_count(backtrace(trace_elems, 20));
  backtrace_symbols_fd(trace_elems, trace_elem_count, 2);
}

static void handler() {
  // Display message of last thrown exception if any
  try { throw; }
  catch(const std::exception& e) {
    int     status;
    size_t  n = 0;
    char   *name = abi::__cxa_demangle(typeid(e).name(), 0, &n, &status);
    std::cerr << "terminate called after throwing an instance of '"
              << (status < 0 ? "UNKNOWN" : name)
              << "'\n  what(): " << e.what() << "\n";
    if(n)
      free(name);
  }
  catch(...) {}

  print_backtrace();
  abort();
}

void show_backtrace() {
  std::set_terminate(handler);
}
