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

#include <jellyfish/dbg.hpp>
#include <jellyfish/time.hpp>
#include <sys/syscall.h>

namespace dbg {
  pthread_mutex_t print_t::_lock      = PTHREAD_MUTEX_INITIALIZER;
  volatile pid_t  print_t::_print_tid = 0;

#ifdef DEBUG
  Time _tic_time;
#endif

  void tic() {
#ifdef DEBUG
    _tic_time.now();
#endif
 }
  Time toc() {
#ifdef DEBUG
    Time t; 
    return t - _tic_time;
#else
    return Time::zero;
#endif
  }

#ifdef SYS_gettid
  pid_t gettid() { return (pid_t)syscall(SYS_gettid); }
#else
  pid_t gettid() { return getpid(); }
#endif

  int print_t::set_signal(int signum) {
    struct sigaction act;
    memset(&act, '\0', sizeof(act));
    act.sa_sigaction = signal_handler;
    act.sa_flags     = SA_SIGINFO;
    return sigaction(signum, &act, 0);
  }

  void print_t::signal_handler(int signum, siginfo_t *info, void *context) {
#ifdef  HAVE_SI_INT
    if(info->si_code != SI_QUEUE)
      return;
    _print_tid = info->si_int;
#endif
  }
}
