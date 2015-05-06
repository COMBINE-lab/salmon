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

#include <jellyfish/thread_exec.hpp>

void jellyfish::thread_exec::exec(int nb_threads) {
  struct thread_info empty = {0, 0, 0};
  infos.resize(nb_threads, empty);

  for(int i = 0; i < nb_threads; i++) {
    infos[i].id   = i;
    infos[i].self = this;
    int err = pthread_create(&infos[i].thid, NULL, start_routine, &infos[i]);
    if(err)
      eraise(Error) << "Can't create thread" << err::no;
  }
}

void jellyfish::thread_exec::join() {
  for(unsigned int i = 0; i < infos.size(); i++) {
    int err = pthread_join(infos[i].thid, NULL);
    if(err)
      eraise(Error) << "Can't join thread '" << infos[i].thid << "'" << err::no;
  }
}

void *jellyfish::thread_exec::start_routine(void *_info) {
  struct thread_info *info = (struct thread_info *)_info;
  info->self->start(info->id);
  return 0;
}
