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
#include <assert.h>
#include <jellyfish/allocators_mmap.hpp>

#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

#ifdef HAVE_VALGRIND
#include <valgrind.h>
// TODO: this should really come from the valgrind switch
// --redzone-size. Don't know how to get access to that yet!
static size_t redzone_size = 128;
#endif

void *allocators::mmap::realloc(size_t new_size) {
  void *new_ptr = MAP_FAILED;
  const size_t asize = new_size
#ifdef HAVE_VALGRIND
    + 2 * redzone_size
#endif
    ;

  if(ptr_ == MAP_FAILED) {
    new_ptr     = ::mmap(NULL, asize, PROT_WRITE|PROT_READ,
			 MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
  }
  // mremap is Linux specific
  // TODO: We must do something if it is not supported
#ifdef MREMAP_MAYMOVE
  else {
    new_ptr = ::mremap(ptr_, size_, new_size, MREMAP_MAYMOVE);
  }
#endif
  if(new_ptr == MAP_FAILED)
    return NULL;

#ifdef HAVE_VALGRIND
  new_ptr = (char*)new_ptr + redzone_size;
  if(ptr_ == MAP_FAILED)
    VALGRIND_MALLOCLIKE_BLOCK(new_ptr, new_size, redzone_size, 1);
  // TODO: resize not yet supported
#endif


  size_ = new_size;
  ptr_  = new_ptr;

  fast_zero();

  return ptr_;
}

size_t allocators::mmap::round_to_page(size_t _size) {
  static const long                                   pg_size = sysconf(_SC_PAGESIZE);
  return (_size / pg_size + (_size % pg_size != 0)) * pg_size;
}

void allocators::mmap::fast_zero() {
  tinfo  info[nb_threads];
  size_t pgsize        = round_to_page(1);
  size_t nb_pages      = size_ / pgsize + (size_ % pgsize != 0);
  int    total_threads = 0;

  for(size_t i = 0; i < (size_t)nb_threads; ++i, ++total_threads) {
    info[i].start = (char *)ptr_ + pgsize * ((i * nb_pages) / nb_threads);
    info[i].end   = (char *)ptr_ + std::min(pgsize * (((i + 1) * nb_pages) / nb_threads), size_);
    info[i].pgsize = pgsize;
    if(pthread_create(&info[i].thid, NULL, _fast_zero, &info[i]))
      break;
  }

  for(int i = 0; i < total_threads; i++)
    pthread_join(info[i].thid, NULL);
}

void * allocators::mmap::_fast_zero(void *_info) {
  tinfo *info = (tinfo *)_info;

  for(char *cptr = info->start; cptr < info->end; cptr += info->pgsize) {
    *cptr = 0;
  }

  return NULL;
}

void allocators::mmap::free() {
  if(ptr_ == MAP_FAILED)
    return;
#ifdef HAVE_VALGRIND
  VALGRIND_FREELIKE_BLOCK(ptr_, redzone_size);
  ptr_   = (char*)ptr_ - redzone_size;
  size_ += 2 * redzone_size;
#endif
  assert(::munmap(ptr_, size_) == 0);
  ptr_  = MAP_FAILED;
  size_ = 0;
}
