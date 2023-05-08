//
// RapMap - Rapid and accurate mapping of short reads to transcriptomes using
// quasi-mapping.
// Copyright (C) 2015, 2016 Rob Patro, Avi Srivastava, Hirak Sarkar
//
// This file is part of RapMap.
//
// RapMap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RapMap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RapMap.  If not, see <http://www.gnu.org/licenses/>.
//

#include "PufferFS.hpp"
#include <cstring>
#include <cstdlib>
#include <sys/stat.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <libgen.h>

namespace puffer {
namespace fs {

// Taken from
// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
bool FileExists(const char* path) {
  struct stat fileStat;
  if (stat(path, &fileStat)) {
    return false;
  }
  if (!S_ISREG(fileStat.st_mode)) {
    return false;
  }
  return true;
}

// Taken from
// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
bool DirExists(const char* path) {
  struct stat fileStat;
  if (stat(path, &fileStat)) {
    return false;
  }
  if (!S_ISDIR(fileStat.st_mode)) {
    return false;
  }
  return true;
}

int MakeDir(const char* path) { return mkdir(path, ACCESSPERMS); }

  /**
   * from https://github.com/troglobit/libite/blob/02d9d9f8c111a3f3e18131494dbd1be2617a8420/src/makepath.c
   * license below
   * mkpath - Like makepath() but takes a mode_t argument
   * @dir:  Directory to created, relative or absolute
   * @mode: A &mode_t mode to create @dir with
   *
   * Returns:
   * POSIX OK(0) on success, otherwise -1 with @errno set.
   */
  int MakePath(const char *dir) {
    struct stat sb;

    if (!dir) {
      errno = EINVAL;
      return -1;
    }

    if (!stat(dir, &sb))
      return 0;

    char* pname = strdup(dir);
    MakePath(dirname(pname));
    free(pname);

    return mkdir(dir, ACCESSPERMS);
  }

}
}

/* mkpath() -- Create all components leading up to a given directory
 *
 * Copyright (c) 2013-2016  Joachim Nilsson <troglobit@gmail.com>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */
