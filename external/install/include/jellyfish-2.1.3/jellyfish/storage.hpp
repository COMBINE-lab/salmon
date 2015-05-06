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

#ifndef __JELLYFISH_STORAGE_HPP__
#define __JELLYFISH_STORAGE_HPP__

#include <stdlib.h>
#include <stdint.h>
#include <jellyfish/misc.hpp>

namespace jellyfish {

  class storage_t {
  public:
    storage_t() {}
    virtual ~storage_t() {}
  };

  // Entry 0 is used only when switching to a large field
  extern const size_t *quadratic_reprobes;

}

#endif // __STORAGE_HPP__
