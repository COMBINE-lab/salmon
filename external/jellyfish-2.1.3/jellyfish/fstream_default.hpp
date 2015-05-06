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

#ifndef __JELLYFISH_FSTREAM_WITH_DEFAULT_HPP__
#define __JELLYFISH_FSTREAM_WITH_DEFAULT_HPP__

#include <iostream>
#include <fstream>

template<typename Base, std::ios_base::openmode def_mode>
class fstream_default : public Base {
  typedef Base super;
  static std::streambuf* open_file(const char* str, std::ios_base::openmode mode) {
    std::filebuf* fb = new std::filebuf;
    return fb->open(str, mode);
  }

  static std::streambuf* get_streambuf(const char* str, Base& def, 
                                       std::ios_base::openmode mode) {
    return (str != 0) ? open_file(str, mode) : def.rdbuf();
  }
  static std::streambuf* get_streambuf(const char* str, std::streambuf* buf, 
                                       std::ios_base::openmode mode) {
    return (str != 0) ? open_file(str, mode) : buf;
  }

  bool do_close;
public:
  fstream_default(const char* str, Base& def, std::ios_base::openmode mode = def_mode) :
    Base(get_streambuf(str, def, mode)), do_close(str != 0) { 
    if(Base::rdbuf() == 0)
      Base::setstate(std::ios_base::badbit);
  }
  fstream_default(const char* str, std::streambuf* def, std::ios_base::openmode mode = def_mode) :
    Base(get_streambuf(str, def, mode)), do_close(str != 0) {
    if(Base::rdbuf() == 0)
      Base::setstate(std::ios_base::badbit);
  }

  ~fstream_default() {
    if(do_close) {
      delete Base::rdbuf(0);
      do_close = false;
    }
  }
  // Close is a noop at this point as GCC 4.4 has a problem with
  // Base::rdbuf in methods (breaks strict aliasing). Beats me! I
  // think it is a false positive.
  void close() {} 
};

typedef fstream_default<std::ostream, std::ios_base::ios_base::out> ofstream_default;
typedef fstream_default<std::istream, std::ios_base::ios_base::in>  ifstream_default;

#endif // __JELLYFISH_FSTREAM_WITH_DEFAULT_HPP__
