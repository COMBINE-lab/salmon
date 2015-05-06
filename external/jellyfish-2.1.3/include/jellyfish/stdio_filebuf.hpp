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

#ifndef __JELLYFISH_STDIO_FILEBUF_HPP__
#define __JELLYFISH_STDIO_FILEBUF_HPP__

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#include <streambuf>
#include <vector>

// Attempt to be (mostly) compatible with GCC ext/stdio_filebuf.h
// class. Contains code from stdio_filbuf.hpp and
// http://www.mr-edd.co.uk/blog/beginners_guide_streambuf. It is only
// meant as a quick replacement when stdio_filebuf is not available.

namespace jellyfish {
template<typename _CharT, typename _Traits = std::char_traits<_CharT> >
class stdio_filebuf : public std::basic_streambuf<_CharT, _Traits>
{
  const int                     fd_;
  FILE* const                   file_;
  const std::ios_base::openmode mode_;
  const size_t                  put_back_;
  std::vector<_CharT>           buffer_;

public:
  // Types:
  typedef _CharT			 char_type;
  typedef _Traits			 traits_type;
  typedef typename traits_type::int_type int_type;
  typedef typename traits_type::pos_type pos_type;
  typedef typename traits_type::off_type off_type;
  //  typedef std::size_t                    size_t;

  /**
   *  @param  __fd  An open file descriptor.
   *  @param  __mode  Same meaning as in a standard filebuf.
   *  @param  __size Optimal or preferred size of internal buffer,
   *                 in chars.
   *
   *  This constructor associates a file stream buffer with an open
   *  POSIX file descriptor. The file descriptor will be automatically
   *  closed when the stdio_filebuf is closed/destroyed.
   */
  stdio_filebuf(int __fd, std::ios_base::openmode __mode,
                size_t __size = static_cast<size_t>(BUFSIZ),
                size_t put_back = 1) :
    fd_(__fd),
    file_(0),
    mode_(__mode),
    put_back_(std::max(put_back, (size_t)1)),
    buffer_(std::max(__size, put_back_) + put_back_)
  {
    _CharT* end = buffer_.data() + buffer_.size();
    this->setg(end, end, end);
  }

  /**
   *  @param  __f  An open @c FILE*.
   *  @param  __mode  Same meaning as in a standard filebuf.
   *  @param  __size Optimal or preferred size of internal buffer,
   *                 in chars.  Defaults to system's @c BUFSIZ.
   *
   *  This constructor associates a file stream buffer with an open
   *  C @c FILE*.  The @c FILE* will not be automatically closed when the
   *  stdio_filebuf is closed/destroyed.
   */
  stdio_filebuf(FILE* __f, std::ios_base::openmode __mode,
                size_t __size = static_cast<size_t>(BUFSIZ),
                size_t put_back = 1) :
    fd_(-1),
    file_(__f),
    mode_(__mode),
    put_back_(std::max(put_back, (size_t)1)),
    buffer_(std::max(__size, put_back_) + put_back_)
  {
    _CharT* end = buffer_.data() + buffer_.size();
    this->setg(end, end, end);
  }

  /**
   *  Closes the external data stream if the file descriptor constructor
   *  was used.
   */
  virtual ~stdio_filebuf() {
    if(fd_ != -1)
      close(fd_);
  }

  /**
   *  @return  The underlying file descriptor.
   *
   *  Once associated with an external data stream, this function can be
   *  used to access the underlying POSIX file descriptor.  Note that
   *  there is no way for the library to track what you do with the
   *  descriptor, so be careful.
   */
  int
  fd() { return fd_ != -1 ? fd_ : fileno(file_); }

  /**
   *  @return  The underlying FILE*.
   *
   *  This function can be used to access the underlying "C" file pointer.
   *  Note that there is no way for the library to track what you do
   *  with the file, so be careful.
   */
  FILE*
  file() {
    if(file_) return file_;
    const char* str_mode;
    if(mode_ & std::ios_base::app) {
      str_mode = "a+";
    } else if(mode_ & std::ios_base::ate) {
      str_mode = "a";
    } else if(mode_ & std::ios_base::in) {
      str_mode = (mode_ & std::ios_base::out) ? "r+" : "r";
    } else if(mode_ & std::ios_base::out) {
      str_mode = "w";
    }
    return fdopen(fd_, str_mode);
  }

private:
  int_type underflow() {
    if(this->gptr() >= this->egptr()) {
      _CharT *base = buffer_.data();
      _CharT *start = base;

      if (this->eback() == base) {
        // Make arrangements for putback characters
        std::memcpy(base, this->egptr() - put_back_, put_back_ * sizeof(_CharT));
        start += put_back_;
      }

      // start is now the start of the buffer, proper.
      // Read from fptr_ in to the provided buffer
      const ssize_t n =
        (fd_ != -1) ?
        read(fd_, start, (buffer_.size() - (start - base)) * sizeof(_CharT)) :
        std::fread(start, sizeof(_CharT), buffer_.size() - (start - base), file_);
      if (n <= 0)
        return _Traits::eof();

      // Set buffer pointers
      this->setg(base, start, start + n);
    }
    return _Traits::to_int_type(*this->gptr());
  }

};
} // namespace jellyfish {
#endif // __JELLYFISH_STDIO_FILEBUF_HPP__
