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

#ifndef _INT128_H_
#define _INT128_H_

//#ifdef HAVE_CONFIG_H
//#include <config.h>
//#ifndef HAVE_INT128
//#error "The type __int128 is not supported"
//#endif
//#endif

#include <unistd.h>
#include <stdint.h>
#include <iostream>
#include <sstream>
#include <limits>
#include <cstring>

// Output of __int128: this might be slow
namespace __int128_ns {
template<int base>
void __print_digits(std::ostream& os, unsigned __int128 x,
                    bool lower = true) {
  char buf[50];
  char* ptr = buf + sizeof(buf);
  do {
    int o  = x % base;
    if(o < 10) {
      *--ptr = '0' + o;
    } else {
      *--ptr = (lower ? 'a' : 'A') + o - 10;
    }
    x     /= base;
  } while (x > 0);
  os.write(ptr, buf + sizeof(buf) - ptr);
}

inline bool is_negative(unsigned __int128 /*x*/) { return false; }
inline bool is_negative(__int128 x) { return x < 0; }

template<typename T>
void __print_decimal(std::ostream& prefix, std::ostream& os, T x,
                     const std::ios::fmtflags& ff) {
  if((ff & std::ios::showpos) && x > 0)
    prefix << "+";
  if(x == 0) {
    os << "0";
    return;
  }
  if(is_negative(x)) {
    prefix << "-";
    x = -x;
  }
  __print_digits<10>(os, x);
}

void __print_bases(std::ostream& prefix, std::ostream& os,
                   unsigned __int128 x, 
                   const std::ios::fmtflags& ff);

template<typename T>
void __print_buf(std::ostream& prefix, std::ostream& os, T x,
                 const std::ios::fmtflags& ff) {
  if(ff & std::ios::dec)
    __print_decimal(prefix, os, x, ff);
  else
    __print_bases(prefix, os, (unsigned __int128)x, ff);
}

template<typename T>
void __print(std::ostream&os, T x) {
  const std::ios_base::fmtflags ff = os.flags();

  if(!(ff & std::ios::adjustfield))
    return __print_buf(os, os, x, ff);

  std::ostringstream prefix;
  std::ostringstream buf;
  __print_buf(prefix, buf, x, ff);
  ssize_t nb_padding = os.width() - (prefix.str().size() + buf.str().size());
  if(nb_padding <= 0) {
    os.write(prefix.str().c_str(), prefix.tellp());
    os.write(buf.str().c_str(), buf.tellp());
    return;
  }

  char padding[nb_padding];
  memset(padding, os.fill(), nb_padding);
  if(ff & std::ios::right)
    os.write(padding, nb_padding);
  os.write(prefix.str().c_str(), prefix.tellp());
  if(ff & std::ios::internal)
    os.write(padding, nb_padding);
  os.write(buf.str().c_str(), buf.tellp());
  if(ff & std::ios::left)
    os.write(padding, nb_padding);
}
}

inline
std::ostream& operator<<(std::ostream& os, __int128 x) {
  __int128_ns::__print(os, x);
  return os;
}

inline
std::ostream& operator<<(std::ostream& os, unsigned __int128 x) {
  __int128_ns::__print(os, x);
  return os;
}

#ifndef HAVE_NUMERIC_LIMITS128
namespace std {
template<>
class numeric_limits<__int128> {
public:
  static const bool is_specialized = true;
  static __int128 max() { return (unsigned __int128)-1 >> 1; }
  static __int128 min() { return  max() + 1; }
  static const int  digits     = 127;
  static const int  digits10   = 38;
#define NLS64 numeric_limits<int64_t>
  static const bool is_signed  = NLS64::is_signed;
  static const bool is_integer = NLS64::is_integer;
  static const bool is_exact   = NLS64::is_exact;
  static const int  radix      = NLS64::radix;
  static __int128 epsilon() { return NLS64::epsilon(); }
  static __int128 round_error() { return NLS64::round_error(); }
  static const int                min_exponent      = NLS64::min_exponent;
  static const int                min_exponent10    = NLS64::min_exponent10;
  static const int                max_exponent      = NLS64::max_exponent;
  static const int                max_exponent10    = NLS64::max_exponent10;
  static const bool               has_infinity      = NLS64::has_infinity;
  static const bool               has_quiet_NaN     = NLS64::has_quiet_NaN;
  static const bool               has_signaling_NaN = NLS64::has_signaling_NaN;
  static const float_denorm_style has_denorm        = NLS64::has_denorm;
  static const bool               has_denorm_loss   = NLS64::has_denorm_loss;
  static __int128 infinity() { return NLS64::infinity(); }
  static __int128 quiet_NaN() { return NLS64::quiet_NaN(); }
  static __int128 signaling_NaN() { return NLS64::signaling_NaN(); }
  static __int128 denorm_min() { return NLS64::denorm_min(); }
  static const bool              is_iec559       = NLS64::is_iec559;
  static const bool              is_bounded      = NLS64::is_bounded;
  static const bool              is_modulo       = NLS64::is_modulo;
  static const bool              traps           = NLS64::traps;
  static const bool              tinyness_before = NLS64::tinyness_before;
  static const float_round_style round_style     = NLS64::round_style;
};

template<>
class numeric_limits<unsigned __int128> {
public:
  static const bool is_specialized = true;
  static __int128 max() { return (unsigned __int128)-1; }
  static __int128 min() { return  0; }
  static const int  digits     = 128;
  static const int  digits10   = 39;
#define NLU64 numeric_limits<uint64_t>
  static const bool is_signed  = NLU64::is_signed;
  static const bool is_integer = NLU64::is_integer;
  static const bool is_exact   = NLU64::is_exact;
  static const int  radix      = NLU64::radix;
  static __int128 epsilon() { return NLU64::epsilon(); }
  static __int128 round_error() { return NLU64::round_error(); }
  static const int                min_exponent      = NLU64::min_exponent;
  static const int                min_exponent10    = NLU64::min_exponent10;
  static const int                max_exponent      = NLU64::max_exponent;
  static const int                max_exponent10    = NLU64::max_exponent10;
  static const bool               has_infinity      = NLU64::has_infinity;
  static const bool               has_quiet_NaN     = NLU64::has_quiet_NaN;
  static const bool               has_signaling_NaN = NLU64::has_signaling_NaN;
  static const float_denorm_style has_denorm        = NLU64::has_denorm;
  static const bool               has_denorm_loss   = NLU64::has_denorm_loss;
  static __int128 infinity() { return NLU64::infinity(); }
  static __int128 quiet_NaN() { return NLU64::quiet_NaN(); }
  static __int128 signaling_NaN() { return NLU64::signaling_NaN(); }
  static __int128 denorm_min() { return NLU64::denorm_min(); }
  static const bool              is_iec559       = NLU64::is_iec559;
  static const bool              is_bounded      = NLU64::is_bounded;
  static const bool              is_modulo       = NLU64::is_modulo;
  static const bool              traps           = NLU64::traps;
  static const bool              tinyness_before = NLU64::tinyness_before;
  static const float_round_style round_style     = NLU64::round_style;
};
} // namespace std
#endif /* HAVE_NUMERIC_LIMITS128 */

#endif /* _INT128_H_ */
