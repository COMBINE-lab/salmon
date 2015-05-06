/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
n    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <signal.h>
#include <config.h>
#include <jellyfish/misc.hpp>
#include <jellyfish/err.hpp>
#include <jellyfish/backtrace.hpp>

namespace jellyfish {
uint64_t bogus_sum(void *data, size_t len) {
  uint64_t res = 0, tmp = 0;
  uint64_t *ptr = (uint64_t *)data;

  while(len >= sizeof(uint64_t)) {
    res ^= *ptr++;
    len -= sizeof(uint64_t);
  }
  
  if(len > 0) {
    memcpy(&tmp, ptr, len);
    res ^= tmp;
  }
  return res;
}

void disabled_misaligned_mem_access() {
#if defined(__GNUC__)
# if defined(__i386__)
  /* Enable Alignment Checking on x86 */ 
  __asm__("pushf\norl $0x40000,(%esp)\npopf");
# elif defined(__x86_64__)
  /* Enable Alignment Checking on x86_64 */
  __asm__("pushf\norl $0x40000,(%rsp)\npopf");
# endif
#endif 
}

// Return -1 if size cannot be obtained.
std::streamoff get_file_size(std::istream& is) {
  if(!is.good()) return -1;
  std::streampos cpos = is.tellg();
  if(!is.good()) { is.clear(); return -1; }
  is.seekg(0, std::ios::end);
  if(!is.good()) { is.clear(); return -1; }
  std::streamoff res = is.tellg() - cpos;
  if(!is.good()) { is.clear(); return -1; }
  is.seekg(cpos);
  return res;
}

template<long int n>
struct ConstFloorLog2 {
  static const int val = ConstFloorLog2<n / 2>::val + 1;
};
template<>
struct ConstFloorLog2<1> {
  static const int val = 0;
};
// Return length random bits
uint64_t random_bits(int length) {
  uint64_t res = 0;
  for(int i = 0; i < length; i += ConstFloorLog2<RAND_MAX>::val) {
    res ^= (uint64_t)random() << i;
  }
  return res & ((uint64_t)-1 >> (bsizeof(uint64_t) - length));
}

bool isblunt(char c) {
  return isalnum(c) || c == '_' || c == '-' || c == '/' || c == '.';
}
std::string quote_arg(const std::string& arg) {
  if(std::all_of(arg.begin(), arg.end(), isblunt))
    return arg;

  std::string res("'");
  size_t pos = 0;
  while(true) {
    size_t qpos = arg.find_first_of("'", pos);
    res += arg.substr(pos, qpos - pos);
    if(qpos == std::string::npos) break;
    res += "'\\''";
    pos = qpos + 1;
  }
  res += "'";
  return res;
}

} // namespace jellyfish
