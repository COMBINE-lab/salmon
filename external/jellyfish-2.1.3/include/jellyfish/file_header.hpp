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

#ifndef __JELLYFISH_FILE_HEADER_HPP__
#define __JELLYFISH_FILE_HEADER_HPP__

#include <string>
#include <jellyfish/generic_file_header.hpp>
#include <jellyfish/rectangular_binary_matrix.hpp>

namespace jellyfish {
/// A header with jellyfish hash specific entries: size, matrix, etc.
class file_header : public generic_file_header {
public:
  file_header() : generic_file_header(sizeof(uint64_t)) { }
  file_header(std::istream& is) : generic_file_header(sizeof(uint64_t)) {
    this->read(is);
  }

  template<typename storage>
  void update_from_ary(const storage& ary) {
    this->size(ary.size());
    this->key_len(ary.key_len());
    this->val_len(ary.val_len());
    this->matrix(ary.matrix());
    this->max_reprobe(ary.max_reprobe());
    this->set_reprobes(ary.reprobes());
  }

  RectangularBinaryMatrix matrix(int i = 1) const {
    std::string name("matrix");
    name += std::to_string((long long int)i); // Cast to make gcc4.4 happy!
    const unsigned int r = root_[name]["r"].asUInt();
    const unsigned int c = root_[name]["c"].asUInt();
    std::vector<uint64_t> raw(c, (uint64_t)0);
    for(unsigned int i = 0; i < c; ++i)
      raw[i] = root_[name]["columns"][i].asUInt64();
    return RectangularBinaryMatrix(raw.data(), r, c);
  }

  void matrix(const RectangularBinaryMatrix& m, int i = 1) {
    std::string name("matrix");
    name += std::to_string((long long int)i);
    root_[name].clear();
    root_[name]["r"] = m.r();
    root_[name]["c"] = m.c();
    for(unsigned int i = 0; i < m.c(); ++i) {
      Json::UInt64 x = m[i];
      root_[name]["columns"].append(x);
    }
  }

  size_t size() const { return root_["size"].asLargestUInt(); }
  void size(size_t s) { root_["size"] = (Json::UInt64)s; }

  unsigned int key_len() const { return root_["key_len"].asUInt(); }
  void key_len(unsigned int k) { root_["key_len"] = (Json::UInt)k; }

  unsigned int val_len() const { return root_["val_len"].asUInt(); }
  void val_len(unsigned int k) { root_["val_len"] = (Json::UInt)k; }

  unsigned int max_reprobe() const { return root_["max_reprobe"].asUInt(); }
  void max_reprobe(unsigned int m) { root_["max_reprobe"] = (Json::UInt)m; }

  size_t max_reprobe_offset() const { return root_["reprobes"][max_reprobe()].asLargestUInt(); }

  double fpr() const { return root_["fpr"].asDouble(); }
  void fpr(double f) { root_["fpr"] = f; }

  unsigned long nb_hashes() const { return root_["nb_hashes"].asUInt(); }
  void nb_hashes(unsigned long nbh) { root_["nb_hashes"] = (Json::UInt)nbh; }

  bool canonical() const { return root_.get("canonical", false).asBool(); }
  void canonical(bool v) { root_["canonical"] = v; }

  /// reprobes must be at least max_reprobe() + 1 long
  void get_reprobes(size_t* reprobes) const {
    for(unsigned int i = 0; i <= max_reprobe(); ++i)
      reprobes[i] = root_["reprobes"][i].asLargestUInt();
  }

  /// This must be call after max_reprobe has been set. reprobes must
  /// be at least max_reprobe() + 1 long.
  void set_reprobes(const size_t* reprobes) {
    root_["reprobes"].clear();
    for(unsigned int i = 0; i <= max_reprobe(); ++i)
      root_["reprobes"].append((Json::UInt64)reprobes[i]);
  }

  /// Length of counter field in binary/sorted format
  unsigned int counter_len() const { return root_["counter_len"].asUInt(); }
  void counter_len(unsigned int l) { root_["counter_len"] = (Json::UInt)l; }

  std::string format() const { return root_["format"].asString(); }
  void format(const std::string& s) { root_["format"] = s; }
};
} // namespace jellyfish

#endif /* __JELLYFISH_FILE_HEADER_HPP__ */
