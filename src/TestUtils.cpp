/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Salmon.

    Salmon is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Salmon is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Salmon.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/

#include <boost/range/irange.hpp>
#include <fstream>
#include <memory>
#include <vector>

#include "matrix_tools.hpp"
#include "utils.hpp"

int main(int argc, char* argv[]) {
  using std::vector;
  using std::string;
  using std::ifstream;

  ifstream ifile(argv[1]);
  auto t2g = utils::readTranscriptToGeneMap(ifile);
  std::cerr << "read " << t2g.numTranscripts() << " transcripts\n";
  std::cerr << "which mapped to " << t2g.numGenes() << " genes\n";
  ifile.close();

  /*
   *  A = [0 1 0 0 0 1 1]
   *      [0 0 1 1 0 1 1]
   *      [0 0 0 1 1 0 0]
   *  counts = [ 0 4 4 8 5 3 2 ]
   */
  std::vector<std::vector<double>> A{
      {0, 1, 0, 0, 0, 1, 1}, {0, 0, 1, 1, 0, 1, 1}, {0, 0, 0, 1, 1, 0, 0}};
  std::vector<double> counts{0, 4, 4, 8, 5, 3, 2};

  std::unique_ptr<std::vector<std::vector<double>>> collapsedA{nullptr};
  std::unique_ptr<std::vector<double>> collapsedCounts{nullptr};
  matrix_tools::collapseIntoCategories(A, counts, collapsedA, collapsedCounts);

  for (size_t i = 0; i < collapsedA->size(); ++i) {
    std::cerr << "[";
    for (size_t j = 0; j < (*collapsedA)[0].size(); ++j) {
      std::cerr << (*collapsedA)[i][j] << ", ";
    }
    std::cerr << "]\n";
  }
  std::cerr << "\n[";
  for (auto c : *collapsedCounts) {
    std::cerr << c << ", ";
  }
  std::cerr << "]\n";

  for (auto i : boost::irange(0, 500)) {
    std::cerr << i << ", ";
  }
}