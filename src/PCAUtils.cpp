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

#include "pca.h"
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace stats {
namespace utils {

arma::Mat<double> make_covariance_matrix(const arma::Mat<double>& data) {
  return std::move((data.t() * data) * (1. / (data.n_rows - 1)));
}

arma::Mat<double> make_shuffled_matrix(const arma::Mat<double>& data) {
  const long n_rows = data.n_rows;
  const long n_cols = data.n_cols;
  arma::Mat<double> shuffle(n_rows, n_cols);
  for (long j = 0; j < n_cols; ++j) {
    for (long i = 0; i < n_rows; ++i) {
      shuffle(i, j) = data(std::rand() % n_rows, j);
    }
  }
  return std::move(shuffle);
}

arma::Col<double> compute_column_means(const arma::Mat<double>& data) {
  const long n_cols = data.n_cols;
  arma::Col<double> means(n_cols);
  for (long i = 0; i < n_cols; ++i)
    means(i) = arma::mean(data.col(i));
  return std::move(means);
}

void remove_column_means(arma::Mat<double>& data,
                         const arma::Col<double>& means) {
  if (data.n_cols != means.n_elem)
    throw std::range_error("Number of elements of means is not equal to the "
                           "number of columns of data");
  for (long i = 0; i < long(data.n_cols); ++i)
    data.col(i) -= means(i);
}

arma::Col<double> compute_column_rms(const arma::Mat<double>& data) {
  const long n_cols = data.n_cols;
  arma::Col<double> rms(n_cols);
  for (long i = 0; i < n_cols; ++i) {
    const double dot = arma::dot(data.col(i), data.col(i));
    rms(i) = std::sqrt(dot / (data.col(i).n_rows - 1));
  }
  return std::move(rms);
}

void normalize_by_column(arma::Mat<double>& data,
                         const arma::Col<double>& rms) {
  if (data.n_cols != rms.n_elem)
    throw std::range_error("Number of elements of rms is not equal to the "
                           "number of columns of data");
  for (long i = 0; i < long(data.n_cols); ++i) {
    if (rms(i) == 0)
      throw std::runtime_error(
          "At least one of the entries of rms equals to zero");
    data.col(i) *= 1. / rms(i);
  }
}

void enforce_positive_sign_by_column(arma::Mat<double>& data) {
  for (long i = 0; i < long(data.n_cols); ++i) {
    const double max = arma::max(data.col(i));
    const double min = arma::min(data.col(i));
    bool change_sign = false;
    if (std::abs(max) >= std::abs(min)) {
      if (max < 0)
        change_sign = true;
    } else {
      if (min < 0)
        change_sign = true;
    }
    if (change_sign)
      data.col(i) *= -1;
  }
}

std::vector<double> extract_column_vector(const arma::Mat<double>& data,
                                          long index) {
  if (index < 0 || index >= long(data.n_cols))
    throw std::range_error(join("Index out of range: ", index));
  const long n_rows = data.n_rows;
  const double* memptr = data.colptr(index);
  std::vector<double> result(memptr, memptr + n_rows);
  return std::move(result);
}

std::vector<double> extract_row_vector(const arma::Mat<double>& data,
                                       long index) {
  if (index < 0 || index >= long(data.n_rows))
    throw std::range_error(join("Index out of range: ", index));
  const arma::Row<double> row(data.row(index));
  const double* memptr = row.memptr();
  std::vector<double> result(memptr, memptr + row.n_elem);
  return std::move(result);
}

void assert_file_good(const bool& is_file_good, const std::string& filename) {
  if (!is_file_good)
    throw std::ios_base::failure(join("Cannot open file: ", filename));
}

double get_mean(const std::vector<double>& iter) {
  const double init = 0;
  return std::accumulate(iter.begin(), iter.end(), init) / iter.size();
}

double get_sigma(const std::vector<double>& iter) {
  const double mean = get_mean(iter);
  double sum = 0;
  for (auto v = iter.begin(); v != iter.end(); ++v)
    sum += std::pow(*v - mean, 2.);
  return std::sqrt(sum / (iter.size() - 1));
}

} // namespace utils
} // namespace stats
