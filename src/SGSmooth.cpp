// From :
// https://raw.githubusercontent.com/thatchristoph/vmd-cvs-github/master/plugins/signalproc/src/sgsmooth.C
//!
// Sliding window signal processing (and linear algebra toolkit).
//
// supported operations:
// <ul>
// <li> Savitzky-Golay smoothing.
// <li> computing a numerical derivative based of Savitzky-Golay smoothing.
// <li> required linear algebra support for SG smoothing using STL based
//      vector/matrix classes
// </ul>
//
// \brief Linear Algebra "Toolkit".
//
// modified by Rob Patro, 2016

// system headers
#include <cmath>   // for fabs
#include <cstddef> // for size_t
#include <cstdio>
#include <vector>

//! default convergence
static const double TINY_FLOAT = 1.0e-300;

//! comfortable array of doubles
using float_vect = std::vector<double>;
//! comfortable array of ints;
using int_vect = std::vector<int>;

/*! matrix class.
 *
 * This is a matrix class derived from a vector of float_vects.  Note that
 * the matrix elements indexed [row][column] with indices starting at 0 (c
 * style). Also note that because of its design looping through rows should
 * be faster than looping through columns.
 *
 * \brief two dimensional floating point array
 */
class float_mat : public std::vector<float_vect> {
private:
  //! disable the default constructor
  explicit float_mat(){};
  //! disable assignment operator until it is implemented.
  float_mat& operator=(const float_mat&) { return *this; };

public:
  //! constructor with sizes
  float_mat(const size_t rows, const size_t cols, const double def = 0.0);
  //! copy constructor for matrix
  float_mat(const float_mat& m);
  //! copy constructor for vector
  float_mat(const float_vect& v);

  //! use default destructor
  // ~float_mat() {};

  //! get size
  size_t nr_rows(void) const { return size(); };
  //! get size
  size_t nr_cols(void) const { return front().size(); };
};

// constructor with sizes
float_mat::float_mat(const size_t rows, const size_t cols, const double defval)
    : std::vector<float_vect>(rows) {
    int sizedRows = static_cast<int>(rows);
    int sizedCols = static_cast<int>(cols);
  int i;
  for (i = 0; i < sizedRows; ++i) {
    (*this)[i].resize(cols, defval);
  }
  if ((rows < 1) || (cols < 1)) {
    char buffer[1024];

    sprintf(buffer, "cannot build matrix with %d rows and %d columns\n", sizedRows, sizedCols);
    // sgs_error(buffer);
  }
}

// copy constructor for matrix
float_mat::float_mat(const float_mat& m) : std::vector<float_vect>(m.size()) {

  float_mat::iterator inew = begin();
  float_mat::const_iterator iold = m.begin();
  for (/* empty */; iold < m.end(); ++inew, ++iold) {
    const size_t oldsz = iold->size();
    inew->resize(oldsz);
    const float_vect oldvec(*iold);
    *inew = oldvec;
  }
}

// copy constructor for vector
float_mat::float_mat(const float_vect& v) : std::vector<float_vect>(1) {

  const size_t oldsz = v.size();
  front().resize(oldsz);
  front() = v;
}

//////////////////////
// Helper functions //
//////////////////////

//! permute() orders the rows of A to match the integers in the index array.
void permute(float_mat& A, int_vect& idx) {
  int_vect i(idx.size());
  int j, k;
  int nrows = static_cast<int>(A.nr_rows());
  for (j = 0; j < nrows; ++j) {
    i[j] = j;
  }

  // loop over permuted indices
  for (j = 0; j < nrows; ++j) {
    if (i[j] != idx[j]) {

      // search only the remaining indices
      for (k = j + 1; k < nrows; ++k) {
        if (i[k] == idx[j]) {
          std::swap(A[j], A[k]); // swap the rows and
          i[k] = i[j];           // the elements of
          i[j] = idx[j];         // the ordered index.
          break;                 // next j
        }
      }
    }
  }
}

/*! \brief Implicit partial pivoting.
 *
 * The function looks for pivot element only in rows below the current
 * element, A[idx[row]][column], then swaps that row with the current one in
 * the index map. The algorithm is for implicit pivoting (i.e., the pivot is
 * chosen as if the max coefficient in each row is set to 1) based on the
 * scaling information in the vector scale. The map of swapped indices is
 * recorded in swp. The return value is +1 or -1 depending on whether the
 * number of row swaps was even or odd respectively. */
static int partial_pivot(float_mat& A, const size_t row, const size_t col,
                         float_vect& scale, int_vect& idx, double tol) {
  if (tol <= 0.0)
    tol = TINY_FLOAT;

  int swapNum = 1;

  // default pivot is the current position, [row,col]
  int pivot = row;
  double piv_elem = fabs(A[idx[row]][col]) * scale[idx[row]];

  // loop over possible pivots below current
  int j;
  int nrows = static_cast<int>(A.nr_rows());
  for (j = row + 1; j < nrows; ++j) {

    const double tmp = fabs(A[idx[j]][col]) * scale[idx[j]];

    // if this elem is larger, then it becomes the pivot
    if (tmp > piv_elem) {
      pivot = j;
      piv_elem = tmp;
    }
  }

#if 0
    if(piv_elem < tol) {
      //sgs_error("partial_pivot(): Zero pivot encountered.\n")
#endif
  int srow = static_cast<int>(row);
  if (pivot > srow) { // bring the pivot to the diagonal
    j = idx[row];    // reorder swap array
    idx[row] = idx[pivot];
    idx[pivot] = j;
    swapNum = -swapNum; // keeping track of odd or even swap
  }
  return swapNum;
}

/*! \brief Perform backward substitution.
 *
 * Solves the system of equations A*b=a, ASSUMING that A is upper
 * triangular. If diag==1, then the diagonal elements are additionally
 * assumed to be 1.  Note that the lower triangular elements are never
 * checked, so this function is valid to use after a LU-decomposition in
 * place.  A is not modified, and the solution, b, is returned in a. */
static void lu_backsubst(float_mat& A, float_mat& a, bool diag = false) {
  int r, c, k;
  int nrows = static_cast<int>(A.nr_rows());
  int ncols = static_cast<int>(A.nr_cols());
  for (r = (nrows - 1); r >= 0; --r) {
    for (c = (ncols - 1); c > r; --c) {
      for (k = 0; k < ncols; ++k) {
        a[r][k] -= A[r][c] * a[c][k];
      }
    }
    if (!diag) {
      for (k = 0; k < ncols; ++k) {
        a[r][k] /= A[r][r];
      }
    }
  }
}

/*! \brief Perform forward substitution.
 *
 * Solves the system of equations A*b=a, ASSUMING that A is lower
 * triangular. If diag==1, then the diagonal elements are additionally
 * assumed to be 1.  Note that the upper triangular elements are never
 * checked, so this function is valid to use after a LU-decomposition in
 * place.  A is not modified, and the solution, b, is returned in a. */
static void lu_forwsubst(float_mat& A, float_mat& a, bool diag = true) {
  int r, k, c;
  int nrows = static_cast<int>(A.nr_rows());
  int ncols = static_cast<int>(A.nr_cols());
  for (r = 0; r < nrows; ++r) {
    for (c = 0; c < r; ++c) {
      for (k = 0; k < ncols; ++k) {
        a[r][k] -= A[r][c] * a[c][k];
      }
    }
    if (!diag) {
      for (k = 0; k < ncols; ++k) {
        a[r][k] /= A[r][r];
      }
    }
  }
}

/*! \brief Performs LU factorization in place.
 *
 * This is Crout's algorithm (cf., Num. Rec. in C, Section 2.3).  The map of
 * swapped indeces is recorded in idx. The return value is +1 or -1
 * depending on whether the number of row swaps was even or odd
 * respectively.  idx must be preinitialized to a valid set of indices
 * (e.g., {1,2, ... ,A.nr_rows()}). */
static int lu_factorize(float_mat& A, int_vect& idx, double tol = TINY_FLOAT) {
  if (tol <= 0.0)
    tol = TINY_FLOAT;

  if ((A.nr_rows() == 0) || (A.nr_rows() != A.nr_cols())) {
    // sgs_error("lu_factorize(): cannot handle empty "
    //           "or nonsquare matrices.\n");

    return 0;
  }

  int nrows = static_cast<int>(A.nr_rows());
  int ncols = static_cast<int>(A.nr_cols());

  float_vect scale(A.nr_rows()); // implicit pivot scaling
  int i, j;
  for (i = 0; i < nrows; ++i) {
    double maxval = 0.0;
    for (j = 0; j < ncols; ++j) {
      if (fabs(A[i][j]) > maxval)
        maxval = fabs(A[i][j]);
    }
    if (maxval == 0.0) {
      // sgs_error("lu_factorize(): zero pivot found.\n");
      return 0;
    }
    scale[i] = 1.0 / maxval;
  }

  int swapNum = 1;
  int c, r;
  for (c = 0; c < ncols; ++c) { // loop over columns
    swapNum *=
        partial_pivot(A, c, c, scale, idx, tol); // bring pivot to diagonal
    for (r = 0; r < nrows; ++r) {          //  loop over rows
      int lim = (r < c) ? r : c;
      for (j = 0; j < lim; ++j) {
        A[idx[r]][c] -= A[idx[r]][j] * A[idx[j]][c];
      }
      if (r > c)
        A[idx[r]][c] /= A[idx[c]][c];
    }
  }
  permute(A, idx);
  return swapNum;
}

/*! \brief Solve a system of linear equations.
 * Solves the inhomogeneous matrix problem with lu-decomposition. Note that
 * inversion may be accomplished by setting a to the identity_matrix. */
static float_mat lin_solve(const float_mat& A, const float_mat& a,
                           double tol = TINY_FLOAT) {
  float_mat B(A);
  float_mat b(a);
  int_vect idx(B.nr_rows());
  int j;
  int nrows = static_cast<int>(B.nr_rows());
  int ncols = static_cast<int>(B.nr_cols());

  for (j = 0; j < nrows; ++j) {
    idx[j] = j; // init row swap label array
  }
  lu_factorize(B, idx, tol); // get the lu-decomp.
  permute(b, idx);           // sort the inhomogeneity to match the lu-decomp
  lu_forwsubst(B, b);        // solve the forward problem
  lu_backsubst(B, b);        // solve the backward problem
  return b;
}

///////////////////////
// related functions //
///////////////////////

//! Returns the inverse of a matrix using LU-decomposition.
static float_mat invert(const float_mat& A) {
  const int n = A.size();
  float_mat E(n, n, 0.0);
  float_mat B(A);
  int i;

  for (i = 0; i < n; ++i) {
    E[i][i] = 1.0;
  }

  return lin_solve(B, E);
}

//! returns the transposed matrix.
static float_mat transpose(const float_mat& a) {
  float_mat res(a.nr_cols(), a.nr_rows());
  int i, j;
  int nrows = static_cast<int>(a.nr_rows());
  int ncols = static_cast<int>(a.nr_cols());

  for (i = 0; i < nrows; ++i) {
    for (j = 0; j < ncols; ++j) {
      res[j][i] = a[i][j];
    }
  }
  return res;
}

//! matrix multiplication.
float_mat operator*(const float_mat& a, const float_mat& b) {
  float_mat res(a.nr_rows(), b.nr_cols());
  if (a.nr_cols() != b.nr_rows()) {
    // sgs_error("incompatible matrices in multiplication\n");
    return res;
  }

  int i, j, k;
  int arows = static_cast<int>(a.nr_rows());
  int acols = static_cast<int>(a.nr_cols());
  int bcols = static_cast<int>(b.nr_cols());
  for (i = 0; i < arows; ++i) {
    for (j = 0; j < bcols; ++j) {
      double sum(0.0);
      for (k = 0; k < acols; ++k) {
        sum += a[i][k] * b[k][j];
      }
      res[i][j] = sum;
    }
  }
  return res;
}

//! calculate savitzky golay coefficients.
static float_vect sg_coeff(const float_vect& b, const size_t deg) {
  const size_t rows(b.size());
  const size_t cols(deg + 1);
  float_mat A(rows, cols);
  float_vect res(rows);

  // generate input matrix for least squares fit
  int i, j;
  int srows = static_cast<int>(rows);
  int scols = static_cast<int>(cols);
  int sdeg = static_cast<int>(deg);
  for (i = 0; i < srows; ++i) {
    for (j = 0; j < scols; ++j) {
      A[i][j] = pow(double(i), double(j));
    }
  }

  float_mat c(invert(transpose(A) * A) * (transpose(A) * transpose(b)));
  int bsize = static_cast<int>(b.size());
  for (i = 0; i < bsize; ++i) {
    res[i] = c[0][0];
    for (j = 1; j <= sdeg; ++j) {
      res[i] += c[j][0] * pow(double(i), double(j));
    }
  }
  return res;
}

/*! \brief savitzky golay smoothing.
 *
 * This method means fitting a polynome of degree 'deg' to a sliding window
 * of width 2w+1 throughout the data.  The needed coefficients are
 * generated dynamically by doing a least squares fit on a "symmetric" unit
 * vector of size 2w+1, e.g. for w=2 b=(0,0,1,0,0). evaluating the polynome
 * yields the sg-coefficients.  at the border non symmectric vectors b are
 * used. */
float_vect sg_smooth(const float_vect& v, const int width, const int deg) {
  float_vect res(v.size(), 0.0);
  int vsize = static_cast<int>(v.size());
  if ((width < 1) || (deg < 0) || (vsize < (2 * width + 2))) {
    // sgs_error("sgsmooth: parameter error.\n");
    return res;
  }

  const int window = 2 * width + 1;
  const int endidx = vsize - 1;

  // do a regular sliding window average
  int i, j;
  if (deg == 0) {
// handle border cases first because we need different coefficients
#if defined(_OPENMP)
#pragma omp parallel for private(i, j) schedule(static)
#endif
    for (i = 0; i < width; ++i) {
      const double scale = 1.0 / double(i + 1);
      const float_vect c1(width, scale);
      for (j = 0; j <= i; ++j) {
        res[i] += c1[j] * v[j];
        res[endidx - i] += c1[j] * v[endidx - j];
      }
    }

    // now loop over rest of data. reusing the "symmetric" coefficients.
    const double scale = 1.0 / double(window);
    const float_vect c2(window, scale);
#if defined(_OPENMP)
#pragma omp parallel for private(i, j) schedule(static)
#endif
    for (i = 0; i <= (vsize - window); ++i) {
      for (j = 0; j < window; ++j) {
        res[i + width] += c2[j] * v[i + j];
      }
    }
    return res;
  }

// handle border cases first because we need different coefficients
#if defined(_OPENMP)
#pragma omp parallel for private(i, j) schedule(static)
#endif
  for (i = 0; i < width; ++i) {
    float_vect b1(window, 0.0);
    b1[i] = 1.0;

    const float_vect c1(sg_coeff(b1, deg));
    for (j = 0; j < window; ++j) {
      res[i] += c1[j] * v[j];
      res[endidx - i] += c1[j] * v[endidx - j];
    }
  }

  // now loop over rest of data. reusing the "symmetric" coefficients.
  float_vect b2(window, 0.0);
  b2[width] = 1.0;
  const float_vect c2(sg_coeff(b2, deg));

#if defined(_OPENMP)
#pragma omp parallel for private(i, j) schedule(static)
#endif
  for (i = 0; i <= (vsize - window); ++i) {
    for (j = 0; j < window; ++j) {
      res[i + width] += c2[j] * v[i + j];
    }
  }
  return res;
}

/*! least squares fit a polynome of degree 'deg' to data in 'b'.
 *  then calculate the first derivative and return it. */
static float_vect lsqr_fprime(const float_vect& b, const int deg) {
  const int rows(b.size());
  const int cols(deg + 1);
  float_mat A(rows, cols);
  float_vect res(rows);

  // generate input matrix for least squares fit
  int i, j;
  for (i = 0; i < rows; ++i) {
    for (j = 0; j < cols; ++j) {
      A[i][j] = pow(double(i), double(j));
    }
  }

  float_mat c(invert(transpose(A) * A) * (transpose(A) * transpose(b)));

  int bsize = static_cast<int>(b.size());
  for (i = 0; i < bsize; ++i) {
    res[i] = c[1][0];
    for (j = 1; j < deg; ++j) {
      res[i] += c[j + 1][0] * double(j + 1) * pow(double(i), double(j));
    }
  }
  return res;
}

/*! \brief savitzky golay smoothed numerical derivative.
 *
 * This method means fitting a polynome of degree 'deg' to a sliding window
 * of width 2w+1 throughout the data.
 *
 * In contrast to the sg_smooth function we do a brute force attempt by
 * always fitting the data to a polynome of degree 'deg' and using the
 * result. */
float_vect sg_derivative(const float_vect& v, const int width, const int deg,
                         const double h) {
  float_vect res(v.size(), 0.0);
  int vsize = static_cast<int>(v.size());
  if ((width < 1) || (deg < 1) || (vsize < (2 * width + 2))) {
    // sgs_error("sgsderiv: parameter error.\n");
    return res;
  }

  const int window = 2 * width + 1;

  // handle border cases first because we do not repeat the fit
  // lower part
  float_vect b(window, 0.0);
  int i, j;

  for (i = 0; i < window; ++i) {
    b[i] = v[i] / h;
  }
  const float_vect c(lsqr_fprime(b, deg));
  for (j = 0; j <= width; ++j) {
    res[j] = c[j];
  }
  // upper part. direction of fit is reversed
  for (i = 0; i < window; ++i) {
    b[i] = v[v.size() - 1 - i] / h;
  }
  const float_vect d(lsqr_fprime(b, deg));
  for (i = 0; i <= width; ++i) {
    res[v.size() - 1 - i] = -d[i];
  }

// now loop over rest of data. wasting a lot of least squares calcs
// since we only use the middle value.
#if defined(_OPENMP)
#pragma omp parallel for private(i, j) schedule(static)
#endif
  for (i = 1; i < (vsize - window); ++i) {
    for (j = 0; j < window; ++j) {
      b[j] = v[i + j] / h;
    }
    res[i + width] = lsqr_fprime(b, deg)[width];
  }
  return res;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// fill-column: 76
// indent-tabs-mode: nil
// End:
