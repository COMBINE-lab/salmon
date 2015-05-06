/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Sailfish.

    Sailfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sailfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sailfish.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/


#ifndef MATRIX_TOOLS_HPP
#define MATRIX_TOOLS_HPP

#include <unordered_set>
#include <unordered_map>
#include <map>
#include <memory>

#include "nnls.h"

// taucs.h will define min/max macros if it's not already done (e.g. by Windows.h).
#ifndef min
    #define CGAL_TAUCS_DEFINES_MIN
#endif
#ifndef max
    #define CGAL_TAUCS_DEFINES_MAX
#endif

// TAUCS is a C library
extern "C" {
    #include "tsnnls.h"
}

// Undefine Taucs' min/max macros to avoid an error
// with std::min()/std::max() calls in standard C++ headers.
#ifdef CGAL_TAUCS_DEFINES_MIN
    #undef min
#endif
#ifdef CGAL_TAUCS_DEFINES_MAX
    #undef max
#endif


#include <vigra/regression.hxx>
#include <Eigen/SparseCore>
#include <boost/dynamic_bitset.hpp>
#include "shotgun_lasso.hpp"

namespace matrix_tools {

/**
*  This function collapses the rate matrix "A" and the corresponding vector of read
*  counts "counts" into a set of unique categories as defined in the paper
*  (Salzman, Jiang and Wong 2011).

*  Example
*  A = [0 1 0 0 0 1 1]
*      [0 0 1 1 0 1 1]
*      [0 0 0 1 1 0 0]
*  counts = [ 0 4 4 8 5 3 2 ]
*
*  might (b/c output column order is undefined) be collapsed to:
*  A = [1 0 0 0 2]
*      [0 1 1 0 2]
*      [0 0 1 1 0]
*  counts = [4 4 8 5 5]
*/
void collapseIntoCategories(
    std::vector<std::vector<double>> &A,
    std::vector<double> &counts,
    std::unique_ptr<std::vector<std::vector<double>>> &collapsedA,
    std::unique_ptr<std::vector<double>> &collapsedCounts
) {

    using boost::dynamic_bitset;
using Count = size_t;
using Index = size_t;

    size_t numRows = A.size();
    size_t numCols = A[0].size();
    std::map<dynamic_bitset<>, std::vector<Index> > categoryCounts;

    for ( size_t i = 0; i < numCols; ++i ) {
        dynamic_bitset<> category(numRows);

        for ( size_t j = 0; j < numRows; ++j ) {
            category[j] = ( A[j][i] > 0.0 ) ? 1 : 0;
        }
        // The trivial category contains all 0s; we simply omit it
        bool trivial = !category.any();
        if ( !trivial ) {
            categoryCounts[ category ].push_back(i);
        }
    }

    size_t numCategories = categoryCounts.size();
    std::unique_ptr<std::vector<std::vector<double>>> captr(new std::vector<std::vector<double>>(
                numRows, std::vector<double>( numCategories, 0.0 ) ));
    collapsedA = std::move(captr);

    std::unique_ptr<std::vector<double>> ccptr( new std::vector<double>(numCategories, 0.0) );
    collapsedCounts = std::move(ccptr);

    size_t j = 0;
    for ( auto catIt = categoryCounts.begin(); catIt != categoryCounts.end(); ++catIt, ++j ) {
        auto category = catIt->first;
        auto origCols = catIt->second;
        for ( auto origColID : origCols ) {
            (*collapsedCounts)[j] += counts[origColID];
        }
        for ( size_t i = 0; i < numRows; ++i ) {
            auto accumCount = 0.0;
            for ( auto col : origCols ) { accumCount += A[i][col]; }
            (*collapsedA)[i][j] = accumCount * category[i];
        }
    }

}



template< typename MatT, typename VecT >
std::vector<double> LARSSolve( MatT &Ain, VecT& bin )
{

        int m = Ain.rows(); 
        int n = Ain.cols();
        vigra::linalg::Matrix<double> A(m, n);
        vigra::linalg::Matrix<double> b(m, 1);
        
        // fill A and b
        for (int k=0; k < Ain.outerSize(); ++k) {

          for (typename MatT::InnerIterator it(Ain,k); it; ++it) {
            auto val = it.value();
            auto i = it.row();   // row index
            auto j = it.col();   // col index (here it is equal to k)
            A(i,j) = val;
         }
        }

        for (auto i : boost::irange(size_t(0), bin.size())) {
            b(i,0) = bin[i];
        }

        // normalize the input
        vigra::linalg::Matrix<double> offset(1,n), scaling(1,n);
        vigra::linalg::prepareColumns(A, A, offset, scaling,  vigra::linalg::DataPreparationGoals(vigra::linalg::ZeroMean|vigra::linalg::UnitVariance));
        vigra::linalg::prepareColumns(b, b, vigra::linalg::DataPreparationGoals(vigra::linalg::ZeroMean));

        // arrays to hold the output
        vigra::ArrayVector<vigra::ArrayVector<int> > activeSets;
        vigra::ArrayVector<vigra::linalg::Matrix<double> > solutions;

        // run leastAngleRegression() in non-negative LASSO mode
        int numSolutions = vigra::linalg::leastAngleRegression(A, b, activeSets, solutions,
                                    vigra::linalg::LeastAngleRegressionOptions().nnlasso());

        // print results
        std::vector<double> x(n);
        vigra::linalg::Matrix<double> denseSolution(1, n);
        for (vigra::MultiArrayIndex k = 0; k < numSolutions; ++k) {
            // transform the sparse solution into a dense vector
            denseSolution.init(0.0); // ensure that inactive variables are zero
            for (unsigned int i = 0; i < activeSets[k].size(); ++i) {
                // set the values of the active variables;
                // activeSets[k][i] is the true index of the i-th variable in the active set
                denseSolution(0, activeSets[k][i]) = solutions[k](i,0);
            }

            // invert the input normalization
            denseSolution = denseSolution * pointWise(scaling);
            for (auto i : boost::irange(int(0), n)) {
              x[i] = denseSolution(0,i);
            }
        }

        return x;
}

template< typename MatT, typename VecT >
std::vector<double> shotgunSolve( MatT &A, VecT& b ) {

    // Try shotgun LASSO solver
    shotgun_data prob;
    //std::cerr << "trying to reserve " << A.nonZeros() << " nonzeros in prob rows/cols\n";
    prob.A_cols.resize( A.nonZeros(), sparse_array() );
    prob.A_rows.resize( A.nonZeros(), sparse_array() );

    std::cerr << "done\n";

    auto nnz = A.nonZeros();
    size_t ctr{0};
    for (int k=0; k < A.outerSize(); ++k) {
      //std::cerr << "column " << k << " of " << A.outerSize() << "\n";
      for (typename MatT::InnerIterator it(A,k); it; ++it) {
        auto val = it.value();
        auto i = it.row();   // row index
        auto j = it.col();   // col index (here it is equal to k)
        //std::cerr << "adding " << i << ", " << j << " : " << val << "\n";
        try {
        prob.A_cols[j].add(i, val);
        prob.A_rows[i].add(j, val);
        ++ctr;
        if ( ctr > nnz ) { std::cerr << "only reserved space for " << nnz <<
          " elements, but you're pushing on the " << ctr << "th element\n";
        }
        } catch ( std::exception& e ) {
          std::cerr << "died trying to add " << i << ", " << j << " : " << val << "\n";
          std::cerr << e.what();
        }
      }
    }

    prob.nx = A.cols();
    prob.ny = A.rows();

    double me = *std::max_element(b.begin(), b.end());
    auto scale = 1.0 / me;

    prob.y.reserve(b.size());
    for ( auto e : b ) { prob.y.push_back(e * scale); };
        
    double lambda = 0.1;
    int K = 200;
    int maxiter = 5000000;
    int verbose = 3;
    double threshold = 1e-15;
    assert( prob.ny == prob.y.size() );
    std::cerr << "before solve lasso\n";
    LassoProblem lprob(&prob, lambda, K, threshold, maxiter, verbose);
    lprob.solve();
    std::cerr << "after solve lasso\n";

    std::vector<double> myX(prob.x);
    for ( auto& e : myX ) { e *= me; }
    return myX;
}

/**
  * Given the matrix A (m x n) and the right-hand-side vector b (m x 1) solve
  * for the solution vector x ( n x 1 ) satisfying
  * min_{x} || Ax - b ||
  * subject to x[i] >= 0 for all i
  */
template< typename MatT, typename VecT >
std::vector<double> nnlsSolve( MatT &A, VecT &b ) {
    //std::cerr << "in least squares problem\n";
    
    Eigen::MatrixXd ADense = Eigen::MatrixXd(A);//.transpose();

    int    mda = ADense.rows();
    int    m = ADense.rows(), n = ADense.cols();

    std::vector<double> w(n, 0.0);
    std::vector<double> zz(m, 0.0);
    std::vector<double> x(n, 1.0 );
    std::vector<int> indx(n, 0);
    assert(b.size() == m);
    double rnorm = 0.0;
    int mode = 0;

    nnls( ADense.data(), mda, m, n, &b[0], &x[0], &rnorm, &w[0], &zz[0], &indx[0], &mode);
        
    return x;
    /*
    size_t errCtr = 0;
    // Create a new TAUCS matrix structure
    taucs_ccs_matrix* mat;
    mat = new taucs_ccs_matrix;

    mat->n = A.cols();
    mat->m = A.rows();
    mat->flags = TAUCS_DOUBLE;

    std::cerr << "inited TAUCS matrix\n";

    // Compress this matrix so we can steal / share the data
    A.makeCompressed();

    std::cerr << "compressed matrix\n";

    mat->colptr = A.outerIndexPtr();
    mat->rowind = A.innerIndexPtr();
    mat->values.d = A.valuePtr();

    double residualNorm{0.0};

    std::cerr << "called NNLS solver\n";
    auto x = t_snnls_fallback( mat, &b[0], &residualNorm, 0.0, 0);
    char *errString;
    tsnnls_error(&errString);
    std::cerr << "NNLS solver returned\n";
    std::cerr << "ERROR STRING " << errString << "\n";

    if( x == NULL ) {
      std::cerr << "AHHH, x is STILL NULL\n";

      std::stringstream ss;
      ss << "errorA_" << errCtr << ".mtx";
      std::ofstream Afile(ss.str());

      Afile << "SPARSE\n";
      Afile << A.rows() << '\t' << A.cols() << '\n' << A.nonZeros() << '\n';
      for (int k=0; k < A.outerSize(); ++k) {
        for (typename MatT::InnerIterator it(A,k); it; ++it) {
          Afile << it.row()+1 << '\t' << it.col()+1 << '\t' << it.value() << '\n';
        }
      }
      Afile.close();

      ss.clear();
      ss << "errorb_"<< errCtr << ".mat";
      std::ofstream bfile("errorb.mtx");
      bfile << "# error right hand side\n";
      bfile << "# rows: " << b.size() << "\n";
      bfile << "# columns: 1\n";
      for (size_t i = 0; i < b.size(); ++i) {
          bfile << b[i] << '\n';
      }
      bfile.close();
      ++errCtr;

      delete mat;
      return std::vector<double>();
      //std::abort();
    } else {

      std::vector<double> myX( A.cols(), 0.0 );

      for( size_t i = 0; i < myX.size(); ++i ) { myX[i] = x[i]; std::cerr << "x[" << i << "]  = " << x[i] << "\n"; }

      free(x);
      delete mat;

      return myX;
    }
    */

    /*
    // Try shotgun LASSO solver
    shotgun_data prob;
    std::cerr << "trying to reserve " << A.nonZeros() << " nonzeros in prob rows/cols\n";
    prob.A_cols.resize( A.nonZeros(), sparse_array() );
    prob.A_rows.resize( A.nonZeros(), sparse_array() );

    std::cerr << "done\n";

    auto nnz = A.nonZeros();
    size_t ctr{0};
    for (int k=0; k < A.outerSize(); ++k) {
      std::cerr << "column " << k << " of " << A.outerSize() << "\n";
      for (typename MatT::InnerIterator it(A,k); it; ++it) {
        auto val = it.value();
        auto i = it.row();   // row index
        auto j = it.col();   // col index (here it is equal to k)
        //std::cerr << "adding " << i << ", " << j << " : " << val << "\n";
        try {
        prob.A_cols[j].add(i, val);
        prob.A_rows[i].add(j, val);
        ++ctr;
        if ( ctr > nnz ) { std::cerr << "only reserved space for " << nnz <<
          " elements, but you're pushing on the " << ctr << "th element\n";
        }
        } catch ( std::exception& e ) {
          std::cerr << "died trying to add " << i << ", " << j << " : " << val << "\n";
          std::cerr << e.what();
        }
      }
    }

    prob.nx = A.cols();
    prob.ny = A.rows();

    prob.y.reserve(b.size());
    for ( auto e : b ) { prob.y.push_back(e); };
    double lambda = 0.0;
    int K = 10;
    int maxiter = 5000;
    int verbose = 10;
    double threshold = 1e-8;
    assert( prob.ny == prob.y.size() );
    std::cerr << "before solve lasso\n";
    solveLasso(&prob, lambda, K, threshold, maxiter, verbose);
    std::cerr << "after solve lasso\n";

    std::vector<double> myX(prob.x);
    return myX;
    */

}





template <typename MatT>
std::vector<size_t> markDuplicateColumns( MatT &M ) {

    std::vector<size_t> removeCols;
    for ( size_t i = 0; i < M.cols(); ++i ) {
        Eigen::SparseVector<typename MatT::Scalar> colI(M.middleCols(i, 1));
        auto normI = colI.dot(colI);
        auto stillValid = true;

        for ( size_t j = i + 1; j < M.cols(); ++j ) {
            if ( stillValid ) {
                Eigen::SparseVector<typename MatT::Scalar> colJ(M.middleCols(j, 1));
                auto normIJ = colI.dot(colJ);
                if ( normIJ > normI ) {
                    removeCols.push_back(i);
                    stillValid = false;
                }
            }
        } // end j
    } // end i

    return removeCols;
}

/**
* Project vector v onto vector u
*/
template <typename VecT>
VecT projectOnto( VecT &v, VecT &u ) {
    return (v.dot(u) / u.dot(u) ) * u;
}

template <typename ColT>
bool classicalGramSchmidt( std::vector<ColT> &basis, typename ColT::Scalar tol = 1e-10 ) {

    // Orthogonalize the basis
    for ( size_t j = 0; j < basis.size(); ++j ) {
        // vj = aj
        ColT vj(basis[j]);
        for ( size_t i = 0; i < j; ++i ) {
            // rij = qi * aj
            //auto rij = basis[i].dot(basis[j]);
            // vj = vj - rii * qi
            vj = vj - projectOnto( basis[j], basis[i] );
            //(basis[i] * rij);
        }
        // rjj = ||vj||_{2}
        auto rjj = vj.dot(vj);
        // If this vector isn't orthogonal, then return false now
        if ( rjj <  tol ) {
            return false;
        }
        // qj = vj / rjj
        basis[j] = vj * (1.0 /  rjj);
    }

    return true;
}

template <typename ColT>
bool addColumnToBasis( std::vector<ColT> &basis, ColT &c, typename ColT::Scalar tol = 1e-10 ) {

    ColT vj(c);//basis[j]);
    for ( size_t i = 0; i < basis.size(); ++i ) {
        // rij = qi * aj
        //auto rij = basis[i].dot(c);
        // vj = vj - rii * qi
        vj = vj - projectOnto( c, basis[i] );
        //vj = vj - (basis[i] * rij);
    }
    // rjj = ||vj||_{2}
    auto rjj = vj.dot(vj);

    if ( rjj >= tol ) {
        // qj = vj / rjj
        basis.push_back( vj / rjj);
        return true;
    } else {
        return false;
    }
}

template <typename MatT>
std::vector<size_t> markDependentColumns( MatT &M ) {

using MatT::Index = typename;
using MatT::Scalar = typename;

    // Reference from MathOverflow:
    // http://mathoverflow.net/questions/109868/finding-linearly-independent-columns-of-a-large-sparse-rectangular-matrix

    // Structure that holds the column index and it's height
    struct ColHeightT {
        IndexT col;
        IndexT height;
    };

    // Compute the "height" of each column, where the "height" of a column is the
    // index of the row where the first non-zero entry occurs
    std::vector< ColHeightT > heights;
    heights.reserve(M.cols());

    for ( size_t i = 0; i < M.cols(); ++i ) {
        typename MatT::InnerIterator it(M, i);
        heights.push_back( {i, it.row()} );
    }

    // Sort the columns by height
    std::sort( heights.begin(), heights.end(),
               []( const ColHeightT & a, const ColHeightT & b) -> bool { return a.height < b.height; } );

    // We'll always keep the first column
    std::unordered_set<size_t> independentCols { heights.front().col };
    IndexT currHeight { heights.front().height };

    // Pick at least one column from each height class
    for ( auto & ch : heights ) {
        if ( ch.height > currHeight ) {
            independentCols.insert(ch.col);
            currHeight = ch.height;
        }
    }

    std::vector< Eigen::SparseVector< ScalarT > > basis;
    for ( auto c : independentCols ) {
        Eigen::SparseVector<ScalarT> col(M.middleCols(c, 1));
        assert(col.size() == M.rows());
        basis.push_back( col );
    }

    // orthogonalize the current basis
    auto isIndependent = classicalGramSchmidt(basis);
    if (!isIndependent) {
        std::cerr << "Something went horribly wrong; the original set of vectors was not linearly independent!\n";
        std::cerr << "Doing the safe thing (start with just one col and build up\n";
        independentCols = {0};
        Eigen::SparseVector<ScalarT> col(M.middleCols(0, 1));
        basis = { col };
        classicalGramSchmidt(basis);
        //std::abort();
    }

    std::vector<size_t> removeCols;
    // For each column that's not already in our basis
    for ( size_t i = 0; i < M.cols(); ++i ) {
        if ( independentCols.find(i) == independentCols.end() ) {
            // Try to add it to our basis
            Eigen::SparseVector<ScalarT> c( M.middleCols(i, 1) );
            assert(c.size() == M.rows());
            auto independent = addColumnToBasis(basis, c);
            // If we were able to add it, it's independent
            if ( independent ) {
                independentCols.insert(i);
            } else { // Otherwise it's dependent
                removeCols.push_back(i);
            }
        } // don't mess with already independent columsn
    } // done loop over columns

    return removeCols;
}
}


#endif // MATRIX_TOOLS_HPP
