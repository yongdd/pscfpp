/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "LuSolver.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_machine.h>
#include <cmath>

namespace Pscf
{
  
   LuSolver::LuSolver()
    : luPtr_(0),
      gMatInverse_(0),
      permPtr_(0),
      maxAbsA_(0.0),
      n_(0)
   {
      // Initialize gs_vector b_ 
      b_.size = 0;
      b_.stride = 1;
      b_.data = 0;
      b_.block = 0;
      b_.owner = 0;

      // Initialize gsl_vector x_ 
      x_.size = 0;
      x_.stride = 1;
      x_.data = 0;
      x_.block = 0;
      x_.owner = 0;
   }

   LuSolver::~LuSolver()
   {
      if (n_ > 0)  {
         if (luPtr_) gsl_matrix_free(luPtr_);
         if (gMatInverse_) gsl_matrix_free(gMatInverse_);
         if (permPtr_) gsl_permutation_free(permPtr_);
      }
   }

   /*
   * Allocate memory and set n.
   */
   void LuSolver::allocate(int n)
   {
      UTIL_CHECK(n > 0);
      UTIL_CHECK(n_ == 0);
      UTIL_CHECK(luPtr_ == 0);
      UTIL_CHECK(permPtr_ == 0);
      luPtr_ = gsl_matrix_alloc(n, n);
      gMatInverse_ = gsl_matrix_alloc(n, n);
      permPtr_ = gsl_permutation_alloc(n);
      b_.size = n;
      x_.size = n;
      n_ = n;
   }

   /*
   * Compute the LU decomposition for later use.
   */
   void LuSolver::computeLU(const Matrix<double>& A)
   {
      UTIL_CHECK(n_ > 0);
      UTIL_CHECK(A.capacity1() == n_);
      UTIL_CHECK(A.capacity2() == n_);

      int i, j;
      int k = 0;
      maxAbsA_ = 0.0;
      for (i = 0; i < n_;  ++i) {
         for (j = 0; j < n_; ++j) {
            luPtr_->data[k] = A(i,j);
            if (std::abs(A(i,j)) > maxAbsA_) {
               maxAbsA_ = std::abs(A(i,j));
            }
            ++k;
         }
      }

      // LU decomposition with partial pivoting, P A = L U, stored in
      // the format of gsl_linalg_LU_decomp. It is not called because,
      // since GSL 2.7, its recursive algorithm (used for n > 24) returns
      // early at an exactly zero pivot and then reads pivot indices
      // that were never set, which crashes for a singular matrix. The
      // unblocked algorithm below (that of GSL for n <= 24, with the
      // same BLAS calls) leaves a zero pivot in place and continues, so
      // that isSingular() can report it.
      gsl_permutation_init(permPtr_);
      signum_ = 1;
      for (j = 0; j < n_; ++j) {

         // Pivot: the element of largest magnitude on or below the
         // diagonal in column j
         gsl_vector_view column;
         column = gsl_matrix_subcolumn(luPtr_, j, j, n_ - j);
         size_t jMax = gsl_blas_idamax(&column.vector);
         if (jMax >= (size_t)(n_ - j)) jMax = 0;
         int jPivot = j + (int)jMax;
         if (jPivot != j) {
            gsl_matrix_swap_rows(luPtr_, j, jPivot);
            gsl_permutation_swap(permPtr_, j, jPivot);
            signum_ = -signum_;
         }

         // Column j of L and update of the remaining submatrix (both
         // are skipped for a zero pivot, whose column below is zero)
         double ajj = gsl_matrix_get(luPtr_, j, j);
         if (ajj != 0.0 && j < n_ - 1) {
            int m = n_ - j - 1;
            gsl_vector_view l = gsl_matrix_subcolumn(luPtr_, j, j+1, m);
            if (std::abs(ajj) >= GSL_DBL_MIN) {
               gsl_blas_dscal(1.0/ajj, &l.vector);
            } else {
               for (i = 0; i < m; ++i) {
                  *gsl_vector_ptr(&l.vector, i) /= ajj;
               }
            }
            gsl_vector_view u = gsl_matrix_subrow(luPtr_, j, j+1, m);
            gsl_matrix_view a22;
            a22 = gsl_matrix_submatrix(luPtr_, j+1, j+1, m, m);
            gsl_blas_dger(-1.0, &l.vector, &u.vector, &a22.matrix);
         }

      }
   }

   /*
   * Is the matrix of the last computeLU numerically singular?
   */
   bool LuSolver::isSingular(double epsilon) const
   {
      UTIL_CHECK(n_ > 0);
      if (maxAbsA_ == 0.0) return true;
      for (int i = 0; i < n_; ++i) {
         double pivot = std::abs(gsl_matrix_get(luPtr_, i, i));
         if (!(pivot > epsilon*maxAbsA_)) return true;
      }
      return false;
   }

   /*
   * Solve Ax = b.
   */
   void LuSolver::solve(Array<double>& b, Array<double>& x)
   {
      UTIL_CHECK(n_ > 0);
      UTIL_CHECK(b.capacity() == n_);
      UTIL_CHECK(x.capacity() == n_);

      // Associate gsl_vectors b_ and x_ with Arrays b and x
      b_.data = b.cArray();
      x_.data = x.cArray();

      // Solve system of equations
      gsl_linalg_LU_solve(luPtr_, permPtr_, &b_, &x_);

      // Destroy temporary associations
      b_.data = 0;
      x_.data = 0;
   }

   /*
   * Solve Ax = b by a basic solution, with a zero pivot allowed.
   */
   void LuSolver::solveBasic(Array<double> const & b, Array<double>& x) 
   const
   {
      UTIL_CHECK(n_ > 0);
      UTIL_CHECK(b.capacity() == n_);
      UTIL_CHECK(x.capacity() == n_);
      UTIL_CHECK(b.cArray() != x.cArray());

      int i, j;
      double sum, pivot;

      // Row permutation, x = P b
      for (i = 0; i < n_; ++i) {
         x[i] = b[gsl_permutation_get(permPtr_, i)];
      }

      // Forward substitution, L y = P b (L has a unit diagonal)
      for (i = 1; i < n_; ++i) {
         sum = x[i];
         for (j = 0; j < i; ++j) {
            sum -= gsl_matrix_get(luPtr_, i, j)*x[j];
         }
         x[i] = sum;
      }

      // Back substitution, U x = y, with x[i] = 0 for a zero pivot
      for (i = n_ - 1; i >= 0; --i) {
         pivot = gsl_matrix_get(luPtr_, i, i);
         if (pivot == 0.0) {
            x[i] = 0.0;
         } else {
            sum = x[i];
            for (j = i + 1; j < n_; ++j) {
               sum -= gsl_matrix_get(luPtr_, i, j)*x[j];
            }
            x[i] = sum/pivot;
         }
      }
   }

   /*  
   * Find inverse 
   */  
   void LuSolver::inverse(Matrix<double>& inv)
   {   
      UTIL_CHECK(n_ > 0); 

      gMatInverse_->data = inv.cArray();
      gsl_linalg_LU_invert(luPtr_, permPtr_, gMatInverse_);

   } 

}
