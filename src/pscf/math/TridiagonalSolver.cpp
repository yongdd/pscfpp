/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <string>
#include <iostream>

#include "TridiagonalSolver.h"
#include <util/global.h>
#include <algorithm>
#include <cmath>

namespace Pscf
{
  
   /*
   * Constructor.
   */
   TridiagonalSolver::TridiagonalSolver()
    : n_(0), isDecomposed_(false)
   {}

   /*
   * Destructor.
   */
   TridiagonalSolver::~TridiagonalSolver()
   {}

   /*
   * Allocate memory.
   */
   void TridiagonalSolver::allocate(int n)
   {
      UTIL_CHECK(n > 0);
      d_.allocate(n);
      u_.allocate(std::max(1, n-1));
      u2_.allocate(std::max(1, n-2));
      l_.allocate(std::max(1, n-1));
      pivot_.allocate(std::max(1, n-1));
      originalD_.allocate(n);
      originalU_.allocate(std::max(1, n-1));
      originalL_.allocate(std::max(1, n-1));
      y_.allocate(n);
      n_ = n;
   }

   /*
   * Compute the LU decomposition of a symmetric matrix for later use.
   */
   void TridiagonalSolver::computeLU(const DArray<double>& d, 
                                     const DArray<double>& u)
   {  computeLU(d, u, u); }

   /*
   * Compute the LU decomposition of a symmetric matrix for later use.
   */
   void TridiagonalSolver::computeLU(const DArray<double>& d, 
                                     const DArray<double>& u,
                                     const DArray<double>& l)
   {
      isDecomposed_ = false;
      UTIL_CHECK(n_ > 0);
      UTIL_CHECK(d.capacity() >= n_);
      UTIL_CHECK(u.capacity() >= n_-1);
      UTIL_CHECK(l.capacity() >= n_-1);

      // Copy to local arrays
      for (int i = 0; i < n_; ++i) {
         UTIL_CHECK(std::isfinite(d[i]));
         originalD_[i] = d_[i] = d[i];
      }
      for (int i = 0; i < n_ - 1; ++i) {
         UTIL_CHECK(std::isfinite(u[i]));
         UTIL_CHECK(std::isfinite(l[i]));
         originalU_[i] = u_[i] = u[i];
         originalL_[i] = l_[i] = l[i];
         pivot_[i] = i;
      }
      for (int i = 0; i < n_-2; ++i) u2_[i] = 0.0;
      gaussElimination();
      isDecomposed_ = true;
   }

   void TridiagonalSolver::gaussElimination()
   {
      // Gauss elimination
      for (int i = 0; i < n_ - 1; ++i) {
         if (std::abs(d_[i]) >= std::abs(l_[i])) {
            UTIL_CHECK(d_[i] != 0.0);
            double q = l_[i]/d_[i];
            l_[i] = q;
            d_[i+1] -= q*u_[i];
         } else {
            // Swap rows i and i+1 before elimination. Only a second
            // upper diagonal is introduced, so the factors remain banded.
            double q = d_[i]/l_[i];
            d_[i] = l_[i];
            l_[i] = q;
            double upper = u_[i];
            u_[i] = d_[i+1];
            d_[i+1] = upper - q*d_[i+1];
            if (i < n_-2) {
               u2_[i] = u_[i+1];
               u_[i+1] *= -q;
            }
            pivot_[i] = i+1;
         }
         UTIL_CHECK(std::isfinite(d_[i+1]));
      }
      UTIL_CHECK(d_[n_-1] != 0.0);

      #if 0
      std::cout << "\n";
      for (int i=0; i < n_; ++i) {
         std::cout << "  " << d_[i];
      }
      std::cout << "\n";
      for (int i=0; i < n_ - 1; ++i) {
         std::cout << "  " << u_[i];
      }
      std::cout << "\n";
      for (int i=0; i < n_ - 1; ++i) {
         std::cout << "  " << l_[i];
      }
      std::cout << "\n";
      #endif
   }

   /*
   * Compute matrix product Ab = x for x, given known vector b.
   */
   void TridiagonalSolver::multiply(const DArray<double>& b, DArray<double>& x)
   {
      UTIL_CHECK(isDecomposed_);
      UTIL_CHECK(b.capacity() >= n_);
      UTIL_CHECK(x.capacity() >= n_);
      // Stage the product so the input and output may be the same array.
      for (int i = 0; i < n_; ++i) {
         y_[i] = originalD_[i]*b[i];
         if (i > 0) y_[i] += originalL_[i-1]*b[i-1];
         if (i < n_-1) y_[i] += originalU_[i]*b[i+1];
      }
      for (int i = 0; i < n_; ++i) x[i] = y_[i];
   }

   /*
   * Solve Ax = b for x, given known b.
   */
   void TridiagonalSolver::solve(const DArray<double>& b, DArray<double>& x)
   {
      UTIL_CHECK(isDecomposed_);
      UTIL_CHECK(b.capacity() >= n_);
      UTIL_CHECK(x.capacity() >= n_);
      for (int i = 0; i < n_; ++i) y_[i] = b[i];

      // Apply the row permutations and lower-triangular eliminations.
      for (int i = 0; i < n_-1; ++i) {
         if (pivot_[i] == i) {
            y_[i+1] -= l_[i]*y_[i];
         } else {
            double old = y_[i];
            y_[i] = y_[i+1];
            y_[i+1] = old - l_[i]*y_[i];
         }
      }

      // U has two upper diagonals after row interchanges.
      x[n_-1] = y_[n_-1]/d_[n_-1];
      if (n_ > 1) x[n_-2] = (y_[n_-2] - u_[n_-2]*x[n_-1])/d_[n_-2];
      for (int i = n_-3; i >= 0; --i) {
         x[i] = (y_[i] - u_[i]*x[i+1] - u2_[i]*x[i+2])/d_[i];
      }
   }

}
