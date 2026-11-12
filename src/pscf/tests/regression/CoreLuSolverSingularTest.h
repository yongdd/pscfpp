#ifndef PSCF_CORE_LU_SOLVER_SINGULAR_TEST_H
#define PSCF_CORE_LU_SOLVER_SINGULAR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/math/LuSolver.h>
#include <util/containers/DMatrix.h>
#include <util/containers/DArray.h>

#include <cmath>

using namespace Util;
using namespace Pscf;

class CoreLuSolverSingularTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   /*
   * Set a diagonally dominant n x n matrix with fixed pseudo-random
   * off-diagonal elements in (-1, 1).
   */
   void setMatrix(DMatrix<double>& a, int n)
   {
      for (int i = 0; i < n; ++i) {
         for (int j = 0; j < n; ++j) {
            a(i,j) = std::sin(1.0 + 7.0*i + 3.0*j*j);
         }
         a(i,i) += 2.0*std::sqrt(double(n));
      }
   }

   /*
   * Return max_i |(A x - b)_i|.
   */
   double residual(DMatrix<double> const & a, DArray<double> const & x,
                   DArray<double> const & b, int n)
   {
      double r, rMax = 0.0;
      for (int i = 0; i < n; ++i) {
         r = -b[i];
         for (int j = 0; j < n; ++j) {
            r += a(i,j)*x[j];
         }
         if (!(std::abs(r) <= rMax)) rMax = std::abs(r);
      }
      return rMax;
   }

   /*
   * A 100 x 100 matrix (GSL decomposes matrices larger than 24 x 24
   * recursively) whose columns 0 and 50 are e_50 and -e_50, as the
   * columns of the r1 NR jacobian for the w fields at a grid point
   * that changes no concentration. gsl_linalg_LU_decomp (GSL 2.7 and
   * later) stopped at the zero pivot and read pivot indices it had
   * not set, which crashed. computeLU must complete, isSingular must
   * report the matrix, and solveBasic must solve a consistent system.
   */
   void testSingular()
   {
      printMethod(TEST_FUNC);

      const int n = 100;
      DMatrix<double> a;
      a.allocate(n, n);
      setMatrix(a, n);
      for (int i = 0; i < n; ++i) {
         a(i, 0) = 0.0;
         a(i, n/2) = 0.0;
      }
      a(n/2, 0) = 1.0;
      a(n/2, n/2) = -1.0;

      LuSolver solver;
      solver.allocate(n);
      solver.computeLU(a);
      TEST_ASSERT(solver.isSingular());
      TEST_ASSERT(solver.isSingular(0.0));

      // Consistent right-hand side b = A x0
      DArray<double> x0, b, x;
      x0.allocate(n);
      b.allocate(n);
      x.allocate(n);
      for (int i = 0; i < n; ++i) {
         x0[i] = std::cos(0.5 + 1.3*i);
      }
      for (int i = 0; i < n; ++i) {
         b[i] = 0.0;
         for (int j = 0; j < n; ++j) {
            b[i] += a(i,j)*x0[j];
         }
      }
      solver.solveBasic(b, x);
      TEST_ASSERT(residual(a, x, b, n) < 1.0E-10);

      // The solution is unique up to the null vector e_0 + e_50
      for (int i = 0; i < n; ++i) {
         if (i == 0 || i == n/2) continue;
         TEST_ASSERT(std::abs(x[i] - x0[i]) < 1.0E-10);
      }
      TEST_ASSERT(std::abs((x[0] - x[n/2]) - (x0[0] - x0[n/2])) < 1.0E-10);
   }

   /*
   * For a nonsingular 100 x 100 matrix, solve and solveBasic give the
   * solution.
   */
   void testNonsingular()
   {
      printMethod(TEST_FUNC);

      const int n = 100;
      DMatrix<double> a;
      a.allocate(n, n);
      setMatrix(a, n);

      LuSolver solver;
      solver.allocate(n);
      solver.computeLU(a);
      TEST_ASSERT(!solver.isSingular());

      DArray<double> x0, b, x, y;
      x0.allocate(n);
      b.allocate(n);
      x.allocate(n);
      y.allocate(n);
      for (int i = 0; i < n; ++i) {
         x0[i] = std::cos(0.5 + 1.3*i);
      }
      for (int i = 0; i < n; ++i) {
         b[i] = 0.0;
         for (int j = 0; j < n; ++j) {
            b[i] += a(i,j)*x0[j];
         }
      }
      solver.solve(b, x);
      solver.solveBasic(b, y);
      for (int i = 0; i < n; ++i) {
         TEST_ASSERT(std::abs(x[i] - x0[i]) < 1.0E-10);
         TEST_ASSERT(std::abs(y[i] - x0[i]) < 1.0E-10);
      }
   }

};

TEST_BEGIN(CoreLuSolverSingularTest)
TEST_ADD(CoreLuSolverSingularTest, testSingular)
TEST_ADD(CoreLuSolverSingularTest, testNonsingular)
TEST_END(CoreLuSolverSingularTest)

#endif
