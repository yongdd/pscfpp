#ifndef PSCF_CORE_LU_SOLVER_TEST_H
#define PSCF_CORE_LU_SOLVER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/math/LuSolver.h>
#include <util/containers/DMatrix.h>

using namespace Util;
using namespace Pscf;

class CoreLuSolverTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   /*
   * isSingular() reports a singular matrix after computeLU, so that a
   * caller can avoid inverse() and solve(), which abort through the
   * GSL error handler for a singular matrix. The test is relative to
   * the largest element, so it does not depend on the scale of A.
   */
   void testIsSingular()
   {
      printMethod(TEST_FUNC);

      DMatrix<double> a;
      a.allocate(3, 3);
      LuSolver solver;
      solver.allocate(3);
      int i, j;

      // Invertible matrix (as in LuSolverTest), at three scales
      double scales[3] = {1.0, 1.0E-20, 1.0E+20};
      for (int k = 0; k < 3; ++k) {
         a(0,0) = 1.0;
         a(1,1) = 3.0;
         a(2,2) = 4.0;
         a(0,1) = 2.0;
         a(1,0) = 2.0;
         a(1,2) = 5.0;
         a(2,1) = 5.0;
         a(0,2) = 0.0;
         a(2,0) = 0.0;
         for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
               a(i,j) *= scales[k];
            }
         }
         solver.computeLU(a);
         TEST_ASSERT(!solver.isSingular());
      }

      // Two equal rows (a chi matrix with one zero chi parameter)
      for (i = 0; i < 3; ++i) {
         for (j = 0; j < 3; ++j) {
            a(i,j) = 0.0;
         }
      }
      a(0,2) = 10.0;
      a(2,0) = 10.0;
      a(1,2) = 10.0;
      a(2,1) = 10.0;
      solver.computeLU(a);
      TEST_ASSERT(solver.isSingular());

      // Rank 1 matrix with no zero element
      for (i = 0; i < 3; ++i) {
         for (j = 0; j < 3; ++j) {
            a(i,j) = double((i + 1)*(j + 2));
         }
      }
      solver.computeLU(a);
      TEST_ASSERT(solver.isSingular());

      // Zero matrix
      for (i = 0; i < 3; ++i) {
         for (j = 0; j < 3; ++j) {
            a(i,j) = 0.0;
         }
      }
      solver.computeLU(a);
      TEST_ASSERT(solver.isSingular());
   }

};

TEST_BEGIN(CoreLuSolverTest)
TEST_ADD(CoreLuSolverTest, testIsSingular)
TEST_END(CoreLuSolverTest)

#endif
