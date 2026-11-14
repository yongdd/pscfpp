#ifndef PSCF_CORE_TRIDIAGONAL_PIVOT_TEST_H
#define PSCF_CORE_TRIDIAGONAL_PIVOT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/math/TridiagonalSolver.h>
#include <util/global.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Pscf;
using namespace Util;

/*
* TridiagonalSolver with partial pivoting: invertible matrices with zero
* or small elimination pivots, consecutive row swaps, small dimensions,
* scaling, reuse of a factorization and singular matrices.
*/
class CoreTridiagonalPivotTest : public UnitTest
{
public:
   void setUp() {}
   void tearDown() {}

   void check(TridiagonalSolver& solver, std::vector<double> const& diagonal,
              std::vector<double> const& upper, std::vector<double> const& lower,
              bool symmetric, double scale=1.0)
   {
      int n = diagonal.size();
      DArray<double> d, u, l, rhs, x, product;
      d.allocate(n); u.allocate(std::max(1,n-1)); l.allocate(std::max(1,n-1));
      rhs.allocate(n); x.allocate(n); product.allocate(n);
      for (int i = 0; i < n; ++i) d[i] = scale*diagonal[i];
      for (int i = 0; i < n-1; ++i) {
         u[i] = scale*upper[i]; l[i] = scale*lower[i];
      }
      // Construct b from the original matrix and the known x_i=i+1.
      for (int i = 0; i < n; ++i) {
         rhs[i] = d[i]*(i+1);
         if (i > 0) rhs[i] += l[i-1]*i;
         if (i < n-1) rhs[i] += u[i]*(i+2);
      }
      if (symmetric) solver.computeLU(d,u);
      else solver.computeLU(d,u,l);
      solver.solve(rhs,x);
      for (int i = 0; i < n; ++i) {
         TEST_ASSERT(std::abs(x[i]-(i+1)) < 1e-10*(i+1));
         long double ax = static_cast<long double>(d[i])*x[i];
         if (i > 0) ax += static_cast<long double>(l[i-1])*x[i-1];
         if (i < n-1) ax += static_cast<long double>(u[i])*x[i+1];
         TEST_ASSERT(std::abs((ax-rhs[i])/scale) < 1e-10*(1+std::abs(rhs[i]/scale)));
         TEST_ASSERT(d[i] == scale*diagonal[i]);
         if (i < n-1) {
            TEST_ASSERT(u[i] == scale*upper[i]);
            TEST_ASSERT(l[i] == scale*lower[i]);
         }
      }
      for (int i = 0; i < n; ++i) x[i] = i+1;
      solver.multiply(x,product);
      solver.multiply(x,x);
      for (int i = 0; i < n; ++i) {
         TEST_ASSERT(std::abs((product[i]-rhs[i])/scale) < 1e-12);
         TEST_ASSERT(x[i] == product[i]);
      }
      // In-place solve and repeated RHS use the same factorization.
      solver.solve(rhs,rhs);
      for (int i = 0; i < n; ++i)
         TEST_ASSERT(std::abs(rhs[i]-(i+1)) < 1e-10*(i+1));
   }

   void testZeroAndSmallPivots()
   {
      printMethod(TEST_FUNC);
      TridiagonalSolver solver;
      solver.allocate(3);
      for (bool symmetric : {false,true}) {
         check(solver,{1,1,1},{1,1},{1,1},symmetric);
         check(solver,{4,4,4},{1,1},{1,1},symmetric);
         check(solver,{0,2,2},{1,1},{1,1},symmetric);
         check(solver,{1e-20,2,2},{1,1},{1,1},symmetric);
         check(solver,{1,3,4},{2,5},{2,5},symmetric);
      }
   }

   void testConsecutiveRowSwaps()
   {
      printMethod(TEST_FUNC);
      TridiagonalSolver solver;
      solver.allocate(6);
      check(solver,{.25,.25,.25,.25,.25,.25},
                   {.125,.125,.125,.125,.125},{2,2,2,2,2},false);
      check(solver,{4,4,4,4,4,4},{1,1,1,1,1},{.5,.5,.5,.5,.5},false);
      check(solver,{0,0,0,0,0,1},{2,2,2,2,2},{1,1,1,1,1},false);
   }

   void testSmallDimensionsAndScaling()
   {
      printMethod(TEST_FUNC);
      TridiagonalSolver scalar, pair, triple;
      scalar.allocate(1); pair.allocate(2); triple.allocate(3);
      check(scalar,{2},{},{},true);
      check(pair,{0,2},{1},{1},true);
      check(pair,{0,2},{2},{1},false);
      for (double scale : {1e-150, 1e150})
         check(triple,{1,1,1},{1,1},{1,1},false,scale);
   }

   void testSingularFailureAndReuse()
   {
      printMethod(TEST_FUNC);
      TridiagonalSolver solver;
      solver.allocate(2);
      DArray<double> d,u,b,x;
      d.allocate(2); u.allocate(1); b.allocate(2); x.allocate(2);
      d[0]=d[1]=u[0]=1; b[0]=b[1]=1;
      bool thrown=false;
      try { solver.computeLU(d,u); } catch (Exception&) { thrown=true; }
      TEST_ASSERT(thrown);
      thrown=false;
      try { solver.solve(b,x); } catch (Exception&) { thrown=true; }
      TEST_ASSERT(thrown);
      check(solver,{0,2},{1},{1},true);
   }
};

TEST_BEGIN(CoreTridiagonalPivotTest)
TEST_ADD(CoreTridiagonalPivotTest, testZeroAndSmallPivots)
TEST_ADD(CoreTridiagonalPivotTest, testConsecutiveRowSwaps)
TEST_ADD(CoreTridiagonalPivotTest, testSmallDimensionsAndScaling)
TEST_ADD(CoreTridiagonalPivotTest, testSingularFailureAndReuse)
TEST_END(CoreTridiagonalPivotTest)
#endif
