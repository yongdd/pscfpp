#ifndef PSCF_GPU_COMPLEX_EXP_ALIAS_TEST_H
#define PSCF_GPU_COMPLEX_EXP_ALIAS_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/HostSendArray.h>
#include <complex>
#include <cmath>

using namespace Util;
using namespace Pscf;

/* Exact aliasing must agree with the analytic complex exponential. */
class GpuComplexExpAliasTest : public UnitTest
{
public:

   void check(bool slice)
   {
      const int n = 6;
      double real[n] = {-.2, 1.0, 0.0, -.3, 2.0, .4};
      double imag[n] = {.8, .5, 1.0, -.7, 0.0, -1.2};
      HostSendArray<cudaComplex,CUT> original(n), actual(n), distinct(n);
      for (int i = 0; i < n; ++i) {
         original[i].x = real[i];
         original[i].y = imag[i];
      }
      BackendArray<cudaComplex,CUT> input(n), output(n), inplace(n);
      input = original;
      output = original;
      inplace = original;
      int begin = slice ? 1 : 0;
      int count = slice ? n-2 : n;
      if (slice) {
         VecOp::expV(output, input, begin, begin, count);
         VecOp::expV(inplace, inplace, begin, begin, count);
      } else {
         VecOp::expV(output, input);
         VecOp::expV(inplace, inplace);
      }
      actual = inplace;
      distinct = output;
      HostSendArray<cudaComplex,CUT> preserved(n);
      preserved = input;
      #ifdef SINGLE_PRECISION
      const double tolerance = 1e-5;
      #else
      const double tolerance = 1e-12;
      #endif
      for (int i = 0; i < n; ++i) {
         std::complex<double> initial(original[i].x, original[i].y);
         bool changed = (i >= begin && i < begin+count);
         std::complex<double> expected = changed ? std::exp(initial) : initial;
         std::complex<double> a(actual[i].x, actual[i].y);
         std::complex<double> b(distinct[i].x, distinct[i].y);
         TEST_ASSERT(std::abs(a-expected) < tolerance*(1+std::abs(expected)));
         TEST_ASSERT(std::abs(b-expected) < tolerance*(1+std::abs(expected)));
         TEST_ASSERT(std::abs(a-b) < tolerance);
         TEST_ASSERT(preserved[i].x == original[i].x);
         TEST_ASSERT(preserved[i].y == original[i].y);
         if (!changed) {
            TEST_ASSERT(actual[i].x == original[i].x);
            TEST_ASSERT(actual[i].y == original[i].y);
         }
      }
   }

   void testWholeArray()
   {
      printMethod(TEST_FUNC);
      check(false);
   }

   void testExactAliasSlice()
   {
      printMethod(TEST_FUNC);
      check(true);
   }
};

TEST_BEGIN(GpuComplexExpAliasTest)
TEST_ADD(GpuComplexExpAliasTest, testWholeArray)
TEST_ADD(GpuComplexExpAliasTest, testExactAliasSlice)
TEST_END(GpuComplexExpAliasTest)

#endif
