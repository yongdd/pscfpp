#ifndef PSCF_GPU_VEC_RANDOM_TEST_H
#define PSCF_GPU_VEC_RANDOM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/backend/cuda/CUT.h>
#include <pscf/backend/cuda/VecRandom.h>
#include <pscf/backend/cuda/HostSendArray.h>
#include <util/misc/Exception.h>

#include <iostream>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests of the GPU random numbers (VecRandom<CUT>).
*/
class GpuVecRandomTest : public UnitTest 
{

public:

   void setUp()
   {  setVerbose(0); }

   void tearDown()
   {}

   void testNormalOddFloat()
   {
      printMethod(TEST_FUNC);
      checkNormalOdd<float>();
   }

   void testNormalOddDouble()
   {
      printMethod(TEST_FUNC);
      checkNormalOdd<double>();
   }

   /*
   * Normal random numbers for arrays with an odd number of elements.
   *
   * cuRAND generates normal values in pairs. Every element, including
   * the last one, must be set, and the values must have the requested
   * mean and standard deviation.
   */
   template <typename T>
   void checkNormalOdd()
   {
      VecRandom<CUT> random;
      random.setSeed(6712983651284);

      const T mean = 1.0;
      const T stddev = 0.5;
      const T unset = 1.0E+30; // marks elements that were not set

      // Odd sizes, from a single element (no pair at all) up
      int sizes[3] = {1, 3, 100001};
      for (int s = 0; s < 3; ++s) {
         int n = sizes[s];
         HostSendArray<T,CUT> h_(n);
         for (int i = 0; i < n; ++i) {
            h_[i] = unset;
         }
         BackendArray<T,CUT> d_(n);
         d_ = h_;

         bool thrown = false;
         try {
            random.normal(d_, stddev, mean);
         } catch (Exception& e) {
            thrown = true;
         }
         TEST_ASSERT(!thrown);

         h_ = d_;
         double ave = 0.0;
         double var = 0.0;
         double val = 0.0;
         for (int i = 0; i < n; ++i) {
            val = h_[i] - mean;
            // Deviation beyond 10 stddev has probability below 1E-22
            TEST_ASSERT(fabs(val) < 10.0*stddev);
            ave += val;
            var += val*val;
         }
         if (n > 1000) {
            ave = ave/double(n);
            var = var/double(n);
            var = sqrt(var) - stddev;
            ave = ave/stddev;
            var = var/stddev;
            if (verbose() > 0) {
               std::cout << std::endl;
               std::cout << "n = " << n << std::endl;
               std::cout << "Average  " << ave << std::endl;
               std::cout << "StdDev   " << var << std::endl;
            }
            TEST_ASSERT(fabs(ave) < 0.1);
            TEST_ASSERT(fabs(var) < 0.1);
         }
      }

      // The last element of an odd array is generated separately.
      // Check its distribution over many calls with a 3-element array.
      int m = 4000;
      BackendArray<T,CUT> d_(3);
      HostSendArray<T,CUT> h_(3);
      double ave = 0.0;
      double var = 0.0;
      double val = 0.0;
      for (int j = 0; j < m; ++j) {
         random.normal(d_, stddev, mean);
         h_ = d_;
         val = h_[2] - mean;
         ave += val;
         var += val*val;
      }
      ave = ave/double(m);
      var = var/double(m);
      var = sqrt(var) - stddev;
      ave = ave/stddev;
      var = var/stddev;
      if (verbose() > 0) {
         std::cout << "Last element of 3, " << m << " calls" << std::endl;
         std::cout << "Average  " << ave << std::endl;
         std::cout << "StdDev   " << var << std::endl;
      }
      // Standard errors are 1/sqrt(m) = 0.016 (ave) and
      // 1/sqrt(2m) = 0.011 (StdDev), in units of stddev
      TEST_ASSERT(fabs(ave) < 0.1);
      TEST_ASSERT(fabs(var) < 0.1);
   }

};

TEST_BEGIN(GpuVecRandomTest)
TEST_ADD(GpuVecRandomTest, testNormalOddFloat)
TEST_ADD(GpuVecRandomTest, testNormalOddDouble)
TEST_END(GpuVecRandomTest)

#endif
