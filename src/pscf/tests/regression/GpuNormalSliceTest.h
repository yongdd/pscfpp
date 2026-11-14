#ifndef PSCF_GPU_NORMAL_SLICE_TEST_H
#define PSCF_GPU_NORMAL_SLICE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/backend/cuda/VecRandom.h>
#include <pscf/backend/cuda/HostSendArray.h>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Normal generation must support odd-offset slices and preserve the stream.
*/
class GpuNormalSliceTest : public UnitTest
{
public:

   void setUp()
   {  setVerbose(0); }

   void tearDown()
   {}

   template <typename Data>
   void checkSlices()
   {
      VecRandom<CUT> referenceRandom;
      VecRandom<CUT> sliceRandom;
      referenceRandom.setSeed(123456789);
      sliceRandom.setSeed(123456789);
      const Data sentinel = Data(-4321.5);

      // Grow and then reuse larger alignment storage for smaller requests.
      // A count1 request also exercises the separate odd-tail pair alone.
      const int sizes[] = {1, 2, 3, 10, 2, 5, 1, 3};
      const int nSize = sizeof(sizes)/sizeof(sizes[0]);
      for (int s = 0; s < nSize; ++s) {
         const int n = sizes[s];
         for (int offset = 0; offset < 3; ++offset) {
            const Data stddev = Data(0.25 + 0.125*(s % 3));
            const Data mean = Data(-0.5 + 0.25*offset);
            BackendArray<Data,CUT> reference(n), parent(n+4), slice;
            HostSendArray<Data,CUT> expected(n), actual(n+4);
            for (int i = 0; i < actual.capacity(); ++i) {
               actual[i] = sentinel;
            }
            parent = actual;
            slice.associate(parent, offset, n);

            // Identical requests must consume identical streams regardless
            // of destination alignment. Neither generator is reseeded here.
            referenceRandom.normal(reference, stddev, mean);
            sliceRandom.normal(slice, stddev, mean);
            expected = reference;
            actual = parent;
            for (int i = 0; i < actual.capacity(); ++i) {
               if (i < offset || i >= offset+n) {
                  TEST_ASSERT(actual[i] == sentinel);
               } else {
                  TEST_ASSERT(std::isfinite(actual[i]));
                  TEST_ASSERT(actual[i] == expected[i-offset]);
               }
            }
         }
      }

      // Extra normal draws from scratch padding would change this stream.
      const int nUniform = 9;
      BackendArray<Data,CUT> reference(nUniform), parent(nUniform+2), slice;
      HostSendArray<Data,CUT> expected(nUniform), actual(nUniform+2);
      for (int i = 0; i < actual.capacity(); ++i) {
         actual[i] = sentinel;
      }
      parent = actual;
      slice.associate(parent, 1, nUniform);
      referenceRandom.uniform(reference);
      sliceRandom.uniform(slice);
      expected = reference;
      actual = parent;
      TEST_ASSERT(actual[0] == sentinel);
      TEST_ASSERT(actual[nUniform+1] == sentinel);
      for (int i = 0; i < nUniform; ++i) {
         TEST_ASSERT(actual[i+1] == expected[i]);
         TEST_ASSERT(actual[i+1] > Data(0));
         TEST_ASSERT(actual[i+1] <= Data(1));
      }
   }

   void testFloatSlices()
   {
      printMethod(TEST_FUNC);
      checkSlices<float>();
   }

   void testDoubleSlices()
   {
      printMethod(TEST_FUNC);
      checkSlices<double>();
   }
};

TEST_BEGIN(GpuNormalSliceTest)
TEST_ADD(GpuNormalSliceTest, testFloatSlices)
TEST_ADD(GpuNormalSliceTest, testDoubleSlices)
TEST_END(GpuNormalSliceTest)

#endif
