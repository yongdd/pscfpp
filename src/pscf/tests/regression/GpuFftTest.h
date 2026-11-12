#ifndef PSCF_GPU_FFT_TEST_H
#define PSCF_GPU_FFT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/field/cuda/FFT.h>
#include <pscf/field/cuda/RField.h>
#include <pscf/field/cuda/RFieldDft.h>

#include <pscf/field/cpp/FFT.h>
#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldDft.h>

#include <pscf/backend/cuda/BackendArray.h>
#include <pscf/backend/cuda/HostRecvArray.h>
#include <pscf/backend/cuda/HostSendArray.h>
#include <pscf/math/IntVec.h>

#include <util/containers/DArray.h>
#include <util/misc/Exception.h>

#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests of the GPU FFT (FFT<D,CUT>).
*/
class GpuFftTest : public UnitTest 
{

private:

   // Error tolerance for equality
   #ifdef SINGLE_PRECISION
   constexpr static cudaReal tolerance_ = 1E-5;
   #else
   #ifdef DOUBLE_PRECISION
   constexpr static cudaReal tolerance_ = 1E-10;
   #endif
   #endif

public:

   void setUp() {}
   void tearDown() {}

   void testTransformRealOdd1D();
   void testTransformRealOdd2D();
   void testTransformRealOdd3D();

   template <int D>
   void checkTransformRealOdd(IntVec<D> const & d);

};

void GpuFftTest::testTransformRealOdd1D()
{
   printMethod(TEST_FUNC);
   IntVec<1> d;
   d[0] = 15;
   checkTransformRealOdd<1>(d);
}

void GpuFftTest::testTransformRealOdd2D()
{
   printMethod(TEST_FUNC);
   IntVec<2> d;
   d[0] = 9;
   d[1] = 7;
   checkTransformRealOdd<2>(d);
}

void GpuFftTest::testTransformRealOdd3D()
{
   printMethod(TEST_FUNC);
   IntVec<3> d;
   d[0] = 5;
   d[1] = 7;
   d[2] = 9;
   checkTransformRealOdd<3>(d);
}

/*
* Real transforms on a mesh with an odd number of points, compared with
* the CPU FFT.
*
* The field is transformed as a separate array, and as the middle one of
* three consecutive fields in one array, like one contour step q(i) of a
* propagator. That field starts rSize elements after the beginning of the
* array, which for odd rSize is not aligned to the size of a complex
* element. The other two fields hold guard values that must not change.
*/
template <int D>
void GpuFftTest::checkTransformRealOdd(IntVec<D> const & d)
{
   int rSize = 1;
   for (int i = 0; i < D; ++i) {
      rSize *= d[i];
   }
   TEST_ASSERT(rSize % 2 == 1);

   // Input data, on host
   DArray<double> data;
   data.allocate(rSize);
   double x;
   for (int i = 0; i < rSize; ++i) {
      x = double(i);
      data[i] = 0.4 + cos(0.37*x) + 0.5*sin(0.013*x*x);
   }

   // Reference transform, by the CPU FFT
   Pscf::FFT<D,CPT> cpuFft;
   cpuFft.setup(d);
   Pscf::RField<D,CPT> rCpu;
   rCpu.allocate(d);
   for (int i = 0; i < rSize; ++i) {
      rCpu[i] = data[i];
   }
   Pscf::RFieldDft<D,CPT> kCpu;
   kCpu.allocate(d);
   cpuFft.forwardTransform(rCpu, kCpu);
   int kSize = kCpu.capacity();

   // Host arrays: input field, three consecutive fields (input field
   // in the middle), and three fields with zero in the middle
   HostSendArray<cudaReal,CUT> rField_h(rSize);
   HostSendArray<cudaReal,CUT> slices_h(3*rSize);
   HostSendArray<cudaReal,CUT> slicesOut_h(3*rSize);
   for (int i = 0; i < rSize; ++i) {
      rField_h[i] = data[i];
      slices_h[i] = -2.0 - 0.001*i;
      slices_h[rSize + i] = data[i];
      slices_h[2*rSize + i] = 3.0 + 0.001*i;
      slicesOut_h[i] = slices_h[i];
      slicesOut_h[rSize + i] = 0.0;
      slicesOut_h[2*rSize + i] = slices_h[2*rSize + i];
   }
   HostSendArray<cudaComplex,CUT> kCpu_h(kSize);
   for (int i = 0; i < kSize; ++i) {
      kCpu_h[i].x = kCpu[i][0];
      kCpu_h[i].y = kCpu[i][1];
   }

   // Device arrays
   Pscf::RField<D,CUT> rField(d);
   rField = rField_h;
   BackendArray<cudaReal,CUT> slices(3*rSize);
   slices = slices_h;
   Pscf::RField<D,CUT> slice;
   slice.associate(slices, rSize, d);
   BackendArray<cudaReal,CUT> slicesOut(3*rSize);
   slicesOut = slicesOut_h;
   Pscf::RField<D,CUT> sliceOut;
   sliceOut.associate(slicesOut, rSize, d);
   Pscf::RField<D,CUT> rFieldOut(d);
   Pscf::RFieldDft<D,CUT> kField(d);
   Pscf::RFieldDft<D,CUT> kSlice(d);
   Pscf::RFieldDft<D,CUT> kIn(d);
   kIn = kCpu_h;
   TEST_ASSERT(kField.capacity() == kSize);

   Pscf::FFT<D,CUT> v;
   v.setup(d);

   // Forward and inverse transforms of the separate and middle fields
   bool thrown = false;
   try {
      v.forwardTransform(rField, kField);
      v.forwardTransform(slice, kSlice);
      v.inverseTransformSafe(kIn, rFieldOut);
      v.inverseTransformUnsafe(kIn, sliceOut);
   } catch (Exception& e) {
      thrown = true;
   }
   TEST_ASSERT(!thrown);

   // Forward transforms equal the CPU transform
   HostRecvArray<cudaComplex,CUT> k_h(kSize);
   k_h = kField;
   for (int i = 0; i < kSize; ++i) {
      TEST_ASSERT(fabs(k_h[i].x - kCpu[i][0]) < tolerance_);
      TEST_ASSERT(fabs(k_h[i].y - kCpu[i][1]) < tolerance_);
   }
   k_h = kSlice;
   for (int i = 0; i < kSize; ++i) {
      TEST_ASSERT(fabs(k_h[i].x - kCpu[i][0]) < tolerance_);
      TEST_ASSERT(fabs(k_h[i].y - kCpu[i][1]) < tolerance_);
   }

   // Inverse transforms of the CPU transform recover the input field,
   // and the input of the forward transform is unchanged
   HostRecvArray<cudaReal,CUT> r_h(rSize);
   r_h = rFieldOut;
   for (int i = 0; i < rSize; ++i) {
      TEST_ASSERT(fabs(r_h[i] - data[i]) < tolerance_);
   }
   HostRecvArray<cudaReal,CUT> slices2_h(3*rSize);
   slices2_h = slices;
   for (int i = 0; i < 3*rSize; ++i) {
      TEST_ASSERT(slices2_h[i] == slices_h[i]);
   }

   // The inverse transform writes the middle field only
   slices2_h = slicesOut;
   for (int i = 0; i < rSize; ++i) {
      TEST_ASSERT(slices2_h[i] == slicesOut_h[i]);
      TEST_ASSERT(fabs(slices2_h[rSize + i] - data[i]) < tolerance_);
      TEST_ASSERT(slices2_h[2*rSize + i] == slicesOut_h[2*rSize + i]);
   }
}

TEST_BEGIN(GpuFftTest)
TEST_ADD(GpuFftTest, testTransformRealOdd1D)
TEST_ADD(GpuFftTest, testTransformRealOdd2D)
TEST_ADD(GpuFftTest, testTransformRealOdd3D)
TEST_END(GpuFftTest)

#endif
