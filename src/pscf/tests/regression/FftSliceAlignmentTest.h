#ifndef PSCF_FFT_SLICE_ALIGNMENT_TEST_H
#define PSCF_FFT_SLICE_ALIGNMENT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/chem/Edge.h>
#include <pscf/chem/PolymerModel.h>
#include <pscf/chem/PolymerType.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/field/cpp/FFT.h>
#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldDft.h>
#include <pscf/field/cpp/WaveList.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/containers/FSArray.h>
#include <util/misc/Exception.h>
#include <cmath>
#include <sstream>

using namespace Util;
using namespace Pscf;

class FftSliceAlignmentTest : public UnitTest
{
public:
   void setUp() {}
   void tearDown() { PolymerModel::setModel(PolymerModel::Thread); }

   // Rounded right angles are physical perturbations, not exact ties.
   // Every inverse pair must still have the same minimum squared norm.

   template <int D>
   void checkFftSlice(IntVec<D> const& dims)
   {
      int n = 1;
      for (int j = 0; j < D; ++j) n *= dims[j];
      BackendArray<double,CPT> source(n+2), result(n+2);
      for (int i = 0; i < n+2; ++i) {
         source[i] = -79.0;
         result[i] = 83.0;
      }
      RField<D,CPT> slice, back, aligned, alignedBack;
      slice.associate(source, 1, dims);
      back.associate(result, 1, dims);
      aligned.allocate(dims);
      alignedBack.allocate(dims);
      for (int i = 0; i < n; ++i) {
         aligned[i] = 5.0*cos(0.37*i) + sin(0.013*i*i);
         slice[i] = aligned[i];
      }
      RFieldDft<D,CPT> ref, transformed, copy;
      ref.allocate(dims);
      transformed.allocate(dims);
      copy.allocate(dims);
      FFT<D,CPT> fft;
      fft.setup(dims);
      fft.forwardTransform(aligned, ref);
      fft.forwardTransform(slice, transformed);
      for (int i = 0; i < ref.capacity(); ++i) {
         TEST_ASSERT(std::abs(ref[i][0] - transformed[i][0]) < 1e-12);
         TEST_ASSERT(std::abs(ref[i][1] - transformed[i][1]) < 1e-12);
      }
      copy = transformed;
      fft.inverseTransformSafe(transformed, back);
      for (int i = 0; i < copy.capacity(); ++i) {
         TEST_ASSERT(transformed[i][0] == copy[i][0]);
         TEST_ASSERT(transformed[i][1] == copy[i][1]);
      }
      fft.inverseTransformUnsafe(copy, alignedBack);
      for (int i = 0; i < n; ++i) {
         TEST_ASSERT(slice[i] == aligned[i]);
         TEST_ASSERT(std::abs(back[i] - aligned[i]) < 1e-11);
         TEST_ASSERT(std::abs(alignedBack[i] - aligned[i]) < 1e-11);
      }
      copy = transformed;
      fft.inverseTransformUnsafe(copy, back);
      for (int i = 0; i < n; ++i) {
         TEST_ASSERT(std::abs(back[i] - aligned[i]) < 1e-11);
      }
      TEST_ASSERT(source[0] == -79.0 && source[n+1] == -79.0);
      TEST_ASSERT(result[0] == 83.0 && result[n+1] == 83.0);
   }

   void testCpuFftSlices1D()
   {
      printMethod(TEST_FUNC);
      IntVec<1> dims;
      dims[0] = 64; checkFftSlice(dims);
      dims[0] = 63; checkFftSlice(dims);
   }

   void testCpuFftSlices2D3D()
   {
      printMethod(TEST_FUNC);
      IntVec<2> d2;
      d2[0] = 8; d2[1] = 8; checkFftSlice(d2);
      d2[0] = 7; d2[1] = 9; checkFftSlice(d2);
      IntVec<3> d3;
      d3[0] = 4; d3[1] = 4; d3[2] = 4; checkFftSlice(d3);
      d3[0] = 3; d3[1] = 5; d3[2] = 7; checkFftSlice(d3);
   }
};

TEST_BEGIN(FftSliceAlignmentTest)
TEST_ADD(FftSliceAlignmentTest, testCpuFftSlices1D)
TEST_ADD(FftSliceAlignmentTest, testCpuFftSlices2D3D)
TEST_END(FftSliceAlignmentTest)

#endif
