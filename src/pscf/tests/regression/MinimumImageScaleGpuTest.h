#ifndef PSCF_MINIMUM_IMAGE_SCALE_GPU_TEST_H
#define PSCF_MINIMUM_IMAGE_SCALE_GPU_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/field/cuda/WaveList.h>
#include <pscf/backend/cuda/HostRecvArray.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/containers/FSArray.h>
#include <util/math/Constants.h>
#include <cmath>

using namespace Util;
using namespace Pscf;

class MinimumImageScaleGpuTest : public UnitTest
{
public:
   void setUp() {}
   void tearDown() {}

   template <int D>
   void checkOrthogonal(UnitCell<D> const& cell, double length)
   {
      IntVec<D> dims(8);
      Mesh<D> mesh(dims);
      WaveList<D,CUT> waves;
      waves.allocate(mesh, cell);
      waves.computeMinimumImages();
      HostRecvArray<int,CUT> images;
      HostRecvArray<cudaReal,CUT> values;
      images = waves.minImages_d();
      values = waves.kSq();
      TEST_ASSERT(values[0] == 0.0);
      IntVec<D> kdims = dims;
      kdims[D-1] = dims[D-1]/2 + 1;
      int ksize = values.capacity();
      MeshIterator<D> iter(kdims);
      #ifdef SINGLE_PRECISION
      const double tol = 2e-5;
      #else
      const double tol = 1e-12;
      #endif
      for (iter.begin(); !iter.atEnd(); ++iter) {
         IntVec<D> v = iter.position();
         double norm = 0.0;
         for (int j = 0; j < D; ++j) {
            int expected = v[j] > dims[j]/2 ? v[j]-dims[j] : v[j];
            TEST_ASSERT(images[iter.rank()+j*ksize] == expected);
            norm += double(expected)*expected;
         }
         double exact = 4*Constants::Pi*Constants::Pi*norm;
         double actual = values[iter.rank()]*length*length;
         TEST_ASSERT(std::abs(actual-exact) < tol*(1.0+exact));
      }
   }

   void testMinimumImages2DUnits()
   {
      printMethod(TEST_FUNC);
      double lengths[] = {1e-6, 1.0, 1e6};
      for (double length : lengths) {
         FSArray<double,6> p;
         p.append(length);
         UnitCell<2> cell;
         cell.set(UnitCell<2>::Square, p);
         checkOrthogonal(cell, length);
      }
   }

   void testMinimumImages3DUnits()
   {
      printMethod(TEST_FUNC);
      double lengths[] = {1e-6, 1.0, 1e6};
      for (double length : lengths) {
         FSArray<double,6> p;
         p.append(length);
         UnitCell<3> cell;
         cell.set(UnitCell<3>::Cubic, p);
         checkOrthogonal(cell, length);
      }
   }
};

TEST_BEGIN(MinimumImageScaleGpuTest)
TEST_ADD(MinimumImageScaleGpuTest, testMinimumImages2DUnits)
TEST_ADD(MinimumImageScaleGpuTest, testMinimumImages3DUnits)
TEST_END(MinimumImageScaleGpuTest)

#endif
