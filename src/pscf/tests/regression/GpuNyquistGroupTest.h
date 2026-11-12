#ifndef PSCF_GPU_NYQUIST_GROUP_TEST_H
#define PSCF_GPU_NYQUIST_GROUP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/field/cuda/WaveList.h>
#include <pscf/field/cpp/WaveList.h>
#include <pscf/backend/cuda/cudaTypes.h>
#include <pscf/backend/cuda/HostRecvArray.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

#include <util/containers/FSArray.h>
#include <util/math/Constants.h>

#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* WaveList<D,CUT>::allowNyquistAverage: with averaging forbidden, the
* GPU images, kSq and dKSq of a variable-angle lattice equal those of
* the CPU WaveList with averaging forbidden (minimum-image values), and
* allowing it again restores the averaged values.
*/
class GpuNyquistGroupTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   template <int D>
   void compare(WaveList<D,CUT>& gpu, WaveList<D,CPT>& cpu, int nParam)
   {
      gpu.computeMinimumImages();
      gpu.computeKSq();
      gpu.computedKSq();
      cpu.computeMinimumImages();
      cpu.computeKSq();
      cpu.computedKSq();
      int kSize = cpu.kSq().capacity();
      DArray< IntVec<D> > const & images = gpu.minImages_h();
      HostRecvArray<cudaReal,CUT> kSq(gpu.kSq());
      TEST_ASSERT(kSq.capacity() == kSize);
      for (int i = 0; i < kSize; ++i) {
         TEST_ASSERT(images[i] == cpu.minImages()[i]);
         double scale = 1.0 + std::abs(cpu.kSq()[i]);
         TEST_ASSERT(std::abs(kSq[i] - cpu.kSq()[i]) < 1.0E-10*scale);
      }
      for (int n = 0; n < nParam; ++n) {
         HostRecvArray<cudaReal,CUT> dKSq(gpu.dKSq(n));
         for (int i = 0; i < kSize; ++i) {
            double scale = 1.0 + std::abs(cpu.kSq()[i]);
            TEST_ASSERT(std::abs(dKSq[i] - cpu.dKSq(n)[i])
                        < 1.0E-10*scale);
         }
      }
   }

   template <int D>
   void checkCell(UnitCell<D> const & cell, IntVec<D> const & d)
   {
      Mesh<D> mesh(d);
      WaveList<D,CUT> gpu(true);
      WaveList<D,CPT> cpuAverage(true);
      WaveList<D,CPT> cpuMinimum(true);
      gpu.allocate(mesh, cell);
      cpuAverage.allocate(mesh, cell);
      cpuMinimum.allocate(mesh, cell);
      cpuMinimum.allowNyquistAverage(false);
      TEST_ASSERT(gpu.hasNyquistAverage());
      compare<D>(gpu, cpuAverage, cell.nParameter());

      gpu.allowNyquistAverage(false);
      TEST_ASSERT(!gpu.hasNyquistAverage());
      TEST_ASSERT(!gpu.hasMinImages());
      TEST_ASSERT(!gpu.hasKSq());
      TEST_ASSERT(!gpu.hasdKSq());
      compare<D>(gpu, cpuMinimum, cell.nParameter());

      gpu.allowNyquistAverage(true);
      TEST_ASSERT(gpu.hasNyquistAverage());
      compare<D>(gpu, cpuAverage, cell.nParameter());
   }

   void testAllowNyquistAverage()
   {
      printMethod(TEST_FUNC);
      double twoThirdsPi = 2.0*Constants::Pi/3.0;

      FSArray<double, 6> p2;
      p2.append(2.0);
      p2.append(twoThirdsPi);
      UnitCell<2> rhombic;
      rhombic.set(UnitCell<2>::Rhombic, p2);
      checkCell<2>(rhombic, IntVec<2>(8));

      FSArray<double, 6> p3;
      p3.append(2.0);
      p3.append(2.0);
      p3.append(1.7);
      p3.append(0.0);
      p3.append(0.0);
      p3.append(twoThirdsPi);
      UnitCell<3> tric;
      tric.set(UnitCell<3>::Triclinic, p3);
      IntVec<3> d(8);
      d[2] = 6;
      checkCell<3>(tric, d);
   }

};

TEST_BEGIN(GpuNyquistGroupTest)
TEST_ADD(GpuNyquistGroupTest, testAllowNyquistAverage)
TEST_END(GpuNyquistGroupTest)

#endif
