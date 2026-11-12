#ifndef PSCF_GPU_NYQUIST_TEST_H
#define PSCF_GPU_NYQUIST_TEST_H

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
* Wavevectors with Nyquist indices (even meshes) in the GPU WaveList:
* images, kSq and dKSq equal those of the CPU WaveList (kSq averaged
* over the signs of Nyquist indices for lattice systems with variable
* angles), and dKSq is the derivative of kSq at right angles.
*/
class GpuNyquistTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   template <int D>
   void makeCell(UnitCell<D>& cell,
                 typename UnitCell<D>::LatticeSystem lattice,
                 double const * p, int nParam)
   {
      FSArray<double, 6> parameters;
      for (int i = 0; i < nParam; ++i) {
         parameters.append(p[i]);
      }
      cell.set(lattice, parameters);
   }

   /*
   * Compare GPU and CPU WaveList (real fields): images, kSq and dKSq,
   * first after allocation and then after a change of the cell.
   */
   template <int D>
   void compareCpu(UnitCell<D>& cell, IntVec<D> const & d)
   {
      Mesh<D> mesh(d);
      WaveList<D,CUT> gpu(true);
      WaveList<D,CPT> cpu(true);
      gpu.allocate(mesh, cell);
      cpu.allocate(mesh, cell);
      for (int pass = 0; pass < 2; ++pass) {
         if (pass == 1) {
            FSArray<double, 6> p = cell.parameters();
            for (int n = 0; n < p.size(); ++n) {
               p[n] *= 1.01;
            }
            cell.setParameters(p);
            gpu.clearUnitCellData();
            cpu.clearUnitCellData();
         }
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
         for (int n = 0; n < cell.nParameter(); ++n) {
            HostRecvArray<cudaReal,CUT> dKSq(gpu.dKSq(n));
            for (int i = 0; i < kSize; ++i) {
               double scale = 1.0 + std::abs(cpu.kSq()[i]);
               TEST_ASSERT(std::abs(dKSq[i] - cpu.dKSq(n)[i])
                           < 1.0E-10*scale);
            }
         }
      }
   }

   void testCompareCpu()
   {
      printMethod(TEST_FUNC);
      double halfPi = 0.5*Constants::Pi;
      IntVec<3> d;
      d[0] = 8; d[1] = 6; d[2] = 10;
      IntVec<3> d8(8);

      double pm[4] = {1.9, 2.1, 2.3, halfPi};
      UnitCell<3> mono;
      makeCell<3>(mono, UnitCell<3>::Monoclinic, pm, 4);
      compareCpu<3>(mono, d);
      pm[3] = 2.0;
      makeCell<3>(mono, UnitCell<3>::Monoclinic, pm, 4);
      compareCpu<3>(mono, d);

      double pr[2] = {2.0, halfPi};
      UnitCell<3> rhomb;
      makeCell<3>(rhomb, UnitCell<3>::Rhombohedral, pr, 2);
      compareCpu<3>(rhomb, d8);
      pr[1] = 1.05;
      makeCell<3>(rhomb, UnitCell<3>::Rhombohedral, pr, 2);
      compareCpu<3>(rhomb, d8);

      double pt[6] = {1.9, 2.1, 2.3, 1.3, 0.5, 1.2};
      UnitCell<3> tric;
      makeCell<3>(tric, UnitCell<3>::Triclinic, pt, 6);
      compareCpu<3>(tric, d);

      double po[3] = {1.9, 2.3, 2.2};
      UnitCell<2> obl;
      makeCell<2>(obl, UnitCell<2>::Oblique, po, 3);
      IntVec<2> d2;
      d2[0] = 8; d2[1] = 6;
      compareCpu<2>(obl, d2);

      double pc[3] = {1.9, 2.1, 2.3};
      UnitCell<3> orth;
      makeCell<3>(orth, UnitCell<3>::Orthorhombic, pc, 3);
      compareCpu<3>(orth, d);
   }

   /*
   * Monoclinic cell at beta = 90 degrees on an even mesh: GPU dKSq
   * for beta equals the central finite difference of GPU kSq.
   */
   void testSmoothAtRightAngle()
   {
      printMethod(TEST_FUNC);
      double halfPi = 0.5*Constants::Pi;
      double h = 1.0E-6;
      IntVec<3> d;
      d[0] = 8; d[1] = 6; d[2] = 10;
      Mesh<3> mesh(d);

      double p[4] = {1.9, 2.1, 2.3, halfPi};
      UnitCell<3> cell, cp, cm;
      makeCell<3>(cell, UnitCell<3>::Monoclinic, p, 4);
      p[3] = halfPi + h;
      makeCell<3>(cp, UnitCell<3>::Monoclinic, p, 4);
      p[3] = halfPi - h;
      makeCell<3>(cm, UnitCell<3>::Monoclinic, p, 4);

      WaveList<3,CUT> w0(true), wp(true), wm(true);
      w0.allocate(mesh, cell);
      wp.allocate(mesh, cp);
      wm.allocate(mesh, cm);
      w0.computedKSq();
      wp.computeKSq();
      wm.computeKSq();
      HostRecvArray<cudaReal,CUT> kSq0(w0.kSq());
      HostRecvArray<cudaReal,CUT> kSqP(wp.kSq());
      HostRecvArray<cudaReal,CUT> kSqM(wm.kSq());
      HostRecvArray<cudaReal,CUT> dKSq(w0.dKSq(3));
      DArray<bool> const & implicit = w0.implicitInverse();
      double maxError = 0.0;
      for (int i = 0; i < kSq0.capacity(); ++i) {
         double fd = (kSqP[i] - kSqM[i])/(2.0*h);
         if (implicit[i]) fd *= 2.0;
         double scale = std::max(1.0, std::abs(kSq0[i]));
         double error = std::abs(dKSq[i] - fd)/scale;
         if (error > maxError) maxError = error;
      }
      TEST_ASSERT(maxError < 1.0E-6);
   }

};

TEST_BEGIN(GpuNyquistTest)
TEST_ADD(GpuNyquistTest, testCompareCpu)
TEST_ADD(GpuNyquistTest, testSmoothAtRightAngle)
TEST_END(GpuNyquistTest)

#endif
