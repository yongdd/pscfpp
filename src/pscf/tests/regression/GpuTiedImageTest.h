#ifndef PSCF_GPU_TIED_IMAGE_TEST_H
#define PSCF_GPU_TIED_IMAGE_TEST_H

#include "CoreTiedImageTest.h"

#include <pscf/field/cuda/WaveList.h>
#include <pscf/backend/cuda/cudaTypes.h>
#include <pscf/backend/cuda/HostRecvArray.h>

/*
* Derivatives dKSq of wavevectors with tied minimum images, GPU WaveList
* (see CoreTiedImageTest): the weighted sums of dKSq of a symmetric cell
* are isotropic, and images, kSq and dKSq equal those of the CPU WaveList.
*/
class GpuTiedImageTest : public CoreTiedImageTest
{

public:

   /*
   * Deviatoric residual of the weighted sums of GPU dKSq (real fields),
   * as CoreTiedImageTest::residual, after comparing the GPU WaveList
   * with the CPU WaveList.
   */
   template <int D>
   double gpuResidual(UnitCell<D> const & cell, IntVec<D> const & d,
                      double const * dLogV, int const * params, int nCheck)
   {
      Mesh<D> mesh(d);
      WaveList<D,CUT> gpu(true);
      WaveList<D,CPT> cpu(true);
      gpu.allocate(mesh, cell);
      cpu.allocate(mesh, cell);
      gpu.computeMinimumImages();
      gpu.computedKSq();
      cpu.computeMinimumImages();
      cpu.computedKSq();
      int kSize = cpu.kSize();
      DArray< IntVec<D> > const & images = gpu.minImages_h();
      HostRecvArray<cudaReal,CUT> kSq(gpu.kSq());
      double mean = 0.0;
      for (int i = 0; i < kSize; ++i) {
         TEST_ASSERT(images[i] == cpu.minImages()[i]);
         double scale = 1.0 + std::abs(cpu.kSq()[i]);
         TEST_ASSERT(std::abs(kSq[i] - cpu.kSq()[i]) < 1.0E-10*scale);
         mean += kSq[i];
      }
      mean /= double(kSize);
      double S[6];
      double scale = 0.0;
      for (int n = 0; n < cell.nParameter(); ++n) {
         HostRecvArray<cudaReal,CUT> dKSq(gpu.dKSq(n));
         S[n] = 0.0;
         for (int i = 0; i < kSize; ++i) {
            double s = 1.0 + std::abs(cpu.kSq()[i]);
            TEST_ASSERT(std::abs(dKSq[i] - cpu.dKSq(n)[i]) < 1.0E-10*s);
            S[n] += std::exp(-kSq[i]/mean)*dKSq[i];
         }
         scale = std::max(scale, std::abs(S[n]));
      }
      int n0 = params[0];
      double maxResidual = 0.0;
      for (int k = 1; k < nCheck; ++k) {
         int n = params[k];
         double r = S[n] - S[n0]*dLogV[n]/dLogV[n0];
         maxResidual = std::max(maxResidual, std::abs(r)/scale);
      }
      return maxResidual;
   }

   void testGpuSymmetricCells()
   {
      printMethod(TEST_FUNC);
      double third = 2.0*Constants::Pi/3.0;
      double tolerance = 1.0E-12;
      int params2[2] = {0, 1};

      // Rhombic at 120 degrees, even and odd meshes
      double a = 1.7;
      double p[2] = {a, third};
      UnitCell<2> rhombic;
      makeCell<2>(rhombic, UnitCell<2>::Rhombic, p, 2);
      double dLogV[2];
      dLogV[0] = 2.0/a;
      dLogV[1] = std::cos(third)/std::sin(third);
      TEST_ASSERT(gpuResidual<2>(rhombic, IntVec<2>(12), dLogV, params2, 2)
                  < tolerance);
      TEST_ASSERT(gpuResidual<2>(rhombic, IntVec<2>(13), dLogV, params2, 2)
                  < tolerance);

      // Rhombohedral at 60 degrees (fcc)
      double ph[2] = {2.0, 0.5*third};
      UnitCell<3> rhombohedral;
      makeCell<3>(rhombohedral, UnitCell<3>::Rhombohedral, ph, 2);
      double c = std::cos(ph[1]);
      dLogV[0] = 3.0/ph[0];
      dLogV[1] = 3.0*std::sin(ph[1])*c*(1.0 - c)/(1.0 - 3.0*c*c + 2.0*c*c*c);
      TEST_ASSERT(gpuResidual<3>(rhombohedral, IntVec<3>(8), dLogV,
                                 params2, 2) < tolerance);
   }

};

TEST_BEGIN(GpuTiedImageTest)
TEST_ADD(GpuTiedImageTest, testGpuSymmetricCells)
TEST_END(GpuTiedImageTest)

#endif
