#ifndef PSCF_GPU_NYQUIST_WEIGHT_TEST_H
#define PSCF_GPU_NYQUIST_WEIGHT_TEST_H

#include "CoreNyquistWeightTest.h"

#include <pscf/field/cuda/WaveList.h>
#include <pscf/backend/cuda/cudaTypes.h>
#include <pscf/backend/cuda/HostRecvArray.h>

#include <vector>

/*
* Average of kSq over the signs of Nyquist indices, GPU WaveList (see
* CoreNyquistWeightTest): kSq is equal for wavevectors related by the
* symmetry operations of the lattice, and images, kSq and dKSq equal
* those of the CPU WaveList, also in the transition between the average
* (near orthogonal geometries) and the minimum image.
*/
class GpuNyquistWeightTest : public CoreNyquistWeightTest
{

public:

   /*
   * Compare a GPU WaveList with a CPU WaveList (complex fields): images
   * exactly, kSq and dKSq relative to 1 + |kSq|. Returns the largest
   * difference and sets kSq to the GPU values.
   */
   template <int D>
   double compareCpu(UnitCell<D> const & cell, IntVec<D> const & d,
                     std::vector<double>& kSq)
   {
      Mesh<D> mesh(d);
      WaveList<D,CUT> gpu(false);
      WaveList<D,CPT> cpu(false);
      gpu.allocate(mesh, cell);
      cpu.allocate(mesh, cell);
      gpu.computeMinimumImages();
      gpu.computedKSq();
      cpu.computeMinimumImages();
      cpu.computedKSq();
      int kSize = cpu.kSize();
      DArray< IntVec<D> > const & images = gpu.minImages_h();
      HostRecvArray<cudaReal,CUT> kSqGpu(gpu.kSq());
      kSq.resize(kSize);
      double maxDiff = 0.0;
      for (int i = 0; i < kSize; ++i) {
         TEST_ASSERT(images[i] == cpu.minImages()[i]);
         kSq[i] = kSqGpu[i];
         double scale = 1.0 + std::abs(cpu.kSq()[i]);
         maxDiff = std::max(maxDiff, std::abs(kSq[i] - cpu.kSq()[i])/scale);
      }
      for (int n = 0; n < cell.nParameter(); ++n) {
         HostRecvArray<cudaReal,CUT> dKSq(gpu.dKSq(n));
         for (int i = 0; i < kSize; ++i) {
            double scale = 1.0 + std::abs(cpu.kSq()[i]);
            maxDiff = std::max(maxDiff,
                               std::abs(dKSq[i] - cpu.dKSq(n)[i])/scale);
         }
      }
      return maxDiff;
   }

   /*
   * Largest relative difference of GPU kSq between wavevectors related
   * by symmetry operations of the lattice.
   */
   template <int D>
   double gpuOrbitDifference(UnitCell<D> const & cell, IntVec<D> const & d)
   {
      std::vector<double> kSq;
      TEST_ASSERT(compareCpu<D>(cell, d, kSq) < 1.0E-10);
      std::vector< std::vector<int> > ops;
      TEST_ASSERT(latticeOperations<D>(cell, d, ops) > 0);
      Mesh<D> mesh(d);
      double maxDiff = 0.0;
      MeshIterator<D> iter(d);
      IntVec<D> w;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         IntVec<D> v = iter.position();
         double kSqV = kSq[iter.rank()];
         for (size_t k = 0; k < ops.size(); ++k) {
            for (int j = 0; j < D; ++j) {
               int s = 0;
               for (int i = 0; i < D; ++i) {
                  s += v[i]*ops[k][D*i + j];
               }
               s = s % d[j];
               if (s < 0) s += d[j];
               w[j] = s;
            }
            double diff = std::abs(kSq[mesh.rank(w)] - kSqV)/(1.0 + kSqV);
            maxDiff = std::max(maxDiff, diff);
         }
      }
      return maxDiff;
   }

   void testGpuSymmetricCells()
   {
      printMethod(TEST_FUNC);
      double third = 2.0*Constants::Pi/3.0;
      double tolerance = 1.0E-10;

      double pr[2] = {1.7, third};
      UnitCell<2> rhombic;
      makeCell<2>(rhombic, UnitCell<2>::Rhombic, pr, 2);
      TEST_ASSERT(gpuOrbitDifference<2>(rhombic, IntVec<2>(12))
                  < tolerance);

      double ph[2] = {2.0, 0.5*third};
      UnitCell<3> rhombohedral;
      makeCell<3>(rhombohedral, UnitCell<3>::Rhombohedral, ph, 2);
      TEST_ASSERT(gpuOrbitDifference<3>(rhombohedral, IntVec<3>(6))
                  < tolerance);

      double pm[4] = {1.7, 2.1, 1.7, third};
      UnitCell<3> monoclinic;
      makeCell<3>(monoclinic, UnitCell<3>::Monoclinic, pm, 4);
      IntVec<3> dm;
      dm[0] = 6; dm[1] = 4; dm[2] = 6;
      TEST_ASSERT(gpuOrbitDifference<3>(monoclinic, dm) < tolerance);
   }

   void testGpuTransition()
   {
      printMethod(TEST_FUNC);
      double tolerance = 1.0E-10;
      std::vector<double> kSq;

      double pr[2] = {1.7, std::acos(-0.175)};
      UnitCell<2> rhombic;
      makeCell<2>(rhombic, UnitCell<2>::Rhombic, pr, 2);
      TEST_ASSERT(compareCpu<2>(rhombic, IntVec<2>(10), kSq) < tolerance);

      double pm[4] = {1.9, 2.1, 2.3, std::acos(-0.175)};
      UnitCell<3> monoclinic;
      makeCell<3>(monoclinic, UnitCell<3>::Monoclinic, pm, 4);
      IntVec<3> dm;
      dm[0] = 8; dm[1] = 6; dm[2] = 10;
      TEST_ASSERT(compareCpu<3>(monoclinic, dm, kSq) < tolerance);

      double pt[6] = {1.9, 2.1, 2.3, 0.4, 0.15, 1.62};
      UnitCell<3> triclinic;
      makeCell<3>(triclinic, UnitCell<3>::Triclinic, pt, 6);
      IntVec<3> dt;
      dt[0] = 6; dt[1] = 8; dt[2] = 10;
      TEST_ASSERT(compareCpu<3>(triclinic, dt, kSq) < tolerance);

      double ph[2] = {2.0, 1.3};
      UnitCell<3> rhombohedral;
      makeCell<3>(rhombohedral, UnitCell<3>::Rhombohedral, ph, 2);
      TEST_ASSERT(compareCpu<3>(rhombohedral, IntVec<3>(8), kSq)
                  < tolerance);

      // One angle at 90 degrees, another far from it
      UnitCell<3> plane;
      orthogonalPlaneCell(plane);
      TEST_ASSERT(compareCpu<3>(plane, IntVec<3>(8), kSq) < tolerance);

      // Rhombic cell at 120 degrees, 64x8 mesh (see testTransition)
      double ph120[2] = {1.0, 2.0*Constants::Pi/3.0};
      UnitCell<2> hexagonal;
      makeCell<2>(hexagonal, UnitCell<2>::Rhombic, ph120, 2);
      IntVec<2> da;
      da[0] = 64; da[1] = 8;
      TEST_ASSERT(compareCpu<2>(hexagonal, da, kSq) < tolerance);

      // Minimum image far from orthogonal geometries
      double po[3] = {1.9, 2.3, 2.2};
      UnitCell<2> oblique;
      makeCell<2>(oblique, UnitCell<2>::Oblique, po, 3);
      IntVec<2> d2;
      d2[0] = 8; d2[1] = 6;
      TEST_ASSERT(compareCpu<2>(oblique, d2, kSq) < tolerance);
      Mesh<2> mesh(d2);
      MeshIterator<2> iter(d2);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         double ref = oblique.ksq(shiftToMinimum(iter.position(), d2,
                                                 oblique));
         TEST_ASSERT(std::abs(kSq[iter.rank()] - ref) < tolerance*(1.0 + ref));
      }
   }

};

TEST_BEGIN(GpuNyquistWeightTest)
TEST_ADD(GpuNyquistWeightTest, testGpuSymmetricCells)
TEST_ADD(GpuNyquistWeightTest, testGpuTransition)
TEST_END(GpuNyquistWeightTest)

#endif
