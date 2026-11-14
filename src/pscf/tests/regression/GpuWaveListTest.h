#ifndef PSCF_GPU_WAVE_LIST_TEST_H
#define PSCF_GPU_WAVE_LIST_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/field/cuda/WaveList.h>
#include <pscf/field/cpp/WaveList.h>
#include <pscf/backend/cuda/cudaTypes.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/crystal/hasVariableAngle.h>

#include <pscf/backend/cuda/HostRecvArray.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/math/Constants.h>
#include <util/containers/FSArray.h>
#include <util/misc/Exception.h>

using namespace Util;
using namespace Pscf;

/*
* Regression tests of the GPU minimum images (WaveList<D,CUT>).
*/
class GpuWaveListTest : public UnitTest
{

private:

   // Error tolerance for array equality
   #ifdef SINGLE_PRECISION
   constexpr static cudaReal tolerance_ = 1E-5;
   #else
   #ifdef DOUBLE_PRECISION
   constexpr static cudaReal tolerance_ = 1E-10;
   #endif
   #endif

public:

   void setUp()
   {}

   void tearDown()
   {}

   /*
   * Compare the device minimum images and kSq with shiftToMinimum and
   * UnitCell::ksq, for all wavevectors of a mesh, except those whose
   * kSq is averaged over the signs of their Nyquist indices (lattice
   * systems with a variable angle, see usesNyquistAverage), which are
   * compared with the image and kSq of the CPU WaveList.
   */
   template <int D>
   void checkMinimumImages(UnitCell<D> const & cell, IntVec<D> const & d)
   {
      Mesh<D> mesh(d);
      IntVec<D> kDims = d;
      kDims[D-1] = d[D-1]/2 + 1;
      int kSize = 1;
      for (int i = 0; i < D; ++i) kSize *= kDims[i];

      WaveList<D,CUT> wavelist;
      wavelist.allocate(mesh, cell);
      HostRecvArray<int,CUT> minImages_h;
      HostRecvArray<cudaReal,CUT> ksq_h;
      wavelist.computeMinimumImages();
      minImages_h = wavelist.minImages_d();
      ksq_h = wavelist.kSq();
      WaveList<D,CPT> cpuList;
      cpuList.allocate(mesh, cell);
      cpuList.computeMinimumImages();

      bool variableAngle = hasVariableAngle<D>(cell.lattice());
      IntVec<D> temp, vec;
      double ksq;
      MeshIterator<D> iter;
      iter.setDimensions(kDims);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         if (variableAngle && usesNyquistAverage(temp, d, cell)) {
            vec = cpuList.minImages()[iter.rank()];
            ksq = cpuList.kSq()[iter.rank()];
         } else {
            vec = shiftToMinimum(temp, d, cell);
            ksq = cell.ksq(vec);
         }
         for (int i = 0; i < D; ++i) {
            TEST_ASSERT(vec[i] == minImages_h[iter.rank() + i*kSize]);
         }
         TEST_ASSERT(fabs(ksq - ksq_h[iter.rank()]) < tolerance_);
      }
   }

   /*
   * Minimum images on every 2D lattice system, including square and
   * hexagonal cells with many images of equal |k|, on even and odd
   * meshes. The 32 images of each wavevector are reduced within one
   * warp: the device images must equal those of shiftToMinimum
   * (except for the Nyquist wavevectors of the rhombic and oblique
   * cells, whose images are computed on the host).
   */
   void testComputeMinimumImages2DLattices()
   {
      printMethod(TEST_FUNC);

      UnitCell<2> cells[5];
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      cells[0].set(UnitCell<2>::Square, parameters);
      cells[1].set(UnitCell<2>::Hexagonal, parameters);
      parameters.append(1.2);
      cells[2].set(UnitCell<2>::Rhombic, parameters);
      parameters[1] = 3.1;
      cells[3].set(UnitCell<2>::Rectangular, parameters);
      parameters[0] = 1.8;
      parameters[1] = 1.3;
      parameters.append(1.7);
      cells[4].set(UnitCell<2>::Oblique, parameters);

      int meshes[3][2] = {{64, 64}, {40, 48}, {45, 63}};
      for (int ic = 0; ic < 5; ++ic) {
         for (int im = 0; im < 3; ++im) {
            IntVec<2> d;
            d[0] = meshes[im][0];
            d[1] = meshes[im][1];
            checkMinimumImages<2>(cells[ic], d);
         }
      }
   }

   /*
   * Skewed cells with strongly anisotropic meshes, where the fixed
   * image search windows alone miss minimum images: the device images
   * and kSq must equal those of the host (checkMinimumImages).
   */
   void testComputeMinimumImagesSkewed()
   {
      printMethod(TEST_FUNC);

      double gammas[2] = {2.4, 2.8};
      int meshes2[2][2] = {{4, 64}, {8, 64}};
      for (int ig = 0; ig < 2; ++ig) {
         FSArray<double, 6> parameters;
         parameters.append(2.0);
         parameters.append(2.0);
         parameters.append(gammas[ig]);
         UnitCell<2> cell;
         cell.set(UnitCell<2>::Oblique, parameters);
         for (int im = 0; im < 2; ++im) {
            IntVec<2> d;
            d[0] = meshes2[im][0];
            d[1] = meshes2[im][1];
            checkMinimumImages<2>(cell, d);
         }
      }

      FSArray<double, 6> parameters;
      parameters.append(2.0);
      parameters.append(2.0);
      parameters.append(2.0);
      parameters.append(2.4);
      UnitCell<3> cell;
      cell.set(UnitCell<3>::Monoclinic, parameters);
      IntVec<3> d;
      d[0] = 4;
      d[1] = 2;
      d[2] = 32;
      checkMinimumImages<3>(cell, d);
   }

   /*
   * A 3D mesh with more than 2^24 wavevectors, for which the number
   * of threads used for the minimum images (128 per wavevector)
   * exceeds the range of int. Every wavevector must still get its
   * minimum image and kSq.
   *
   * The cell is orthorhombic, so that the minimum image is the one of
   * smallest magnitude in each direction separately. This is checked
   * for all wavevectors. A sample of images is also compared with
   * shiftToMinimum.
   */
   void testComputeMinimumImagesLarge()
   {
      printMethod(TEST_FUNC);

      IntVec<3> d;
      d[0] = 256;
      d[1] = 258;
      d[2] = 510;
      Mesh<3> mesh(d);
      IntVec<3> kDims = d;
      kDims[2] = d[2]/2 + 1;
      int kSize = kDims[0]*kDims[1]*kDims[2];
      TEST_ASSERT(kSize > (1 << 24));

      FSArray<double, 6> parameters;
      parameters.append(2.0);
      parameters.append(2.5);
      parameters.append(3.0);
      UnitCell<3> cell;
      cell.set(UnitCell<3>::Orthorhombic, parameters);

      WaveList<3,CUT> wavelist;
      wavelist.allocate(mesh, cell);
      HostRecvArray<int,CUT> minImages_h;
      HostRecvArray<cudaReal,CUT> ksq_h;
      cudaError_t error = cudaSuccess;
      bool thrown = false;
      try {
         wavelist.computeMinimumImages();
         error = cudaGetLastError(); // also clears a launch error
         minImages_h = wavelist.minImages_d();
         ksq_h = wavelist.kSq();
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
      TEST_ASSERT(error == cudaSuccess);

      // All wavevectors: each component of the image is congruent with
      // the k-grid position and has the smallest magnitude; kSq is the
      // square magnitude of that image
      double b[3];
      for (int i = 0; i < 3; ++i) {
         b[i] = 2.0*Constants::Pi/parameters[i];
      }
      IntVec<3> pos;
      int rank = 0;
      for (pos[0] = 0; pos[0] < kDims[0]; ++pos[0]) {
         for (pos[1] = 0; pos[1] < kDims[1]; ++pos[1]) {
            for (pos[2] = 0; pos[2] < kDims[2]; ++pos[2]) {
               double kSq = 0.0;
               for (int i = 0; i < 3; ++i) {
                  int m = pos[i];
                  if (2*m > d[i]) m -= d[i];
                  int img = minImages_h[rank + i*kSize];
                  TEST_ASSERT((img - pos[i]) % d[i] == 0);
                  TEST_ASSERT(abs(img) == abs(m));
                  kSq += b[i]*b[i]*double(m*m);
               }
               TEST_ASSERT(fabs(ksq_h[rank] - kSq) < tolerance_*(1.0 + kSq));
               ++rank;
            }
         }
      }
      TEST_ASSERT(rank == kSize);

      // A sample of wavevectors, and the last ones: images equal those
      // of shiftToMinimum (also for images of equal magnitude)
      IntVec<3> vec;
      rank = 0;
      while (rank < kSize) {
         pos[2] = rank % kDims[2];
         pos[1] = (rank / kDims[2]) % kDims[1];
         pos[0] = rank / (kDims[1]*kDims[2]);
         vec = shiftToMinimum(pos, d, cell);
         for (int i = 0; i < 3; ++i) {
            TEST_ASSERT(vec[i] == minImages_h[rank + i*kSize]);
         }
         if (rank < kSize - 2048) {
            rank += 4999;
            if (rank > kSize - 2048) rank = kSize - 2048;
         } else {
            ++rank;
         }
      }
   }

};

TEST_BEGIN(GpuWaveListTest)
TEST_ADD(GpuWaveListTest, testComputeMinimumImages2DLattices)
TEST_ADD(GpuWaveListTest, testComputeMinimumImagesSkewed)
TEST_ADD(GpuWaveListTest, testComputeMinimumImagesLarge)
TEST_END(GpuWaveListTest)

#endif
