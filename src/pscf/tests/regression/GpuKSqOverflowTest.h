#ifndef PSCF_GPU_KSQ_OVERFLOW_TEST_H
#define PSCF_GPU_KSQ_OVERFLOW_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/field/cuda/WaveList.h>
#include <pscf/backend/cuda/cudaTypes.h>
#include <pscf/backend/cuda/HostRecvArray.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

#include <util/containers/FSArray.h>

#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* GPU |k|^2 and its derivatives for wave indices whose products exceed
* the range of int (|k| > 46340, 1D meshes of 92682 or more points).
*/
class GpuKSqOverflowTest : public UnitTest
{

private:

   #ifdef SINGLE_PRECISION
   constexpr static double tolerance_ = 1E-5;
   #else
   constexpr static double tolerance_ = 1E-12;
   #endif

public:

   void setUp()
   {}

   void tearDown()
   {}

   /*
   * Relative difference of a and b.
   */
   double relDiff(double a, double b)
   {
      double scale = std::abs(a) + std::abs(b);
      if (scale == 0.0) return 0.0;
      return std::abs(a - b)/scale;
   }

   /*
   * Make a lamellar cell of length 1.6*n/1000 and a mesh of n points.
   */
   void makeLamellar(int n, UnitCell<1>& cell, Mesh<1>& mesh)
   {
      FSArray<double, 6> p;
      p.append(1.6*n/1000.0);
      cell.set(UnitCell<1>::Lamellar, p);
      IntVec<1> d;
      d[0] = n;
      mesh.setDimensions(d);
   }

   /*
   * Check kSq of a 1D mesh of n points against the exact values, as 
   * given by the minimum-image kernel (first call after allocate) and 
   * by the kSq kernel (after a change of the cell).
   */
   void checkKSq(int n)
   {
      UnitCell<1> cell;
      Mesh<1> mesh;
      makeLamellar(n, cell, mesh);
      int kSize = n/2 + 1;

      WaveList<1,CUT> waveList(true);
      waveList.allocate(mesh, cell);

      // kSq from the minimum-image kernel
      waveList.computeKSq();
      HostRecvArray<cudaReal,CUT> kSq_h(waveList.kSq());
      TEST_ASSERT(kSq_h.capacity() == kSize);
      IntVec<1> k;
      double maxError = 0.0;
      for (int i = 0; i < kSize; ++i) {
         k[0] = i;
         double error = relDiff(kSq_h[i], cell.ksq(k));
         if (error > maxError) maxError = error;
      }
      TEST_ASSERT(maxError < tolerance_);

      // kSq from the kSq kernel, after a change of the cell
      FSArray<double, 6> p = cell.parameters();
      p[0] *= 1.1;
      cell.setParameters(p);
      waveList.clearUnitCellData();
      waveList.computeKSq();
      kSq_h = waveList.kSq();
      maxError = 0.0;
      for (int i = 0; i < kSize; ++i) {
         k[0] = i;
         double error = relDiff(kSq_h[i], cell.ksq(k));
         if (error > maxError) maxError = error;
      }
      TEST_ASSERT(maxError < tolerance_);
   }

   /*
   * Check dKSq of a 1D mesh of n points against the exact values 
   * -2|k|^2/L (times 2 for waves with an implicit inverse).
   */
   void checkDKSq(int n)
   {
      UnitCell<1> cell;
      Mesh<1> mesh;
      makeLamellar(n, cell, mesh);
      double L = cell.parameter(0);
      int kSize = n/2 + 1;

      WaveList<1,CUT> waveList(true);
      waveList.allocate(mesh, cell);
      waveList.computedKSq();
      HostRecvArray<cudaReal,CUT> dKSq_h(waveList.dKSq(0));
      TEST_ASSERT(dKSq_h.capacity() == kSize);
      DArray<bool> const & implicit = waveList.implicitInverse();
      IntVec<1> k;
      double maxError = 0.0;
      for (int i = 0; i < kSize; ++i) {
         k[0] = i;
         double exact = -2.0*cell.ksq(k)/L;
         if (implicit[i]) exact *= 2.0;
         double error = relDiff(dKSq_h[i], exact);
         if (error > maxError) maxError = error;
      }
      TEST_ASSERT(maxError < tolerance_);
   }

   void testKSqBelowThreshold()
   {
      printMethod(TEST_FUNC);
      checkKSq(92680);
   }

   void testKSqAboveThreshold()
   {
      printMethod(TEST_FUNC);
      checkKSq(92682);
      checkKSq(150000);
   }

   void testDKSqBelowThreshold()
   {
      printMethod(TEST_FUNC);
      checkDKSq(92680);
   }

   void testDKSqAboveThreshold()
   {
      printMethod(TEST_FUNC);
      checkDKSq(92682);
      checkDKSq(150000);
   }

};

TEST_BEGIN(GpuKSqOverflowTest)
TEST_ADD(GpuKSqOverflowTest, testKSqBelowThreshold)
TEST_ADD(GpuKSqOverflowTest, testKSqAboveThreshold)
TEST_ADD(GpuKSqOverflowTest, testDKSqBelowThreshold)
TEST_ADD(GpuKSqOverflowTest, testDKSqAboveThreshold)
TEST_END(GpuKSqOverflowTest)

#endif
