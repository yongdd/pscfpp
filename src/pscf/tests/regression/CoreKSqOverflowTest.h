#ifndef PSCF_CORE_KSQ_OVERFLOW_TEST_H
#define PSCF_CORE_KSQ_OVERFLOW_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/crystal/UnitCell.h>
#include <pscf/field/cpp/WaveList.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/containers/FSArray.h>
#include <util/math/Constants.h>

#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Derivatives of |k|^2 for wave indices whose products exceed the
* range of int (|k| > 46340, e.g. on 1D meshes of 92682 or more points).
*/
class CoreKSqOverflowTest : public UnitTest
{

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
   * UnitCell<1>::dksq for indices around and above 46341, against the
   * exact d|k|^2/dL = -2|k|^2/L of a lamellar cell.
   */
   void testDksq1D()
   {
      printMethod(TEST_FUNC);

      double L = 1.6;
      UnitCell<1> cell;
      FSArray<double, 6> p;
      p.append(L);
      cell.set(UnitCell<1>::Lamellar, p);

      int values[5] = {1, 46340, 46341, 50000, 150000};
      IntVec<1> k;
      for (int i = 0; i < 5; ++i) {
         k[0] = values[i];
         double exact = -2.0*cell.ksq(k)/L;
         TEST_ASSERT(exact < 0.0);
         TEST_ASSERT(relDiff(cell.dksq(k, 0), exact) < 1.0E-12);
         k[0] = -values[i];
         TEST_ASSERT(relDiff(cell.dksq(k, 0), exact) < 1.0E-12);
      }
   }

   /*
   * UnitCell<3>::dksq with large indices: cubic cell against the exact
   * -2|k|^2/a, triclinic cell against central finite differences.
   */
   void testDksq3D()
   {
      printMethod(TEST_FUNC);

      IntVec<3> k;
      k[0] = 50000;
      k[1] = -47000;
      k[2] = 3;

      // Cubic
      {
         double a = 2.3;
         UnitCell<3> cell;
         FSArray<double, 6> p;
         p.append(a);
         cell.set(UnitCell<3>::Cubic, p);
         double exact = -2.0*cell.ksq(k)/a;
         TEST_ASSERT(relDiff(cell.dksq(k, 0), exact) < 1.0E-12);
      }

      // Triclinic
      {
         FSArray<double, 6> p;
         p.append(2.0);
         p.append(2.3);
         p.append(2.6);
         p.append(1.4);
         p.append(0.2);
         p.append(0.1);
         UnitCell<3> cell;
         cell.set(UnitCell<3>::Triclinic, p);
         UnitCell<3> cp, cm;
         cp = cell;
         cm = cell;
         double h = 1.0E-5;
         for (int n = 0; n < 6; ++n) {
            FSArray<double, 6> pp = p;
            FSArray<double, 6> pm = p;
            pp[n] += h;
            pm[n] -= h;
            cp.setParameters(pp);
            cm.setParameters(pm);
            double fd = (cp.ksq(k) - cm.ksq(k))/(2.0*h);
            double scale = cell.ksq(k);
            TEST_ASSERT(std::abs(cell.dksq(k, n) - fd) < 1.0E-6*scale);
         }
      }
   }

   /*
   * WaveList<1,CPT>::dKSq on a 1D mesh of 100000 points (indices up to
   * 50000), against the exact -2|k|^2/L (times 2 for implicit inverses).
   */
   void testWaveListDKSq1D()
   {
      printMethod(TEST_FUNC);

      int n = 100000;
      double L = 1.6*n/1000.0;
      UnitCell<1> cell;
      FSArray<double, 6> p;
      p.append(L);
      cell.set(UnitCell<1>::Lamellar, p);
      IntVec<1> d;
      d[0] = n;
      Mesh<1> mesh(d);

      WaveList<1,CPT> waveList(true);
      waveList.allocate(mesh, cell);
      waveList.computeMinimumImages();
      waveList.computeKSq();
      waveList.computedKSq();

      RField<1,CPT> const & kSq = waveList.kSq();
      RField<1,CPT> const & dKSq = waveList.dKSq(0);
      DArray<bool> const & implicit = waveList.implicitInverse();
      int kSize = n/2 + 1;
      TEST_ASSERT(kSq.capacity() == kSize);
      double maxError = 0.0;
      IntVec<1> k;
      for (int i = 0; i < kSize; ++i) {
         k[0] = i;
         double kSqExact = cell.ksq(k);
         TEST_ASSERT(relDiff(kSq[i], kSqExact) < 1.0E-12);
         double exact = -2.0*kSqExact/L;
         if (implicit[i]) exact *= 2.0;
         double error = relDiff(dKSq[i], exact);
         if (error > maxError) maxError = error;
      }
      TEST_ASSERT(maxError < 1.0E-12);
   }

};

TEST_BEGIN(CoreKSqOverflowTest)
TEST_ADD(CoreKSqOverflowTest, testDksq1D)
TEST_ADD(CoreKSqOverflowTest, testDksq3D)
TEST_ADD(CoreKSqOverflowTest, testWaveListDKSq1D)
TEST_END(CoreKSqOverflowTest)

#endif
