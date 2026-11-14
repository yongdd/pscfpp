#ifndef PSCF_CORE_NYQUIST_TEST_H
#define PSCF_CORE_NYQUIST_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/field/cpp/WaveList.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/crystal/Basis.h>
#include <pscf/crystal/SpaceGroup.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/containers/FSArray.h>
#include <util/containers/DArray.h>
#include <util/math/Constants.h>

#include <cmath>
#include <string>

using namespace Util;
using namespace Pscf;

/*
* Wavevectors with Nyquist indices (even meshes) in the CPU WaveList.
*
* For lattice systems with a variable angle, kSq of such a wavevector
* is |k|^2 averaged over the signs of its Nyquist indices near 
* orthogonal geometries, so that kSq and the free energy are smooth 
* where an angle passes through 90 degrees, and the minimum-image |k|^2
* far from them (see ksqNyquist); dKSq is the derivative of kSq. Other
* lattice systems keep the minimum-image |k|^2.
*/
class CoreNyquistTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   /*
   * Make a unit cell from a lattice system and parameters.
   */
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
   * Is wave index u on a Nyquist plane of a mesh dimension n?
   */
   bool isNyquist(int u, int n)
   {
      if (n % 2 != 0) return false;
      int r = u % n;
      if (r < 0) r += n;
      return (2*r == n);
   }

   /*
   * Compare dKSq of a complex-field WaveList with central finite
   * differences of kSq, for all wavevectors and all cell parameters.
   * Waves whose stored image changes between the displaced cells (an
   * image switch, where kSq has a kink) are skipped if skipSwitch is
   * true, and are counted otherwise. Returns the maximum error relative
   * to max(1, |kSq|).
   */
   template <int D>
   double finiteDifferenceError(UnitCell<D> const & cell,
                                IntVec<D> const & d,
                                bool skipSwitch = false)
   {
      Mesh<D> mesh(d);
      WaveList<D,CPT> waveList(false);
      waveList.allocate(mesh, cell);
      waveList.computeMinimumImages();
      waveList.computedKSq();
      int kSize = mesh.size();

      FSArray<double, 6> p = cell.parameters();
      int nParam = cell.nParameter();
      double h = 1.0E-6;
      double maxError = 0.0;
      for (int n = 0; n < nParam; ++n) {
         UnitCell<D> cp, cm;
         cp = cell;
         cm = cell;
         FSArray<double, 6> pp = p;
         FSArray<double, 6> pm = p;
         pp[n] += h;
         pm[n] -= h;
         cp.setParameters(pp);
         cm.setParameters(pm);
         WaveList<D,CPT> wp(false), wm(false);
         wp.allocate(mesh, cp);
         wm.allocate(mesh, cm);
         wp.computeMinimumImages();
         wm.computeMinimumImages();
         RField<D,CPT> const & dKSq = waveList.dKSq(n);
         for (int i = 0; i < kSize; ++i) {
            if (skipSwitch) {
               if (!(wp.minImages()[i] == wm.minImages()[i])) continue;
            }
            double fd = (wp.kSq()[i] - wm.kSq()[i])/(2.0*h);
            double scale = std::max(1.0, std::abs(waveList.kSq()[i]));
            double error = std::abs(dKSq[i] - fd)/scale;
            if (error > maxError) maxError = error;
         }
      }
      return maxError;
   }

   /*
   * Lattice systems without a variable angle: kSq and dKSq are those
   * of the minimum image, exactly (bitwise), also on Nyquist planes.
   */
   void testFixedAngleUnchanged()
   {
      printMethod(TEST_FUNC);

      double p3[3] = {1.9, 2.1, 2.3};
      UnitCell<3> cells3[5];
      makeCell<3>(cells3[0], UnitCell<3>::Cubic, p3, 1);
      makeCell<3>(cells3[1], UnitCell<3>::Tetragonal, p3, 2);
      makeCell<3>(cells3[2], UnitCell<3>::Orthorhombic, p3, 3);
      makeCell<3>(cells3[3], UnitCell<3>::Hexagonal, p3, 2);
      makeCell<3>(cells3[4], UnitCell<3>::Hexagonal, p3, 2);
      IntVec<3> d3[5];
      d3[0] = IntVec<3>(8);
      d3[1][0] = 8; d3[1][1] = 8; d3[1][2] = 6;
      d3[2][0] = 8; d3[2][1] = 6; d3[2][2] = 10;
      d3[3][0] = 6; d3[3][1] = 6; d3[3][2] = 8;
      d3[4][0] = 12; d3[4][1] = 12; d3[4][2] = 4;
      for (int ic = 0; ic < 5; ++ic) {
         checkMinimumImageValues<3>(cells3[ic], d3[ic]);
      }

      double p2[2] = {1.9, 2.3};
      UnitCell<2> cells2[3];
      makeCell<2>(cells2[0], UnitCell<2>::Square, p2, 1);
      makeCell<2>(cells2[1], UnitCell<2>::Rectangular, p2, 2);
      makeCell<2>(cells2[2], UnitCell<2>::Hexagonal, p2, 1);
      IntVec<2> d2[3];
      d2[0] = IntVec<2>(8);
      d2[1][0] = 8; d2[1][1] = 6;
      d2[2] = IntVec<2>(12);
      for (int ic = 0; ic < 3; ++ic) {
         checkMinimumImageValues<2>(cells2[ic], d2[ic]);
      }
   }

   /*
   * Exact comparison of a real-field WaveList with the minimum images.
   */
   template <int D>
   void checkMinimumImageValues(UnitCell<D> const & cell,
                                IntVec<D> const & d)
   {
      Mesh<D> mesh(d);
      WaveList<D,CPT> waveList(true);
      waveList.allocate(mesh, cell);
      waveList.computeMinimumImages();
      waveList.computedKSq();
      DArray<bool> const & implicit = waveList.implicitInverse();
      MeshIterator<D> iter(waveList.kMeshDimensions());
      for (iter.begin(); !iter.atEnd(); ++iter) {
         int rank = iter.rank();
         IntVec<D> v = shiftToMinimum(iter.position(), d, cell);
         TEST_ASSERT(v == waveList.minImages()[rank]);
         TEST_ASSERT(waveList.kSq()[rank] == cell.ksq(v));
         for (int n = 0; n < cell.nParameter(); ++n) {
            double dksq = cell.dksq(v, n);
            if (implicit[rank]) dksq *= 2.0;
            TEST_ASSERT(waveList.dKSq(n)[rank] == dksq);
         }
      }
   }

   /*
   * Cells with variable angles at 90 degrees, on even meshes: dKSq
   * equals the central finite difference of kSq for every wavevector,
   * as kSq is smooth through the right angle. (The minimum-image |k|^2
   * of a Nyquist wavevector has a kink there.)
   */
   void testSmoothAtRightAngles()
   {
      printMethod(TEST_FUNC);
      double halfPi = 0.5*Constants::Pi;
      double tolerance = 1.0E-6;

      // Monoclinic (a, b, c, beta)
      {
         double p[4] = {1.9, 2.1, 2.3, halfPi};
         UnitCell<3> cell;
         makeCell<3>(cell, UnitCell<3>::Monoclinic, p, 4);
         IntVec<3> d;
         d[0] = 8; d[1] = 6; d[2] = 10;
         TEST_ASSERT(finiteDifferenceError<3>(cell, d) < tolerance);
      }

      // Rhombohedral (a, alpha)
      {
         double p[2] = {2.0, halfPi};
         UnitCell<3> cell;
         makeCell<3>(cell, UnitCell<3>::Rhombohedral, p, 2);
         IntVec<3> d(8);
         TEST_ASSERT(finiteDifferenceError<3>(cell, d) < tolerance);
      }

      // Triclinic (a, b, c, phi, theta, gamma), with c tilted by theta
      // from the normal of the (a, b) plane (orthogonal at theta = 0)
      {
         double p[6] = {1.9, 2.1, 2.3, 0.4, 0.0, halfPi};
         UnitCell<3> cell;
         makeCell<3>(cell, UnitCell<3>::Triclinic, p, 6);
         IntVec<3> d;
         d[0] = 6; d[1] = 8; d[2] = 10;
         TEST_ASSERT(finiteDifferenceError<3>(cell, d) < tolerance);
      }

      // Rhombic (a, gamma) and oblique (a, b, gamma)
      {
         double p[2] = {2.0, halfPi};
         UnitCell<2> cell;
         makeCell<2>(cell, UnitCell<2>::Rhombic, p, 2);
         IntVec<2> d(8);
         TEST_ASSERT(finiteDifferenceError<2>(cell, d) < tolerance);
      }
      {
         double p[3] = {1.9, 2.3, halfPi};
         UnitCell<2> cell;
         makeCell<2>(cell, UnitCell<2>::Oblique, p, 3);
         IntVec<2> d;
         d[0] = 8; d[1] = 6;
         TEST_ASSERT(finiteDifferenceError<2>(cell, d) < tolerance);
      }
   }

   /*
   * Skewed cells: dKSq equals the finite difference of kSq, except at
   * image switches, and kSq is the minimum over images of |k|^2 
   * averaged over the sign reversals of the components of Nyquist
   * indices that are decoupled from the other nonzero indices (see
   * ksqNyquist), e.g. the index normal to the others in the monoclinic
   * cell, otherwise the minimum-image |k|^2 (checked by a search over 
   * many images).
   */
   void testSkewedCells()
   {
      printMethod(TEST_FUNC);

      double p6[6] = {1.9, 2.1, 2.3, 1.3, 0.5, 1.2};
      UnitCell<3> tric;
      makeCell<3>(tric, UnitCell<3>::Triclinic, p6, 6);
      IntVec<3> d;
      d[0] = 6; d[1] = 8; d[2] = 10;
      TEST_ASSERT(finiteDifferenceError<3>(tric, d, true) < 1.0E-6);
      checkSearch<3>(tric, d);

      double p4[4] = {1.9, 2.1, 2.3, 2.0};
      UnitCell<3> mono;
      makeCell<3>(mono, UnitCell<3>::Monoclinic, p4, 4);
      TEST_ASSERT(finiteDifferenceError<3>(mono, d, true) < 1.0E-6);
      checkSearch<3>(mono, d);

      double p2[2] = {2.0, 1.05};
      UnitCell<3> rhomb;
      makeCell<3>(rhomb, UnitCell<3>::Rhombohedral, p2, 2);
      IntVec<3> d8(8);
      TEST_ASSERT(finiteDifferenceError<3>(rhomb, d8, true) < 1.0E-6);
      checkSearch<3>(rhomb, d8);

      double p3[3] = {1.9, 2.3, 2.2};
      UnitCell<2> obl;
      makeCell<2>(obl, UnitCell<2>::Oblique, p3, 3);
      IntVec<2> d2;
      d2[0] = 8; d2[1] = 6;
      TEST_ASSERT(finiteDifferenceError<2>(obl, d2, true) < 1.0E-6);
      checkSearch<2>(obl, d2);
   }

   /*
   * |k|^2 of u averaged over the reversals of the signs of the indices
   * of the components comps[c] (bit masks of directions; all indices of
   * a component are reversed together).
   */
   template <int D>
   double groupAveragedKSq(IntVec<D> const & u, int const * comps, 
                           int nComp, UnitCell<D> const & cell)
   {
      double sum = 0.0;
      int nImage = 1 << nComp;
      for (int s = 0; s < nImage; ++s) {
         IntVec<D> w = u;
         for (int c = 0; c < nComp; ++c) {
            if (!(s & (1 << c))) continue;
            for (int i = 0; i < D; ++i) {
               if (comps[c] & (1 << i)) w[i] = -w[i];
            }
         }
         sum += cell.ksq(w);
      }
      return sum/double(nImage);
   }

   /*
   * Components averaged by ksqNyquist, for a cell whose decoupling 
   * weights (nyquistDecouplingWeight) are all 0 or 1: the connected
   * components of the directions of nonzero indices of v, with the 
   * pairs of weight 0 as edges, whose indices are all Nyquist indices,
   * except a component of all nonzero indices. Returns their number and
   * sets their bit masks in comps.
   */
   template <int D>
   int averagedComponents(IntVec<D> const & v, IntVec<D> const & d,
                          UnitCell<D> const & cell, int* comps)
   {
      int comp[D];
      for (int i = 0; i < D; ++i) comp[i] = (1 << i);
      for (int i = 0; i < D; ++i) {
         if (v[i] % d[i] == 0) continue;
         for (int j = i + 1; j < D; ++j) {
            if (v[j] % d[j] == 0) continue;
            double w = nyquistDecouplingWeight<D>(i, j, cell);
            TEST_ASSERT(w == 0.0 || w == 1.0);
            if (w == 0.0) {
               int merged = comp[i] | comp[j];
               for (int k = 0; k < D; ++k) {
                  if (merged & (1 << k)) comp[k] = merged;
               }
            }
         }
      }
      int nonzero = 0;
      for (int i = 0; i < D; ++i) {
         if (v[i] % d[i] != 0) nonzero |= (1 << i);
      }
      int nComp = 0;
      int listed = 0;
      for (int i = 0; i < D; ++i) {
         if (v[i] % d[i] == 0 || (listed & comp[i])) continue;
         bool allNyquist = (comp[i] != nonzero);
         for (int k = 0; k < D; ++k) {
            if ((comp[i] & (1 << k)) && !isNyquist(v[k], d[k])) {
               allNyquist = false;
            }
         }
         listed |= comp[i];
         if (allNyquist) {
            comps[nComp] = comp[i];
            ++nComp;
         }
      }
      return nComp;
   }

   /*
   * Search over images v + s*d with |s[i]| <= 3, for a cell whose
   * decoupling weights are 0 or 1: the WaveList kSq is the minimum of
   * |k|^2 averaged over the sign reversals of the averaged components
   * (the minimum |k|^2 for waves without a Nyquist index), and equal 
   * to the value of the stored image.
   */
   template <int D>
   void checkSearch(UnitCell<D> const & cell, IntVec<D> const & d)
   {
      Mesh<D> mesh(d);
      WaveList<D,CPT> waveList(false);
      waveList.allocate(mesh, cell);
      waveList.computeMinimumImages();
      MeshIterator<D> iter(d);
      int nNyquist = 0;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         IntVec<D> v = iter.position();
         bool hasNyquist = false;
         for (int i = 0; i < D; ++i) {
            if (isNyquist(v[i], d[i])) hasNyquist = true;
         }
         int comps[D];
         int nComp = 0;
         if (hasNyquist) {
            ++nNyquist;
            nComp = averagedComponents<D>(v, d, cell, comps);
         }
         double best = -1.0;
         IntVec<D> s(-3);
         bool done = false;
         while (!done) {
            IntVec<D> u;
            for (int i = 0; i < D; ++i) u[i] = v[i] + s[i]*d[i];
            double g = groupAveragedKSq<D>(u, comps, nComp, cell);
            if (best < 0.0 || g < best) best = g;
            // Next shift vector
            int i = 0;
            while (i < D) {
               if (s[i] < 3) {
                  ++s[i];
                  break;
               }
               s[i] = -3;
               ++i;
            }
            if (i == D) done = true;
         }
         int rank = iter.rank();
         double kSq = waveList.kSq()[rank];
         TEST_ASSERT(std::abs(kSq - best) < 1.0E-10*(1.0 + best));
         IntVec<D> const & image = waveList.minImages()[rank];
         double g = groupAveragedKSq<D>(image, comps, nComp, cell);
         TEST_ASSERT(std::abs(kSq - g) < 1.0E-10*(1.0 + best));
         for (int i = 0; i < D; ++i) {
            TEST_ASSERT((image[i] - v[i]) % d[i] == 0);
         }
      }
      TEST_ASSERT(nNyquist > 0);
   }

   /*
   * The kSq of all waves in each star of a symmetry-adapted basis are
   * equal, for groups of variable-angle lattices on even meshes.
   */
   void testStarsHaveEqualKSq()
   {
      printMethod(TEST_FUNC);

      double pr[2] = {2.0, 1.2};
      UnitCell<3> rhomb;
      makeCell<3>(rhomb, UnitCell<3>::Rhombohedral, pr, 2);
      IntVec<3> d8(8);
      checkStars<3>(rhomb, d8, "R_-3_m:R");

      double pm[4] = {1.9, 2.1, 2.3, 1.9};
      UnitCell<3> mono;
      makeCell<3>(mono, UnitCell<3>::Monoclinic, pm, 4);
      IntVec<3> d;
      d[0] = 8; d[1] = 6; d[2] = 10;
      checkStars<3>(mono, d, "P_1_2%c_1");

      double pc[2] = {2.0, 0.5*Constants::Pi};
      UnitCell<3> cub;
      makeCell<3>(cub, UnitCell<3>::Rhombohedral, pc, 2);
      checkStars<3>(cub, d8, "R_-3_m:R");
   }

   template <int D>
   void checkStars(UnitCell<D> const & cell, IntVec<D> const & d,
                   std::string groupName)
   {
      Mesh<D> mesh(d);
      SpaceGroup<D> group;
      readGroup(groupName, group);
      Basis<D> basis;
      basis.makeBasis(mesh, cell, group);
      WaveList<D,CPT> waveList(false);
      waveList.allocate(mesh, cell);
      waveList.computeMinimumImages();
      RField<D,CPT> const & kSq = waveList.kSq();
      for (int is = 0; is < basis.nStar(); ++is) {
         typename Basis<D>::Star const & star = basis.star(is);
         int r0 = mesh.rank(basis.wave(star.beginId).indicesStd);
         for (int iw = star.beginId; iw < star.endId; ++iw) {
            int r = mesh.rank(basis.wave(iw).indicesStd);
            TEST_ASSERT(std::abs(kSq[r] - kSq[r0])
                        < 1.0E-10*(1.0 + kSq[r0]));
         }
      }
   }

};

TEST_BEGIN(CoreNyquistTest)
TEST_ADD(CoreNyquistTest, testFixedAngleUnchanged)
TEST_ADD(CoreNyquistTest, testSmoothAtRightAngles)
TEST_ADD(CoreNyquistTest, testSkewedCells)
TEST_ADD(CoreNyquistTest, testStarsHaveEqualKSq)
TEST_END(CoreNyquistTest)

#endif
