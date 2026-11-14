#ifndef PSCF_CORE_NYQUIST_WEIGHT_TEST_H
#define PSCF_CORE_NYQUIST_WEIGHT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/field/cpp/WaveList.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/containers/FSArray.h>
#include <util/math/Constants.h>

#include <cmath>
#include <vector>

using namespace Util;
using namespace Pscf;

/*
* Average of kSq over the signs of Nyquist indices, CPU WaveList.
*
* For lattice systems with a variable angle, kSq of a wavevector with a
* Nyquist index (even meshes) was averaged over the signs of its Nyquist
* indices for all cell parameters. This average is invariant only under
* symmetry operations that are signed permutations of the basis. In a
* cell whose metric has other symmetry operations (a rhombic, oblique,
* monoclinic or triclinic cell with a 120 degree angle, a rhombohedral
* cell at 60 or 109.47 degrees), it gave different kSq to wavevectors
* related by these operations, so that a state with the symmetry of the
* cell was not stationary there. The average is now used only near
* orthogonal geometries, where the minimum image has a kink, and the
* minimum image (which has the symmetry of every lattice) far from them,
* with a smooth transition.
*/
class CoreNyquistWeightTest : public UnitTest
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
   * Reciprocal metric element b_i . b_j.
   */
   template <int D>
   double metric(UnitCell<D> const & cell, int i, int j)
   {
      double value = 0.0;
      for (int k = 0; k < D; ++k) {
         value += cell.kBasis(i)[k]*cell.kBasis(j)[k];
      }
      return value;
   }

   /*
   * Symmetry operations of the lattice that map the mesh onto itself:
   * integer matrices M with elements -1, 0 or 1 and |det M| = 1, with
   * M g M^T = g for the reciprocal metric g, acting on the indices of
   * wavevectors as v -> v M. Element (i, j) of operation k is stored in
   * ops[k][D*i + j]. Returns the number of operations that are not
   * signed permutations.
   */
   template <int D>
   int latticeOperations(UnitCell<D> const & cell, IntVec<D> const & d,
                         std::vector< std::vector<int> >& ops)
   {
      double g[3][3];
      double gMax = 0.0;
      for (int i = 0; i < D; ++i) {
         for (int j = 0; j < D; ++j) {
            g[i][j] = metric<D>(cell, i, j);
            gMax = std::max(gMax, std::abs(g[i][j]));
         }
      }
      ops.clear();
      int nOther = 0;
      int nMatrix = 1;
      for (int k = 0; k < D*D; ++k) nMatrix *= 3;
      std::vector<int> m(D*D);
      for (int code = 0; code < nMatrix; ++code) {
         int c = code;
         for (int k = 0; k < D*D; ++k) {
            m[k] = (c % 3) - 1;
            c /= 3;
         }

         // Determinant
         int det;
         if (D == 2) {
            det = m[0]*m[3] - m[1]*m[2];
         } else {
            det = m[0]*(m[4]*m[8] - m[5]*m[7])
                - m[1]*(m[3]*m[8] - m[5]*m[6])
                + m[2]*(m[3]*m[7] - m[4]*m[6]);
         }
         if (det != 1 && det != -1) continue;

         // Mesh: off-diagonal elements only between equal dimensions
         bool ok = true;
         for (int i = 0; i < D; ++i) {
            for (int j = 0; j < D; ++j) {
               if (i != j && m[D*i + j] != 0 && d[i] != d[j]) ok = false;
            }
         }
         if (!ok) continue;

         // Metric: M g M^T = g
         for (int i = 0; i < D && ok; ++i) {
            for (int j = 0; j < D && ok; ++j) {
               double s = 0.0;
               for (int k = 0; k < D; ++k) {
                  for (int l = 0; l < D; ++l) {
                     s += m[D*i + k]*g[k][l]*m[D*j + l];
                  }
               }
               if (std::abs(s - g[i][j]) > 1.0E-9*gMax) ok = false;
            }
         }
         if (!ok) continue;

         // Signed permutation?
         bool isSignedPermutation = true;
         for (int i = 0; i < D; ++i) {
            int nRow = 0;
            int nColumn = 0;
            for (int j = 0; j < D; ++j) {
               if (m[D*i + j] != 0) ++nRow;
               if (m[D*j + i] != 0) ++nColumn;
            }
            if (nRow != 1 || nColumn != 1) isSignedPermutation = false;
         }
         if (!isSignedPermutation) ++nOther;
         ops.push_back(m);
      }
      return nOther;
   }

   /*
   * kSq is equal for all wavevectors related by symmetry operations of
   * the lattice (complex-field WaveList, all wavevectors). Returns the
   * largest relative difference.
   */
   template <int D>
   double orbitDifference(UnitCell<D> const & cell, IntVec<D> const & d)
   {
      std::vector< std::vector<int> > ops;
      int nOther = latticeOperations<D>(cell, d, ops);
      TEST_ASSERT(nOther > 0);

      Mesh<D> mesh(d);
      WaveList<D,CPT> waveList(false);
      waveList.allocate(mesh, cell);
      waveList.computeMinimumImages();
      RField<D,CPT> const & kSq = waveList.kSq();

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
   * |k|^2 averaged over the signs of the Nyquist indices of u.
   */
   template <int D>
   double averagedKSq(IntVec<D> const & u, IntVec<D> const & d,
                      UnitCell<D> const & cell)
   {
      int nyquist[D];
      int m = 0;
      for (int i = 0; i < D; ++i) {
         if (isNyquist(u[i], d[i])) {
            nyquist[m] = i;
            ++m;
         }
      }
      double sum = 0.0;
      int nImage = 1 << m;
      for (int s = 0; s < nImage; ++s) {
         IntVec<D> w = u;
         for (int j = 0; j < m; ++j) {
            if (s & (1 << j)) w[nyquist[j]] = -w[nyquist[j]];
         }
         sum += cell.ksq(w);
      }
      return sum/double(nImage);
   }

   /*
   * Minimum over the images v + s*d, |s[i]| <= 3, of |k|^2 (average
   * is false) or of |k|^2 averaged over the signs of the Nyquist
   * indices (average is true).
   */
   template <int D>
   double searchMinimum(IntVec<D> const & v, IntVec<D> const & d,
                        UnitCell<D> const & cell, bool average)
   {
      double best = -1.0;
      IntVec<D> s(-3);
      bool done = false;
      while (!done) {
         IntVec<D> u;
         for (int i = 0; i < D; ++i) u[i] = v[i] + s[i]*d[i];
         double g = average ? averagedKSq<D>(u, d, cell) : cell.ksq(u);
         if (best < 0.0 || g < best) best = g;
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
      return best;
   }

   /*
   * kSq of every wavevector equals the minimum over images of |k|^2
   * (average false) or of its average over the signs of the Nyquist
   * indices (average true). Returns the largest relative error.
   */
   template <int D>
   double searchError(UnitCell<D> const & cell, IntVec<D> const & d,
                      bool average)
   {
      Mesh<D> mesh(d);
      WaveList<D,CPT> waveList(false);
      waveList.allocate(mesh, cell);
      waveList.computeMinimumImages();
      double maxError = 0.0;
      MeshIterator<D> iter(d);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         double best = searchMinimum<D>(iter.position(), d, cell, average);
         double error = std::abs(waveList.kSq()[iter.rank()] - best);
         maxError = std::max(maxError, error/(1.0 + best));
      }
      return maxError;
   }

   /*
   * Largest error of dKSq relative to the central finite difference of
   * kSq, relative to max(1, |kSq|), over all wavevectors and parameters.
   */
   template <int D>
   double finiteDifferenceError(UnitCell<D> const & cell,
                                IntVec<D> const & d)
   {
      Mesh<D> mesh(d);
      WaveList<D,CPT> waveList(false);
      waveList.allocate(mesh, cell);
      waveList.computeMinimumImages();
      waveList.computedKSq();
      int kSize = mesh.size();

      FSArray<double, 6> p = cell.parameters();
      double h = 1.0E-6;
      double maxError = 0.0;
      for (int n = 0; n < cell.nParameter(); ++n) {
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
         for (int i = 0; i < kSize; ++i) {
            double fd = (wp.kSq()[i] - wm.kSq()[i])/(2.0*h);
            double scale = std::max(1.0, std::abs(waveList.kSq()[i]));
            double error = std::abs(waveList.dKSq(n)[i] - fd)/scale;
            maxError = std::max(maxError, error);
         }
      }
      return maxError;
   }

   /*
   * Cells whose metric has symmetry operations that are not signed
   * permutations, on even meshes: kSq is equal for wavevectors related
   * by any symmetry operation of the lattice.
   */
   void testSymmetricCells()
   {
      printMethod(TEST_FUNC);
      double third = 2.0*Constants::Pi/3.0;
      double tolerance = 1.0E-10;

      // Rhombic at 120 and 60 degrees, oblique with a = b at 120
      double pr[2] = {1.7, third};
      UnitCell<2> rhombic;
      makeCell<2>(rhombic, UnitCell<2>::Rhombic, pr, 2);
      TEST_ASSERT(orbitDifference<2>(rhombic, IntVec<2>(12)) < tolerance);
      TEST_ASSERT(orbitDifference<2>(rhombic, IntVec<2>(8)) < tolerance);
      pr[1] = 0.5*third;
      makeCell<2>(rhombic, UnitCell<2>::Rhombic, pr, 2);
      TEST_ASSERT(orbitDifference<2>(rhombic, IntVec<2>(8)) < tolerance);
      double po[3] = {1.7, 1.7, third};
      UnitCell<2> oblique;
      makeCell<2>(oblique, UnitCell<2>::Oblique, po, 3);
      TEST_ASSERT(orbitDifference<2>(oblique, IntVec<2>(10)) < tolerance);

      // Rhombohedral at 60 degrees (fcc) and 109.47 degrees (bcc)
      double ph[2] = {2.0, 0.5*third};
      UnitCell<3> rhombohedral;
      makeCell<3>(rhombohedral, UnitCell<3>::Rhombohedral, ph, 2);
      TEST_ASSERT(orbitDifference<3>(rhombohedral, IntVec<3>(6))
                  < tolerance);
      ph[1] = std::acos(-1.0/3.0);
      makeCell<3>(rhombohedral, UnitCell<3>::Rhombohedral, ph, 2);
      TEST_ASSERT(orbitDifference<3>(rhombohedral, IntVec<3>(6))
                  < tolerance);

      // Monoclinic with a = c and beta = 120 degrees
      double pm[4] = {1.7, 2.1, 1.7, third};
      UnitCell<3> monoclinic;
      makeCell<3>(monoclinic, UnitCell<3>::Monoclinic, pm, 4);
      IntVec<3> dm;
      dm[0] = 6; dm[1] = 4; dm[2] = 6;
      TEST_ASSERT(orbitDifference<3>(monoclinic, dm) < tolerance);

      // Triclinic with a = b, gamma = 120 degrees, c normal to a and b
      double pt[6] = {1.7, 1.7, 2.3, 0.3, 0.0, third};
      UnitCell<3> triclinic;
      makeCell<3>(triclinic, UnitCell<3>::Triclinic, pt, 6);
      IntVec<3> dt;
      dt[0] = 6; dt[1] = 6; dt[2] = 4;
      TEST_ASSERT(orbitDifference<3>(triclinic, dt) < tolerance);
   }

   /*
   * Far from orthogonal geometries (absolute cosines of the angles
   * between reciprocal basis vectors of at least 0.25), kSq is the
   * minimum-image |k|^2; near them (cosines up to 0.1), the average
   * over the signs of the Nyquist indices; dKSq is the derivative of
   * kSq everywhere, also between these ranges.
   */
   void testTransition()
   {
      printMethod(TEST_FUNC);
      double tolerance = 1.0E-10;

      // Minimum image
      double pr[2] = {1.7, 1.92};
      UnitCell<2> rhombic;
      makeCell<2>(rhombic, UnitCell<2>::Rhombic, pr, 2);
      TEST_ASSERT(searchError<2>(rhombic, IntVec<2>(8), false)
                  < tolerance);
      double po[3] = {1.9, 2.3, 2.2};
      UnitCell<2> oblique;
      makeCell<2>(oblique, UnitCell<2>::Oblique, po, 3);
      IntVec<2> d2;
      d2[0] = 8; d2[1] = 6;
      TEST_ASSERT(searchError<2>(oblique, d2, false) < tolerance);
      double pt[6] = {3.0, 4.0, 5.0, 1.8, 0.5, 1.3};
      UnitCell<3> triclinic;
      makeCell<3>(triclinic, UnitCell<3>::Triclinic, pt, 6);
      IntVec<3> d3;
      d3[0] = 6; d3[1] = 8; d3[2] = 10;
      TEST_ASSERT(searchError<3>(triclinic, d3, false) < tolerance);
      double ph[2] = {2.0, 1.05};
      UnitCell<3> rhombohedral;
      makeCell<3>(rhombohedral, UnitCell<3>::Rhombohedral, ph, 2);
      TEST_ASSERT(searchError<3>(rhombohedral, IntVec<3>(6), false)
                  < tolerance);
      double pm[4] = {1.9, 2.1, 2.3, 2.0};
      UnitCell<3> monoclinic;
      makeCell<3>(monoclinic, UnitCell<3>::Monoclinic, pm, 4);
      IntVec<3> dm;
      dm[0] = 8; dm[1] = 6; dm[2] = 10;
      TEST_ASSERT(searchError<3>(monoclinic, dm, false) < tolerance);

      // Rhombic cell at 120 degrees on a 64x8 mesh: the minimum image of
      // (32,0) is (32,-16), shifted along the direction of zero index
      double ph120[2] = {1.0, 2.0*Constants::Pi/3.0};
      UnitCell<2> hexagonal;
      makeCell<2>(hexagonal, UnitCell<2>::Rhombic, ph120, 2);
      IntVec<2> da;
      da[0] = 64; da[1] = 8;
      TEST_ASSERT(searchError<2>(hexagonal, da, false) < tolerance);

      // Average over the signs of the Nyquist indices
      pr[1] = 0.5*Constants::Pi + 0.08;
      makeCell<2>(rhombic, UnitCell<2>::Rhombic, pr, 2);
      TEST_ASSERT(searchError<2>(rhombic, IntVec<2>(8), true)
                  < tolerance);
      pm[3] = 1.65;
      makeCell<3>(monoclinic, UnitCell<3>::Monoclinic, pm, 4);
      TEST_ASSERT(searchError<3>(monoclinic, dm, true) < tolerance);

      // Derivatives at both ends of and within the transition (cells
      // chosen without tied minimum images, where kSq has a kink)
      double c[3] = {0.1, 0.175, 0.25};
      IntVec<2> d10;
      d10[0] = 10; d10[1] = 8;
      for (int k = 0; k < 3; ++k) {
         double pq[3] = {1.7, 1.9, std::acos(-c[k])};
         makeCell<2>(oblique, UnitCell<2>::Oblique, pq, 3);
         TEST_ASSERT(finiteDifferenceError<2>(oblique, d10) < 1.0E-6);
         pm[3] = std::acos(-c[k]);
         makeCell<3>(monoclinic, UnitCell<3>::Monoclinic, pm, 4);
         TEST_ASSERT(finiteDifferenceError<3>(monoclinic, dm) < 1.0E-6);
      }
      double pt2[6] = {1.9, 2.1, 2.3, 0.4, 0.15, 1.62};
      makeCell<3>(triclinic, UnitCell<3>::Triclinic, pt2, 6);
      TEST_ASSERT(finiteDifferenceError<3>(triclinic, d3) < 1.0E-6);
      ph[1] = 1.3;
      makeCell<3>(rhombohedral, UnitCell<3>::Rhombohedral, ph, 2);
      TEST_ASSERT(finiteDifferenceError<3>(rhombohedral, IntVec<3>(8))
                  < 1.0E-6);
   }

   /*
   * Triclinic cell (a, b, c, phi, theta, gamma) = (2, 2.3, 2.7, 0, 0.45,
   * gamma) whose reciprocal basis vector b_1 is orthogonal to b_0 and
   * b_2, while b_0 and b_2 are not (cosine -0.435), 8x8x8 mesh: the 
   * images of wavevectors (4, k, 0) that differ in the sign of the 
   * Nyquist index, and of (4, k, 4) that differ in the signs of both 
   * Nyquist indices, exchange where gamma passes through this value, 
   * and kSq is smooth there: dKSq equals the central difference of kSq.
   */
   void orthogonalPlaneCell(UnitCell<3>& cell)
   {
      double lo = 1.0;
      double hi = 2.2;
      double p[6] = {2.0, 2.3, 2.7, 0.0, 0.45, 0.0};
      for (int it = 0; it < 80; ++it) {
         double mid = 0.5*(lo + hi);
         p[5] = lo;
         makeCell<3>(cell, UnitCell<3>::Triclinic, p, 6);
         bool signLo = metric<3>(cell, 0, 1) > 0.0;
         p[5] = mid;
         makeCell<3>(cell, UnitCell<3>::Triclinic, p, 6);
         bool signMid = metric<3>(cell, 0, 1) > 0.0;
         if (signLo == signMid) {
            lo = mid;
         } else {
            hi = mid;
         }
      }
      p[5] = 0.5*(lo + hi);
      makeCell<3>(cell, UnitCell<3>::Triclinic, p, 6);
   }

   void testOrthogonalPlane()
   {
      printMethod(TEST_FUNC);
      UnitCell<3> cell;
      orthogonalPlaneCell(cell);
      double c02 = metric<3>(cell, 0, 2)
                 / std::sqrt(metric<3>(cell, 0, 0)*metric<3>(cell, 2, 2));
      TEST_ASSERT(std::abs(c02) > 0.4);
      TEST_ASSERT(finiteDifferenceError<3>(cell, IntVec<3>(8)) < 1.0E-6);
   }

};

TEST_BEGIN(CoreNyquistWeightTest)
TEST_ADD(CoreNyquistWeightTest, testSymmetricCells)
TEST_ADD(CoreNyquistWeightTest, testTransition)
TEST_ADD(CoreNyquistWeightTest, testOrthogonalPlane)
TEST_END(CoreNyquistWeightTest)

#endif
