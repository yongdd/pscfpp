#ifndef PSCF_CORE_TIED_IMAGE_TEST_H
#define PSCF_CORE_TIED_IMAGE_TEST_H

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

using namespace Util;
using namespace Pscf;

/*
* Derivatives dKSq of wavevectors with tied minimum images, CPU WaveList.
*
* In a cell with a variable angle at a symmetric geometry (e.g. a
* rhombic cell at 120 degrees, whose lattice is hexagonal), minimum
* images of some wavevectors lie on the boundary of the Voronoi cell of
* the lattice of image shifts, on meshes of any size, and are tied with
* other images whose |k|^2 has a different derivative. dKSq was taken
* from one image (the largest in the order of IntVec), so that the stress
* of a state with the symmetry of the cell was not symmetric: a flexible
* cell of such a state moved away from the symmetric geometry. dKSq is
* now averaged over the tied images.
*
* Check: for weights w(kSq) (which have the symmetry of the cell), the
* tensor sum of w k k over all wavevectors is isotropic, so that the sum
* S_n of w dKSq_n is proportional to the derivative of the logarithm of
* the cell volume with respect to parameter n.
*/
class CoreTiedImageTest : public UnitTest
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
   * Sums S[n] of w(kSq) dKSq_n over all wavevectors of a real-field
   * WaveList (dKSq contains the factor 2 of implicit wavevectors),
   * with w = exp(-kSq/<kSq>).
   */
   template <int D>
   void weightedSums(UnitCell<D> const & cell, IntVec<D> const & d,
                     double* S)
   {
      Mesh<D> mesh(d);
      WaveList<D,CPT> waveList(true);
      waveList.allocate(mesh, cell);
      waveList.computeMinimumImages();
      waveList.computedKSq();
      int kSize = waveList.kSize();
      double mean = 0.0;
      for (int i = 0; i < kSize; ++i) {
         mean += waveList.kSq()[i];
      }
      mean /= double(kSize);
      for (int n = 0; n < cell.nParameter(); ++n) {
         S[n] = 0.0;
         for (int i = 0; i < kSize; ++i) {
            double w = std::exp(-waveList.kSq()[i]/mean);
            S[n] += w*waveList.dKSq(n)[i];
         }
      }
   }

   /*
   * Deviatoric residual of S relative to max |S[n]|, for a cell whose
   * metric is invariant under the symmetry operations of the cell.
   * dLogV[n] is the derivative of the log of the volume with respect
   * to parameter n, for the parameters in the list (the residual of
   * parameter n is S[n] - S[n0] dLogV[n]/dLogV[n0], with n0 the first).
   */
   template <int D>
   double residual(UnitCell<D> const & cell, IntVec<D> const & d,
                   double const * dLogV, int const * params, int nCheck)
   {
      double S[6];
      weightedSums<D>(cell, d, S);
      double scale = 0.0;
      for (int n = 0; n < cell.nParameter(); ++n) {
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

   /*
   * Generic cell (no tied images with different derivatives, and no
   * Nyquist average, whose weights vanish far from right angles): dKSq
   * is the derivative of |k|^2 of the stored image, bit for bit.
   */
   template <int D>
   void checkUnchanged(UnitCell<D> const & cell, IntVec<D> const & d)
   {
      Mesh<D> mesh(d);
      WaveList<D,CPT> waveList(true);
      waveList.allocate(mesh, cell);
      waveList.computeMinimumImages();
      waveList.computedKSq();
      DArray<bool> const & implicit = waveList.implicitInverse();
      for (int i = 0; i < waveList.kSize(); ++i) {
         IntVec<D> const & v = waveList.minImages()[i];
         for (int n = 0; n < cell.nParameter(); ++n) {
            double dksq = cell.dksq(v, n);
            if (implicit[i]) dksq *= 2.0;
            TEST_ASSERT(waveList.dKSq(n)[i] == dksq);
         }
      }
   }

   /*
   * Rhombic cells at 120 and 60 degrees and an oblique cell with a = b
   * at 120 degrees (hexagonal lattices), on even and odd meshes; a
   * generic rhombic cell is unchanged.
   */
   void testHexagonal()
   {
      printMethod(TEST_FUNC);
      double third = 2.0*Constants::Pi/3.0;
      double tolerance = 1.0E-12;
      int params2[2] = {0, 1};

      // Rhombic (a, gamma): area a^2 sin(gamma)
      double a = 1.7;
      double p[2] = {a, third};
      UnitCell<2> rhombic;
      double dLogV[3];
      int meshes[3] = {12, 13, 14};
      for (int m = 0; m < 3; ++m) {
         makeCell<2>(rhombic, UnitCell<2>::Rhombic, p, 2);
         dLogV[0] = 2.0/a;
         dLogV[1] = std::cos(third)/std::sin(third);
         TEST_ASSERT(residual<2>(rhombic, IntVec<2>(meshes[m]), dLogV,
                                 params2, 2) < tolerance);
      }
      p[1] = 0.5*third;
      makeCell<2>(rhombic, UnitCell<2>::Rhombic, p, 2);
      dLogV[1] = std::cos(0.5*third)/std::sin(0.5*third);
      TEST_ASSERT(residual<2>(rhombic, IntVec<2>(15), dLogV, params2, 2)
                  < tolerance);

      // Oblique (a, b, gamma): area a b sin(gamma)
      double po[3] = {a, a, third};
      UnitCell<2> oblique;
      makeCell<2>(oblique, UnitCell<2>::Oblique, po, 3);
      dLogV[0] = 1.0/a;
      dLogV[1] = 1.0/a;
      dLogV[2] = std::cos(third)/std::sin(third);
      int params3[3] = {0, 1, 2};
      TEST_ASSERT(residual<2>(oblique, IntVec<2>(12), dLogV, params3, 3)
                  < tolerance);

      // Generic rhombic cell
      p[1] = 1.9;
      makeCell<2>(rhombic, UnitCell<2>::Rhombic, p, 2);
      checkUnchanged<2>(rhombic, IntVec<2>(12));
      checkUnchanged<2>(rhombic, IntVec<2>(13));
   }

   /*
   * Rhombohedral cells at 60 degrees (fcc) and 109.47 degrees (bcc),
   * and a monoclinic cell with a = c and beta = 120 degrees; a generic
   * triclinic cell is unchanged.
   */
   void testThreeDimensional()
   {
      printMethod(TEST_FUNC);
      double tolerance = 1.0E-12;

      // Rhombohedral (a, alpha): V = a^3 (1 - 3c^2 + 2c^3)^(1/2)
      double a = 2.0;
      double alpha[2] = {Constants::Pi/3.0, std::acos(-1.0/3.0)};
      int params2[2] = {0, 1};
      for (int k = 0; k < 2; ++k) {
         double p[2] = {a, alpha[k]};
         UnitCell<3> cell;
         makeCell<3>(cell, UnitCell<3>::Rhombohedral, p, 2);
         double c = std::cos(alpha[k]);
         double dLogV[2];
         dLogV[0] = 3.0/a;
         dLogV[1] = 3.0*std::sin(alpha[k])*c*(1.0 - c)
                    /(1.0 - 3.0*c*c + 2.0*c*c*c);
         TEST_ASSERT(residual<3>(cell, IntVec<3>(8), dLogV, params2, 2)
                     < tolerance);
      }

      // Monoclinic (a, b, c, beta) with a = c, beta = 120 degrees: the
      // (a, c) plane is hexagonal, b is normal to it; the residuals of
      // c and beta are relative to a (area of the (a, c) plane)
      double third = 2.0*Constants::Pi/3.0;
      double pm[4] = {1.7, 2.1, 1.7, third};
      UnitCell<3> monoclinic;
      makeCell<3>(monoclinic, UnitCell<3>::Monoclinic, pm, 4);
      double dLogA[4];
      dLogA[0] = 1.0/pm[0];
      dLogA[1] = 0.0;
      dLogA[2] = 1.0/pm[2];
      dLogA[3] = std::cos(third)/std::sin(third);
      int params[3] = {0, 2, 3};
      IntVec<3> d;
      d[0] = 12; d[1] = 4; d[2] = 12;
      TEST_ASSERT(residual<3>(monoclinic, d, dLogA, params, 3)
                  < tolerance);
   }

};

TEST_BEGIN(CoreTiedImageTest)
TEST_ADD(CoreTiedImageTest, testHexagonal)
TEST_ADD(CoreTiedImageTest, testThreeDimensional)
TEST_END(CoreTiedImageTest)

#endif
