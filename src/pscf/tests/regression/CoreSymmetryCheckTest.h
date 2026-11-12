#ifndef PSCF_CORE_SYMMETRY_CHECK_TEST_H
#define PSCF_CORE_SYMMETRY_CHECK_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/fieldIo/rFieldIo.h>
#include <pscf/field/cpp/FFT.h>
#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldDft.h>
#include <pscf/crystal/Basis.h>
#include <pscf/crystal/SpaceGroup.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/DArray.h>
#include <util/containers/FSArray.h>
#include <util/math/Constants.h>
#include <util/random/Random.h>

#include <complex>
#include <string>

using namespace Util;
using namespace Pscf;

/*
* Pscf::hasSymmetry (CHECK_RGRID_SYMMETRY, the warning of KGRID_TO_BASIS
* and RGRID_TO_BASIS) must detect fields that break a space group
* without inversion, whose stars come in open pairs (invertFlag +1/-1).
*/
class CoreSymmetryCheckTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   /*
   * Is a wave stored explicitly in the DFT of a real field?
   */
   template <int D>
   bool isExplicit(Basis<D> const & basis, int iw)
   {  return !basis.wave(iw).implicit; }

   /*
   * Make the DFT of a field with the symmetry of the basis, with
   * random coefficients.
   */
   template <int D>
   void makeSymmetric(Basis<D> const & basis,
                      IntVec<D> const & meshDims,
                      RFieldDft<D,CPT>& kField)
   {
      Random random;
      random.setSeed(73);
      DArray<double> components;
      components.allocate(basis.nBasis());
      for (int i = 0; i < basis.nBasis(); ++i) {
         components[i] = random.uniform(-1.0, 1.0);
      }
      kField.allocate(meshDims);
      convertBasisToKGrid(components, kField, basis,
                          kField.dftDimensions());
   }

   /*
   * Multiply the explicit waves of the first star of each open pair by
   * exp(i*phase), for every pair in which no wave of either star has
   * an explicit inverse (so that the result is still the DFT of a real
   * field). Returns the number of pairs changed.
   */
   template <int D>
   int twistOpenPairs(Basis<D> const & basis,
                      RFieldDft<D,CPT>& kField, double phase)
   {
      Mesh<D> dftMesh(kField.dftDimensions());
      std::complex<double> factor = std::polar(1.0, phase);
      int nTwist = 0;
      for (int is = 0; is < basis.nStar() - 1; ++is) {
         typename Basis<D>::Star const & s0 = basis.star(is);
         typename Basis<D>::Star const & s1 = basis.star(is+1);
         if (s0.cancel || s0.invertFlag != 1) continue;
         TEST_ASSERT(s1.invertFlag == -1);
         bool explicit0 = false;
         bool explicit1 = false;
         bool linked = false;
         for (int iw = s0.beginId; iw < s0.endId; ++iw) {
            if (isExplicit(basis, iw)) explicit0 = true;
            int jw = basis.wave(iw).inverseId;
            if (isExplicit(basis, jw)) explicit1 = true;
            if (isExplicit(basis, iw) && isExplicit(basis, jw)) {
               linked = true;
            }
         }
         if (!explicit0 || !explicit1 || linked) continue;
         for (int iw = s0.beginId; iw < s0.endId; ++iw) {
            if (!isExplicit(basis, iw)) continue;
            int rank = dftMesh.rank(basis.wave(iw).indicesStd);
            std::complex<double> c(kField[rank][0], kField[rank][1]);
            c *= factor;
            kField[rank][0] = c.real();
            kField[rank][1] = c.imag();
         }
         ++nTwist;
      }
      return nTwist;
   }

   /*
   * For a group without inversion: the symmetric field passes, the
   * field with twisted open pairs fails.
   */
   template <int D>
   void checkGroup(std::string groupName,
                   typename UnitCell<D>::LatticeSystem lattice,
                   FSArray<double, 6> const & parameters,
                   IntVec<D> const & meshDims)
   {
      UnitCell<D> cell;
      cell.set(lattice, parameters);
      Mesh<D> mesh(meshDims);
      SpaceGroup<D> group;
      readGroup(groupName, group);
      Basis<D> basis;
      basis.makeBasis(mesh, cell, group);
      TEST_ASSERT(basis.isValid());

      RFieldDft<D,CPT> kField;
      makeSymmetric<D>(basis, meshDims, kField);
      TEST_ASSERT(hasSymmetry(kField, basis, kField.dftDimensions(),
                              1.0E-8, false));

      int nTwist = twistOpenPairs<D>(basis, kField, 0.9);
      TEST_ASSERT(nTwist > 0);
      TEST_ASSERT(!hasSymmetry(kField, basis, kField.dftDimensions(),
                               1.0E-8, false));
   }

   /*
   * Hexagonal p_3 on a 12x12 mesh, with the field
   * w(i,j) = 1 + 0.3 cos(2 pi (i - 2j)/12 + 0.4), which the 3-fold
   * rotation changes by up to 0.6. A p_3 symmetric field passes.
   */
   void testP3Hexagonal()
   {
      printMethod(TEST_FUNC);

      IntVec<2> d;
      d[0] = 12;
      d[1] = 12;
      FSArray<double, 6> p;
      p.append(1.7);
      UnitCell<2> cell;
      cell.set(UnitCell<2>::Hexagonal, p);
      Mesh<2> mesh(d);
      SpaceGroup<2> group;
      readGroup("p_3", group);
      Basis<2> basis;
      basis.makeBasis(mesh, cell, group);
      TEST_ASSERT(basis.isValid());

      // Non-symmetric field
      RField<2,CPT> rField;
      rField.allocate(d);
      double twoPi = 2.0*Constants::Pi;
      MeshIterator<2> iter(d);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         int i = iter.position(0);
         int j = iter.position(1);
         rField[iter.rank()] = 1.0
                             + 0.3*cos(twoPi*(i - 2*j)/12.0 + 0.4);
      }
      FFT<2,CPT> fft;
      fft.setup(d);
      RFieldDft<2,CPT> kField;
      kField.allocate(d);
      fft.forwardTransform(rField, kField);
      TEST_ASSERT(!hasSymmetry(kField, basis, kField.dftDimensions(),
                               1.0E-8, false));

      // Symmetric field
      RFieldDft<2,CPT> kSym;
      makeSymmetric<2>(basis, d, kSym);
      TEST_ASSERT(hasSymmetry(kSym, basis, kSym.dftDimensions(),
                              1.0E-8, false));
   }

   /*
   * Twisted open pairs in 2D and 3D groups without inversion.
   */
   void testTwistedOpenPairs()
   {
      printMethod(TEST_FUNC);

      // 2D: p_3_m_1 hexagonal
      {
         FSArray<double, 6> p;
         p.append(1.7);
         IntVec<2> d(12);
         checkGroup<2>("p_3_m_1", UnitCell<2>::Hexagonal, p, d);
      }

      // 3D cubic: I_41_3_2 and P_21_3
      {
         FSArray<double, 6> p;
         p.append(3.5);
         IntVec<3> d(12);
         checkGroup<3>("I_41_3_2", UnitCell<3>::Cubic, p, d);
         checkGroup<3>("P_21_3", UnitCell<3>::Cubic, p, d);
      }

      // 3D orthorhombic: P_2_2_2
      {
         FSArray<double, 6> p;
         p.append(2.0);
         p.append(2.4);
         p.append(2.8);
         IntVec<3> d;
         d[0] = 10;
         d[1] = 12;
         d[2] = 14;
         checkGroup<3>("P_2_2_2", UnitCell<3>::Orthorhombic, p, d);
      }

      // 3D monoclinic: P_1_2_1
      {
         FSArray<double, 6> p;
         p.append(2.0);
         p.append(2.4);
         p.append(2.8);
         p.append(1.9);
         IntVec<3> d(12);
         checkGroup<3>("P_1_2_1", UnitCell<3>::Monoclinic, p, d);
      }
   }

   /*
   * Centrosymmetric groups (closed stars only) and P_1 (where no
   * pair is unlinked) are unaffected: symmetric fields pass.
   */
   void testSymmetricFieldsPass()
   {
      printMethod(TEST_FUNC);

      FSArray<double, 6> p;
      p.append(3.5);
      std::string names[3] = {"I_a_-3_d", "P_21_3", "I_41_3_2"};
      for (int i = 0; i < 3; ++i) {
         IntVec<3> d(12);
         UnitCell<3> cell;
         cell.set(UnitCell<3>::Cubic, p);
         Mesh<3> mesh(d);
         SpaceGroup<3> group;
         readGroup(names[i], group);
         Basis<3> basis;
         basis.makeBasis(mesh, cell, group);
         RFieldDft<3,CPT> kField;
         makeSymmetric<3>(basis, d, kField);
         TEST_ASSERT(hasSymmetry(kField, basis, kField.dftDimensions(),
                                 1.0E-8, false));
      }
   }

};

TEST_BEGIN(CoreSymmetryCheckTest)
TEST_ADD(CoreSymmetryCheckTest, testP3Hexagonal)
TEST_ADD(CoreSymmetryCheckTest, testTwistedOpenPairs)
TEST_ADD(CoreSymmetryCheckTest, testSymmetricFieldsPass)
TEST_END(CoreSymmetryCheckTest)

#endif
