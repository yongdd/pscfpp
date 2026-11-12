#ifndef PSCF_CORE_NYQUIST_GROUP_TEST_H
#define PSCF_CORE_NYQUIST_GROUP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/field/cpp/WaveList.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/crystal/Basis.h>
#include <pscf/crystal/SpaceGroup.h>
#include <pscf/crystal/groupFile.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/containers/FSArray.h>
#include <util/containers/DArray.h>
#include <util/math/Constants.h>

#include <dirent.h>
#include <cmath>
#include <string>
#include <vector>

using namespace Util;
using namespace Pscf;

/*
* Averaging of kSq over the signs of Nyquist indices and space groups.
*
* SpaceGroup::isSignedPermutationGroup identifies groups with 3- or
* 6-fold rotations in a hexagonal setting, under which this average
* is not invariant, and WaveList::allowNyquistAverage(false) makes kSq
* and dKSq the minimum-image values for a variable-angle lattice.
*/
class CoreNyquistGroupTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   /*
   * Names of all standard space groups of dimension D.
   */
   std::vector<std::string> groupNames(int D)
   {
      std::string dir = makeGroupFileName(D, "");
      std::vector<std::string> names;
      DIR* dp = opendir(dir.c_str());
      UTIL_CHECK(dp);
      struct dirent* ep;
      while ((ep = readdir(dp)) != nullptr) {
         std::string name(ep->d_name);
         if (name[0] != '.') names.push_back(name);
      }
      closedir(dp);
      return names;
   }

   /*
   * Is a group in the hexagonal family in a hexagonal setting? (2D:
   * p_3, p_3_m_1, p_3_1_m, p_6, p_6_m_m; 3D: primitive trigonal and
   * hexagonal groups, and rhombohedral groups in hexagonal setting)
   */
   bool isHexagonalSetting(int D, std::string const & name)
   {
      if (D == 2) {
         return (name.compare(0, 3, "p_3") == 0
                 || name.compare(0, 3, "p_6") == 0);
      }
      if (D == 3) {
         if (name.compare(0, 3, "P_3") == 0) return true;
         if (name.compare(0, 4, "P_-3") == 0) return true;
         if (name.compare(0, 3, "P_6") == 0) return true;
         if (name.compare(0, 4, "P_-6") == 0) return true;
         if (name.size() > 2 && name.compare(name.size()-2, 2, ":H") == 0) {
            return true;
         }
      }
      return false;
   }

   template <int D>
   void checkGroups(int nExpected, int nHexExpected)
   {
      std::vector<std::string> names = groupNames(D);
      int nHex = 0;
      for (size_t i = 0; i < names.size(); ++i) {
         SpaceGroup<D> group;
         readGroup(names[i], group);
         bool hex = isHexagonalSetting(D, names[i]);
         if (hex) ++nHex;
         TEST_ASSERT(group.isSignedPermutationGroup() == !hex);
      }
      TEST_ASSERT((int)names.size() == nExpected);
      TEST_ASSERT(nHex == nHexExpected);
   }

   /*
   * All standard groups: only the hexagonal settings contain rotations
   * that are not signed permutations.
   */
   void testSignedPermutationGroups()
   {
      printMethod(TEST_FUNC);
      checkGroups<1>(2, 0);
      checkGroups<2>(17, 5);
      checkGroups<3>(261, 52);
   }

   /*
   * Rhombic cell at 120 degrees on an 8x8 mesh: kSq of the Nyquist wave
   * (4,1) is the average 17 a*^2 by default, that of its p_6_m_m star
   * partner (-1,-3) is 13 a*^2. With allowNyquistAverage(false), kSq
   * and dKSq are the minimum-image values, and the stars of p_6_m_m
   * have equal kSq; allowing it again restores the average.
   */
   void testAllowNyquistAverage()
   {
      printMethod(TEST_FUNC);

      FSArray<double, 6> p;
      p.append(2.0);
      p.append(2.0*Constants::Pi/3.0);
      UnitCell<2> cell;
      cell.set(UnitCell<2>::Rhombic, p);
      IntVec<2> d(8);
      Mesh<2> mesh(d);
      WaveList<2,CPT> waveList(false);
      waveList.allocate(mesh, cell);
      TEST_ASSERT(waveList.hasNyquistAverage());
      waveList.computedKSq();
      DArray<double> averaged;
      averaged.allocate(waveList.kSize());
      for (int i = 0; i < waveList.kSize(); ++i) {
         averaged[i] = waveList.kSq()[i];
      }

      // Wave (4,1) and its star partner (-1,-3) = (7,5)
      IntVec<2> a, b;
      a[0] = 4; a[1] = 1;
      b[0] = 7; b[1] = 5;
      IntVec<2> e0(0);
      e0[0] = 1;
      double aStarSq = cell.ksq(e0);
      double kSqA = waveList.kSq()[mesh.rank(a)];
      double kSqB = waveList.kSq()[mesh.rank(b)];
      TEST_ASSERT(std::abs(kSqA - 17.0*aStarSq) < 1.0E-10*kSqA);
      TEST_ASSERT(std::abs(kSqB - 13.0*aStarSq) < 1.0E-10*kSqB);

      // Forbid averaging: minimum-image values
      waveList.allowNyquistAverage(false);
      TEST_ASSERT(!waveList.hasNyquistAverage());
      TEST_ASSERT(!waveList.hasMinImages());
      TEST_ASSERT(!waveList.hasKSq());
      TEST_ASSERT(!waveList.hasdKSq());
      waveList.computedKSq();
      MeshIterator<2> iter(d);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         int rank = iter.rank();
         IntVec<2> v = shiftToMinimum(iter.position(), d, cell);
         TEST_ASSERT(waveList.minImages()[rank] == v);
         TEST_ASSERT(waveList.kSq()[rank] == cell.ksq(v));
         for (int n = 0; n < cell.nParameter(); ++n) {
            TEST_ASSERT(waveList.dKSq(n)[rank] == cell.dksq(v, n));
         }
      }
      kSqA = waveList.kSq()[mesh.rank(a)];
      kSqB = waveList.kSq()[mesh.rank(b)];
      TEST_ASSERT(std::abs(kSqA - kSqB) < 1.0E-10*kSqA);

      // Stars of p_6_m_m have equal kSq
      SpaceGroup<2> group;
      readGroup("p_6_m_m", group);
      TEST_ASSERT(!group.isSignedPermutationGroup());
      Basis<2> basis;
      basis.makeBasis(mesh, cell, group);
      for (int is = 0; is < basis.nStar(); ++is) {
         Basis<2>::Star const & star = basis.star(is);
         int r0 = mesh.rank(basis.wave(star.beginId).indicesStd);
         for (int iw = star.beginId; iw < star.endId; ++iw) {
            int r = mesh.rank(basis.wave(iw).indicesStd);
            TEST_ASSERT(std::abs(waveList.kSq()[r] - waveList.kSq()[r0])
                        < 1.0E-10*(1.0 + waveList.kSq()[r0]));
         }
      }

      // Allow again: averaged values are restored
      waveList.allowNyquistAverage(true);
      TEST_ASSERT(waveList.hasNyquistAverage());
      waveList.computeKSq();
      for (int i = 0; i < waveList.kSize(); ++i) {
         TEST_ASSERT(waveList.kSq()[i] == averaged[i]);
      }

      // A lattice without a variable angle never averages
      UnitCell<2> hex;
      FSArray<double, 6> ph;
      ph.append(2.0);
      hex.set(UnitCell<2>::Hexagonal, ph);
      WaveList<2,CPT> hexList(true);
      hexList.allocate(mesh, hex);
      TEST_ASSERT(!hexList.hasNyquistAverage());
      hexList.allowNyquistAverage(true);
      TEST_ASSERT(!hexList.hasNyquistAverage());
   }

};

TEST_BEGIN(CoreNyquistGroupTest)
TEST_ADD(CoreNyquistGroupTest, testSignedPermutationGroups)
TEST_ADD(CoreNyquistGroupTest, testAllowNyquistAverage)
TEST_END(CoreNyquistGroupTest)

#endif
