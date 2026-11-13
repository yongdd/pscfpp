#ifndef PSCF_BASIS_SCALE_TEST_H
#define PSCF_BASIS_SCALE_TEST_H

#include <test/UnitTestRunner.h>
#include <pscf/crystal/Basis.h>
#include <pscf/crystal/SpaceGroup.h>
#include <pscf/crystal/unitCellHeader.h>
#include <util/math/Constants.h>
#include <util/misc/Exception.h>
#include <set>
#include <sstream>

using namespace Util;
using namespace Pscf;

class BasisScaleTest : public UnitTest
{
public:

   void testHexagonalBasisUnits()
   {
      printMethod(TEST_FUNC);
      IntVec<2> dims(24);
      Mesh<2> mesh(dims);
      SpaceGroup<2> group;
      readGroup("p_6_m_m", group);
      // All translations are zero: count independent plane-wave orbits
      // directly with integer group actions, without squared norms.
      std::set<int> roots;
      MeshIterator<2> iter(dims);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         int root = mesh.size();
         for (int j = 0; j < group.size(); ++j) {
            IntVec<2> v = iter.position()*group[j];
            mesh.shift(v);
            int rank = v[0]*dims[1]+v[1];
            root = std::min(root, rank);
         }
         roots.insert(root);
      }
      TEST_ASSERT(roots.size() == 61);
      double scales[] = {1e-6, 1e-4, .01, 1.0, 1e6};
      for (double scale : scales) {
         FSArray<double,6> p;
         p.append(scale);
         UnitCell<2> cell;
         cell.set(UnitCell<2>::Hexagonal, p);
         Basis<2> basis;
         basis.makeBasis(mesh, cell, group);
         TEST_ASSERT(basis.isValid());
         TEST_ASSERT(basis.nBasis() == int(roots.size()));
         TEST_ASSERT(basis.wave(0).indicesStd == IntVec<2>(0));
      }
   }

   void testRoundedAngleBasisUnits()
   {
      printMethod(TEST_FUNC);
      double right = 1.5707963268;
      double down = 1.5707963267;
      double scales[] = {1e-6, 1.0, 25.0, 1e6};
      IntVec<3> dims(4);
      Mesh<3> mesh(dims);
      for (double scale : scales) {
         for (int kind = 0; kind < 4; ++kind) {
            FSArray<double,6> p;
            p.append(2.0*scale);
            p.append((kind == 3 ? 2.0 : 2.5)*scale);
            p.append((kind == 3 ? 2.5 : 2.0)*scale);
            UnitCell<3> cell;
            if (kind == 3) {
               p.append(right);
               cell.set(UnitCell<3>::Monoclinic, p);
            } else {
               p.append(kind == 0 ? right : 0.0);
               p.append(0.25);
               p.append(kind == 2 ? down : right);
               cell.set(UnitCell<3>::Triclinic, p);
            }
            Basis<3> basis;
            basis.makeBasis(mesh, cell, "I");
            TEST_ASSERT(basis.isValid());
            TEST_ASSERT(basis.nBasis() == mesh.size());
            MeshIterator<3> iter(dims);
            for (iter.begin(); !iter.atEnd(); ++iter) {
               IntVec<3> inv;
               inv.negate(iter.position());
               mesh.shift(inv);
               int id = basis.waveId(iter.position());
               TEST_ASSERT(basis.wave(id).inverseId == basis.waveId(inv));
            }
         }
      }
   }

   void testSerializedHexagonalAngle()
   {
      printMethod(TEST_FUNC);
      Mesh<2> mesh(IntVec<2>(24));
      double scales[] = {1e-6, .01, 1.0, 1e6};
      for (double scale : scales) {
         FSArray<double,6> p;
         p.append(scale);
         p.append(2.0*Constants::Pi/3.0);
         UnitCell<2> original, restored;
         original.set(UnitCell<2>::Rhombic, p);
         std::stringstream header;
         writeUnitCellHeader(header, original);
         readUnitCellHeader(header, restored);
         Basis<2> basis;
         basis.makeBasis(mesh, restored, "p_6_m_m");
         TEST_ASSERT(basis.isValid());
         TEST_ASSERT(basis.nBasis() == 61);

         // A resolvably incompatible angle must still be rejected.
         p[1] += 1e-5;
         UnitCell<2> incompatible;
         incompatible.set(UnitCell<2>::Rhombic, p);
         bool rejected = false;
         try {
            Basis<2> invalid;
            invalid.makeBasis(mesh, incompatible, "p_6_m_m");
         } catch (Exception&) { rejected = true; }
         TEST_ASSERT(rejected);
      }
   }
};

TEST_BEGIN(BasisScaleTest)
TEST_ADD(BasisScaleTest, testHexagonalBasisUnits)
TEST_ADD(BasisScaleTest, testRoundedAngleBasisUnits)
TEST_ADD(BasisScaleTest, testSerializedHexagonalAngle)
TEST_END(BasisScaleTest)

#endif
