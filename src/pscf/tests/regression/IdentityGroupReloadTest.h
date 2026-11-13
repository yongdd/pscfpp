#ifndef PSCF_IDENTITY_GROUP_RELOAD_TEST_H
#define PSCF_IDENTITY_GROUP_RELOAD_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/crystal/SpaceGroup.h>
#include <pscf/crystal/Basis.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <util/containers/FSArray.h>

using namespace Util;
using namespace Pscf;

/* Reloading the identity must remove every previous symmetry restriction. */
class IdentityGroupReloadTest : public UnitTest
{
public:

   template <int D>
   void checkReload()
   {
      SpaceGroup<D> group, fresh;
      SpaceSymmetry<D> inversion;
      for (int i = 0; i < D; ++i) inversion.R(i,i) = -1;
      group.add(inversion);
      TEST_ASSERT(group.size() == 2);
      TEST_ASSERT(group.isValid());
      readGroup("I", group);
      readGroup("I", fresh);
      TEST_ASSERT(group.size() == 1);
      TEST_ASSERT(group.isValid());
      TEST_ASSERT(group[0] == fresh[0]);
      // Repeating the read is idempotent.
      readGroup("I", group);
      TEST_ASSERT(group.size() == 1);
   }

   void testIdentityReplacement()
   {
      printMethod(TEST_FUNC);
      checkReload<1>();
      checkReload<2>();
      checkReload<3>();
   }

   int basisSize(SpaceGroup<2> const & group)
   {
      UnitCell<2> cell;
      FSArray<double,6> p;
      p.append(2.0);
      cell.set(UnitCell<2>::Square, p);
      IntVec<2> dimensions(8);
      Mesh<2> mesh(dimensions);
      Basis<2> basis;
      basis.makeBasis(mesh, cell, group);
      TEST_ASSERT(basis.isValid());
      return basis.nBasis();
   }

   void testUnrestrictedBasisAfterReload()
   {
      printMethod(TEST_FUNC);
      SpaceGroup<2> reused, fresh, fileControl;
      readGroup("p_2", reused);
      TEST_ASSERT(basisSize(reused) == 34);
      readGroup("I", reused);
      readGroup("I", fresh);
      readGroup("p_2", fileControl);
      readGroup("p_1", fileControl);
      // Identity allows all 64 independent values on this mesh.
      TEST_ASSERT(basisSize(reused) == 64);
      TEST_ASSERT(basisSize(fresh) == 64);
      TEST_ASSERT(basisSize(fileControl) == 64);
   }
};

TEST_BEGIN(IdentityGroupReloadTest)
TEST_ADD(IdentityGroupReloadTest, testIdentityReplacement)
TEST_ADD(IdentityGroupReloadTest, testUnrestrictedBasisAfterReload)
TEST_END(IdentityGroupReloadTest)

#endif
