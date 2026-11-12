#ifndef PSCF_CORE_BASIS_TEST_H
#define PSCF_CORE_BASIS_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/crystal/Basis.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

#include <util/containers/FSArray.h>
#include <util/math/Constants.h>

#include <string>

using namespace Util;
using namespace Pscf;

class CoreBasisTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   /*
   * Groups R_-3_c:H and R_3_c:H (hexagonal setting), in which waves
   * of a star have phase differences from the fractional translations
   * of the c glide. Phase differences must be wrapped by 2 pi: with a
   * wrap by 1.0, R_3_c:H aborted in makeStars and R_-3_c:H gave a
   * wrong number of basis functions.
   *
   * The expected nStar is the number of orbits of the point group on
   * the 12x12x12 mesh, and nBasis is the number of orbits whose waves
   * have no group element that leaves the wave invariant with a phase
   * exp(2 pi i k.t) other than 1 (counted independently of Basis).
   */
   void testMake3DBasis_R_3b_c_H()
   {
      printMethod(TEST_FUNC);

      UnitCell<3> unitCell;
      FSArray<double, 6> parameters;
      parameters.append(1.777035357);
      parameters.append(2.5);
      unitCell.set(UnitCell<3>::Hexagonal, parameters);

      IntVec<3> d;
      d[0] = 12;
      d[1] = 12;
      d[2] = 12;
      Mesh<3> mesh(d);

      Basis<3> basis;
      std::string groupName = "R_-3_c:H";
      basis.makeBasis(mesh, unitCell, groupName);

      TEST_ASSERT(basis.isInitialized());
      TEST_ASSERT(basis.isValid());
      TEST_ASSERT(eq(basis.nWave(), 1728));
      TEST_ASSERT(eq(basis.nStar(), 193));
      TEST_ASSERT(eq(basis.nBasis(), 57));
   }

   void testMake3DBasis_R_3_c_H()
   {
      printMethod(TEST_FUNC);

      UnitCell<3> unitCell;
      FSArray<double, 6> parameters;
      parameters.append(1.777035357);
      parameters.append(2.5);
      unitCell.set(UnitCell<3>::Hexagonal, parameters);

      IntVec<3> d;
      d[0] = 12;
      d[1] = 12;
      d[2] = 12;
      Mesh<3> mesh(d);

      Basis<3> basis;
      std::string groupName = "R_3_c:H";
      basis.makeBasis(mesh, unitCell, groupName);

      TEST_ASSERT(basis.isInitialized());
      TEST_ASSERT(basis.isValid());
      TEST_ASSERT(eq(basis.nWave(), 1728));
      TEST_ASSERT(eq(basis.nStar(), 372));
      TEST_ASSERT(eq(basis.nBasis(), 100));
   }

};

TEST_BEGIN(CoreBasisTest)
TEST_ADD(CoreBasisTest, testMake3DBasis_R_3b_c_H)
TEST_ADD(CoreBasisTest, testMake3DBasis_R_3_c_H)
TEST_END(CoreBasisTest)

#endif
