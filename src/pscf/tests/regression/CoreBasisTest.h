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

   /*
   * Groups without inversion on skewed cells with anisotropic meshes.
   * Basis construction aborted, first because the minimum image search
   * window was too small, then because the rule that chooses among
   * tied minimum images depended on the representative passed to
   * shiftToMinimum (Basis compares the images of k and -k). Every
   * star of P_1 has one wave; the stars of P_1_2_1 on a 4x4x32 mesh
   * are the (512 + 16)/2 orbits of the 2-fold axis along b.
   */
   void testMake2DBasis_p_1_oblique()
   {
      printMethod(TEST_FUNC);

      UnitCell<2> unitCell;
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      parameters.append(2.0);
      parameters.append(2.8);
      unitCell.set(UnitCell<2>::Oblique, parameters);

      IntVec<2> d;
      d[0] = 8;
      d[1] = 64;
      Mesh<2> mesh(d);

      Basis<2> basis;
      std::string groupName = "p_1";
      basis.makeBasis(mesh, unitCell, groupName);

      TEST_ASSERT(basis.isInitialized());
      TEST_ASSERT(basis.isValid());
      TEST_ASSERT(eq(basis.nWave(), 512));
      TEST_ASSERT(eq(basis.nStar(), 512));
      TEST_ASSERT(eq(basis.nBasis(), 512));
   }

   void testMake3DBasis_P_1_2_1_monoclinic()
   {
      printMethod(TEST_FUNC);

      UnitCell<3> unitCell;
      FSArray<double, 6> parameters;
      parameters.append(1.0);
      parameters.append(1.0);
      parameters.append(1.0);
      parameters.append(2.8);
      unitCell.set(UnitCell<3>::Monoclinic, parameters);

      IntVec<3> d;
      d[0] = 4;
      d[1] = 4;
      d[2] = 32;
      Mesh<3> mesh(d);

      Basis<3> basis;
      std::string groupName = "P_1_2_1";
      basis.makeBasis(mesh, unitCell, groupName);

      TEST_ASSERT(basis.isInitialized());
      TEST_ASSERT(basis.isValid());
      TEST_ASSERT(eq(basis.nWave(), 512));
      TEST_ASSERT(eq(basis.nStar(), 264));
      TEST_ASSERT(eq(basis.nBasis(), 264));
   }

   void testMake3DBasis_P_1_triclinic()
   {
      printMethod(TEST_FUNC);

      UnitCell<3> unitCell;
      FSArray<double, 6> parameters;
      parameters.append(1.0);
      parameters.append(1.0);
      parameters.append(1.0);
      parameters.append(0.0);
      parameters.append(1.2);
      parameters.append(0.5*Constants::Pi);
      unitCell.set(UnitCell<3>::Triclinic, parameters);

      IntVec<3> d;
      d[0] = 4;
      d[1] = 4;
      d[2] = 32;
      Mesh<3> mesh(d);

      Basis<3> basis;
      std::string groupName = "P_1";
      basis.makeBasis(mesh, unitCell, groupName);

      TEST_ASSERT(basis.isInitialized());
      TEST_ASSERT(basis.isValid());
      TEST_ASSERT(eq(basis.nWave(), 512));
      TEST_ASSERT(eq(basis.nStar(), 512));
      TEST_ASSERT(eq(basis.nBasis(), 512));
   }

};

TEST_BEGIN(CoreBasisTest)
TEST_ADD(CoreBasisTest, testMake3DBasis_R_3b_c_H)
TEST_ADD(CoreBasisTest, testMake3DBasis_R_3_c_H)
TEST_ADD(CoreBasisTest, testMake2DBasis_p_1_oblique)
TEST_ADD(CoreBasisTest, testMake3DBasis_P_1_2_1_monoclinic)
TEST_ADD(CoreBasisTest, testMake3DBasis_P_1_triclinic)
TEST_END(CoreBasisTest)

#endif
