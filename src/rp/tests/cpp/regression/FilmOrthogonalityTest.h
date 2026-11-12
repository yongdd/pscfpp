#ifndef RPC_FILM_ORTHOGONALITY_TEST_H
#define RPC_FILM_ORTHOGONALITY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/environment/Environment.h>

#include <util/containers/FSArray.h>
#include <util/math/Constants.h>
#include <util/misc/Exception.h>
#include <util/tests/LogFileUnitTest.h>

#include <sstream>
#include <string>

using namespace Util;
using namespace Pscf;

/*
* Regression tests for the orthogonality check of thin-film unit cells.
*
* The lattice vectors parallel to the walls must be orthogonal to the
* wall normal. The check must not depend on the unit of length: a right
* angle written with the 11 significant digits of a field file header
* (1.5707963268) is accepted at any length scale, and a cell that is
* not orthogonal is rejected at any length scale.
*/
class FilmOrthogonalityTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Generate the fields of a film with a given unit cell.
   *
   * The film (group p_1 or P_1, no Iterator) has walls normal to lattice
   * vector normalVecId. Returns true iff an Exception is thrown.
   */
   template <int D>
   bool generateThrows(std::string const & lattice, int normalVecId,
                       FSArray<double, 6> const & parameters)
   {
      std::string mesh = (D == 2) ? "8 8" : "4 4 4";
      std::string group = (D == 2) ? "p_1" : "P_1";

      std::ostringstream param;
      param << "System{\n"
            << "  Mixture{\n"
            << "    nMonomer  2\n"
            << "    monomers[  1.0  1.0  ]\n"
            << "    nPolymer  1\n"
            << "    Polymer{\n"
            << "      type    linear\n"
            << "      nBlock  2\n"
            << "      blocks[\n"
            << "              0  0.5\n"
            << "              1  0.5\n"
            << "      ]\n"
            << "      phi     1.0\n"
            << "    }\n"
            << "    ds   0.1\n"
            << "  }\n"
            << "  Interaction{\n"
            << "    chi(\n"
            << "         1   0   10.0\n"
            << "    )\n"
            << "  }\n"
            << "  Domain{\n"
            << "    mesh           " << mesh << "\n"
            << "    lattice        " << lattice << "\n"
            << "    groupName      " << group << "\n"
            << "  }\n"
            << "  FilmEnvironment{\n"
            << "    normalVecId        " << normalVecId << "\n"
            << "    interfaceThickness 0.2\n"
            << "    excludedThickness  0.4\n"
            << "    chiBottom[  0.0  0.0  ]\n"
            << "    chiTop[  0.0  0.0  ]\n"
            << "  }\n"
            << "}\n";

      Rp::System<D,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream in(param.str());
      system.readParam(in);
      system.setUnitCell(parameters);
      TEST_ASSERT(system.hasEnvironment());

      bool thrown = false;
      try {
         system.environment().generate();
      } catch (Exception& e) {
         thrown = true;
      }
      return thrown;
   }

   // Cell parameters
   FSArray<double, 6> cell(double p0, double p1, double p2,
                           double p3 = -1.0, double p4 = -1.0,
                           double p5 = -1.0)
   {
      FSArray<double, 6> parameters;
      parameters.append(p0);
      parameters.append(p1);
      parameters.append(p2);
      if (p3 > -1.0) parameters.append(p3);
      if (p4 > -1.0) parameters.append(p4);
      if (p5 > -1.0) parameters.append(p5);
      return parameters;
   }

   /*
   * Oblique 2D films with a right angle rounded to 11 significant digits
   * (cos = -5.1e-12) or given exactly, with lengths of order 1 and 100.
   */
   void testRoundedRightAngle2D()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmOrthogonalityRounded2D.log");
      double pi = Constants::Pi;
      double right = 1.5707963268;

      for (int nv = 0; nv < 2; ++nv) {
         TEST_ASSERT(!generateThrows<2>("oblique", nv, cell(2.0, 2.5, right)));
         TEST_ASSERT(!generateThrows<2>("oblique", nv, cell(100.0, 125.0, right)));
         TEST_ASSERT(!generateThrows<2>("oblique", nv,
                                        cell(100.0, 125.0, 0.5*pi)));
      }
   }

   /*
   * Triclinic and monoclinic 3D films with right angles rounded to 11
   * significant digits (up, cos = -5.1e-12, or down, cos = 9.5e-11),
   * with lengths of order 1 and 50.
   */
   void testRoundedRightAngle3D()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmOrthogonalityRounded3D.log");
      double right = 1.5707963268;
      double down = 1.5707963267;

      for (int k = 0; k < 2; ++k) {
         double s = (k == 0) ? 1.0 : 25.0;

         // Walls normal to a: phi = gamma = pi/2
         TEST_ASSERT(!generateThrows<3>("triclinic", 0,
                     cell(2.0*s, 2.5*s, 2.0*s, right, 0.25, right)));

         // Walls normal to b: phi = 0, gamma = pi/2
         TEST_ASSERT(!generateThrows<3>("triclinic", 1,
                     cell(2.0*s, 2.5*s, 2.0*s, 0.0, 0.25, right)));
         TEST_ASSERT(!generateThrows<3>("triclinic", 1,
                     cell(2.0*s, 2.5*s, 2.0*s, 0.0, 0.25, down)));

         // Walls normal to c: beta = pi/2
         TEST_ASSERT(!generateThrows<3>("monoclinic", 2,
                     cell(2.0*s, 2.0*s, 2.5*s, right)));
      }
   }

   /*
   * Cells whose angles differ from pi/2 by 1e-6 are rejected also when
   * the lattice vectors are short (dot products below 1e-8).
   */
   void testSmallNonOrthogonalCells()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmOrthogonalitySmall.log");
      double pi = Constants::Pi;
      double angle = 0.5*pi + 1.0E-6;
      double s = 0.025;

      TEST_ASSERT(generateThrows<2>("oblique", 0,
                  cell(2.0*s, 2.5*s, angle)));
      TEST_ASSERT(generateThrows<3>("triclinic", 1,
                  cell(2.0*s, 2.5*s, 2.0*s, 0.0, 0.25, angle)));
      TEST_ASSERT(generateThrows<3>("monoclinic", 2,
                  cell(2.0*s, 2.0*s, 2.5*s, angle)));

      // The same cells with lengths of order 1 are rejected
      TEST_ASSERT(generateThrows<2>("oblique", 0, cell(2.0, 2.5, angle)));
      TEST_ASSERT(generateThrows<3>("monoclinic", 2,
                  cell(2.0, 2.0, 2.5, angle)));
   }

};

TEST_BEGIN(FilmOrthogonalityTest)
TEST_ADD(FilmOrthogonalityTest, testRoundedRightAngle2D)
TEST_ADD(FilmOrthogonalityTest, testRoundedRightAngle3D)
TEST_ADD(FilmOrthogonalityTest, testSmallNonOrthogonalCells)
TEST_END(FilmOrthogonalityTest)

#endif
