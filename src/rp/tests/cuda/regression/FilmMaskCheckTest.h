#ifndef RPG_FILM_MASK_CHECK_TEST_H
#define RPG_FILM_MASK_CHECK_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/environment/FilmFieldGenMask.h>
#include <rp/scft/iterator/Iterator.h>
#include <rp/system/System.h>
#include <rp/field/Domain.h>

#include <pscf/crystal/UnitCell.h>

#include <util/containers/FSArray.h>
#include <util/math/Constants.h>
#include <util/misc/Exception.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;

/*
* Regression tests for the compatibility checks of thin-film masks.
*/
class FilmMaskCheckTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   // Read parameter file to create a System object
   template <int D>
   void createSystem(Rp::System<D,CUT>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

   // Read parameter file section to create a FilmFieldGenMask object
   template <int D>
   void createFilmFieldGenMask(Rp::FilmFieldGenMask<D,CUT>& mask,
                               std::string fname)
   {
      std::ifstream in;
      openInputFile(fname, in);
      mask.readParameters(in);
      in.close();
   }

   /*
   * Check the compatibility of a mask with a triclinic system.
   *
   * The system (group P_1) has cell parameters (a, b, c, phi, theta,
   * gamma) = (2.0, 2.2, 2.4, phi, theta, gamma), in which b lies in the
   * x-y plane at angle gamma from a, and c has polar angle theta and
   * azimuthal angle phi. All parameters are flexible, unless rigid is
   * true. The walls are normal to lattice vector normalVecId.
   *
   * Returns true iff checkCompatibility throws an Exception. On return,
   * flexible contains the flexible parameters of the iterator.
   */
   bool checkTriclinicFilm(double phi, double theta, double gamma,
                           int normalVecId, bool rigid,
                           FSArray<bool,6>& flexible)
   {
      Rp::System<3,CUT> system;
      createSystem(system, "in/Film/tric3D_param");
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      parameters.append(2.2);
      parameters.append(2.4);
      parameters.append(phi);
      parameters.append(theta);
      parameters.append(gamma);
      system.setUnitCell(parameters);
      if (rigid) {
         FSArray<bool,6> none;
         for (int i = 0; i < 6; ++i) {
            none.append(false);
         }
         system.iterator().setFlexibleParams(none);
      }

      Rp::FilmFieldGenMask<3,CUT> mask(system);
      if (normalVecId == 0) {
         createFilmFieldGenMask(mask, "in/Film/mask_nv0");
      } else if (normalVecId == 1) {
         createFilmFieldGenMask(mask, "in/Film/mask_nv1");
      } else {
         createFilmFieldGenMask(mask, "in/Film/mask_nv2");
      }
      TEST_ASSERT(mask.normalVecId() == normalVecId);

      bool thrown = false;
      try {
         mask.checkCompatibility();
      } catch (Exception& e) {
         thrown = true;
      }
      flexible = system.iterator().flexibleParams();
      return thrown;
   }

   /*
   * The lattice vectors parallel to the walls must be orthogonal to the
   * wall normal: dot products of either sign must be rejected (3D).
   */
   void testCheckLatticeVectorsSign()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmCheckLatticeVectorsSign.log");
      double pi = Constants::Pi;
      FSArray<bool,6> flex;

      // Walls normal to c, with a.c < 0
      TEST_ASSERT(checkTriclinicFilm(pi, 0.3, 0.5*pi, 2, true, flex));

      // Walls normal to b, with a.b < 0
      TEST_ASSERT(checkTriclinicFilm(0.0, 0.0, 2.0, 1, true, flex));

      // Walls normal to b, with b.c < 0
      TEST_ASSERT(checkTriclinicFilm(-0.5*pi, 0.3, 0.5*pi, 1, true, flex));

      // Compatible cells
      TEST_ASSERT(!checkTriclinicFilm(0.0, 0.3, 0.5*pi, 1, true, flex));
      TEST_ASSERT(!checkTriclinicFilm(0.0, 0.0, 2.0, 2, true, flex));
   }

   /*
   * Flexible parameters of triclinic films. The angle between the two
   * lattice vectors parallel to the walls is flexible: gamma for walls
   * normal to c, and theta for walls normal to a (phi = pi/2) or to b
   * (phi = 0). The length of the vector normal to the walls is rigid
   * (no fBulk).
   */
   void testTriclinicFlexibleAngles()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmTriclinicFlexibleAngles.log");
      double pi = Constants::Pi;
      FSArray<bool,6> flex;

      // Walls normal to a: gamma = phi = pi/2, theta is free
      TEST_ASSERT(!checkTriclinicFilm(0.5*pi, 0.2, 0.5*pi, 0, false, flex));
      TEST_ASSERT(flex.size() == 6);
      TEST_ASSERT(!flex[0] && flex[1] && flex[2]);
      TEST_ASSERT(!flex[3] && flex[4] && !flex[5]);

      // Walls normal to b: gamma = pi/2, phi = 0, theta is free
      TEST_ASSERT(!checkTriclinicFilm(0.0, 0.2, 0.5*pi, 1, false, flex));
      TEST_ASSERT(flex[0] && !flex[1] && flex[2]);
      TEST_ASSERT(!flex[3] && flex[4] && !flex[5]);

      // Walls normal to c: theta = 0, gamma is free
      TEST_ASSERT(!checkTriclinicFilm(0.3, 0.0, 1.9, 2, false, flex));
      TEST_ASSERT(flex[0] && flex[1] && !flex[2]);
      TEST_ASSERT(!flex[3] && !flex[4] && flex[5]);
   }

   /*
   * A flexible theta is the angle between the two lattice vectors
   * parallel to the walls only for phi = pi/2 (walls normal to a) or
   * phi = 0 (walls normal to b). At theta = 0 every phi passes the
   * orthogonality check, so other values of phi must be rejected.
   */
   void testTriclinicPhiPin()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmTriclinicPhiPin.log");
      double pi = Constants::Pi;
      FSArray<bool,6> flex;

      // Walls normal to a, phi = 0
      TEST_ASSERT(checkTriclinicFilm(0.0, 0.0, 0.5*pi, 0, false, flex));

      // Walls normal to b, phi = pi/2
      TEST_ASSERT(checkTriclinicFilm(0.5*pi, 0.0, 0.5*pi, 1, false, flex));

      // The same cells are accepted if the angles are rigid
      TEST_ASSERT(!checkTriclinicFilm(0.0, 0.0, 0.5*pi, 0, true, flex));
      TEST_ASSERT(!checkTriclinicFilm(0.5*pi, 0.0, 0.5*pi, 1, true, flex));
   }

};

TEST_BEGIN(FilmMaskCheckTest)
TEST_ADD(FilmMaskCheckTest, testCheckLatticeVectorsSign)
TEST_ADD(FilmMaskCheckTest, testTriclinicFlexibleAngles)
TEST_ADD(FilmMaskCheckTest, testTriclinicPhiPin)
TEST_END(FilmMaskCheckTest)

#endif
