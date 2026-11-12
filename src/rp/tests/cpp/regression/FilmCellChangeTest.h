#ifndef RPC_FILM_CELL_CHANGE_TEST_H
#define RPC_FILM_CELL_CHANGE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
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
* Regression tests for the film orthogonality check after a change of
* the unit cell.
*
* The film is 2D (oblique lattice, group p_2) with walls normal to lattice
* vector a, and w fields from in/Film/oblique_w.bf (cell 3.1, 1.6, pi/2).
* Lattice vector a = (a, 0) does not depend on the angle gamma, so a
* change of gamma leaves the vector normal to the walls unchanged. A cell
* with gamma != pi/2 must be rejected whether it is the first cell or a
* later one.
*/
class FilmCellChangeTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Parameter file of the film.
   *
   * If sweep is true, the file has a sweep of gamma from its initial
   * value down by 0.3, in 2 steps.
   */
   std::string filmParam(bool sweep = false)
   {
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
            << "    ds   0.01\n"
            << "  }\n"
            << "  Interaction{\n"
            << "    chi(\n"
            << "         1   0   15.0\n"
            << "    )\n"
            << "  }\n"
            << "  Domain{\n"
            << "    mesh           32 16\n"
            << "    lattice        oblique\n"
            << "    groupName      p_2\n"
            << "  }\n"
            << "  FilmEnvironment{\n"
            << "    normalVecId        0\n"
            << "    interfaceThickness 0.2\n"
            << "    excludedThickness  0.4\n"
            << "    chiBottom[  5.0  0.0  ]\n"
            << "    chiTop[  5.0  0.0  ]\n"
            << "  }\n"
            << "  AmIteratorBasis{\n"
            << "    epsilon      1.0e-7\n"
            << "    maxItr       300\n"
            << "    maxHist      30\n"
            << "    isFlexible   0\n"
            << "  }\n";
      if (sweep) {
         param << "  LinearSweep{\n"
               << "    ns            2\n"
               << "    baseFileName  out/FilmCellChangeSweep_\n"
               << "    nParameter    1\n"
               << "    parameters[\n"
               << "                  cell_param  2  -0.3\n"
               << "    ]\n"
               << "  }\n";
      }
      param << "}\n";
      return param.str();
   }

   /*
   * Read the parameter file and the w fields, and set the unit cell
   * (3.1, 1.6, gamma).
   */
   void initSystem(Rp::System<2,CPT>& system, double gamma,
                   bool sweep = false)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream in(filmParam(sweep));
      system.readParam(in);
      TEST_ASSERT(system.hasEnvironment());

      system.w().readBasis("in/Film/oblique_w.bf");
      FSArray<double, 6> parameters;
      parameters.append(3.1);
      parameters.append(1.6);
      parameters.append(gamma);
      system.setUnitCell(parameters);
   }

   /*
   * Execute a command script given as a string.
   *
   * Returns true if a Util::Exception was thrown.
   */
   bool readCommandsThrows(Rp::System<2,CPT>& system,
                           std::string const & commands)
   {
      std::istringstream in(commands);
      try {
         system.readCommands(in);
      } catch (Exception& e) {
         return true;
      }
      return false;
   }

   /*
   * SET_UNIT_CELL to gamma != pi/2 after a first COMPUTE at gamma = pi/2.
   */
   void testSetUnitCell()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmCellChangeSetUnitCell.log");
      double pi = Constants::Pi;

      // A first cell with gamma != pi/2 is rejected
      {
         Rp::System<2,CPT> system;
         initSystem(system, 1.2);
         TEST_ASSERT(readCommandsThrows(system, "COMPUTE\nFINISH\n"));
      }

      // A later cell with gamma != pi/2 is rejected, by COMPUTE ...
      {
         Rp::System<2,CPT> system;
         initSystem(system, 0.5*pi);
         TEST_ASSERT(!readCommandsThrows(system, "COMPUTE\nFINISH\n"));
         TEST_ASSERT(readCommandsThrows(system,
                     "SET_UNIT_CELL oblique 3.1 1.6 1.2\nCOMPUTE\nFINISH\n"));
      }

      // ... and by ITERATE
      {
         Rp::System<2,CPT> system;
         initSystem(system, 0.5*pi);
         TEST_ASSERT(!readCommandsThrows(system, "COMPUTE\nFINISH\n"));
         TEST_ASSERT(readCommandsThrows(system,
                     "SET_UNIT_CELL oblique 3.1 1.6 1.2\nITERATE\nFINISH\n"));
      }

      // A change of gamma back to pi/2 is accepted
      {
         Rp::System<2,CPT> system;
         initSystem(system, 0.5*pi);
         TEST_ASSERT(!readCommandsThrows(system, "COMPUTE\nFINISH\n"));
         FSArray<double, 6> parameters;
         parameters.append(3.1);
         parameters.append(1.6);
         parameters.append(1.2);
         system.setUnitCell(parameters);
         parameters[2] = 0.5*pi;
         system.setUnitCell(parameters);
         TEST_ASSERT(!readCommandsThrows(system, "COMPUTE\nFINISH\n"));
      }
   }

   /*
   * READ_W_BASIS of a file whose header has gamma != pi/2, after a first
   * COMPUTE at gamma = pi/2.
   */
   void testReadWHeader()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmCellChangeReadW.log");
      double pi = Constants::Pi;

      // Write a w field file with gamma = 1.2 in its header
      std::string wFile = "out/FilmCellChangeReadW_w.bf";
      {
         Rp::System<2,CPT> system;
         initSystem(system, 1.2);
         system.w().writeBasis(wFile);
      }

      Rp::System<2,CPT> system;
      initSystem(system, 0.5*pi);
      TEST_ASSERT(!readCommandsThrows(system, "COMPUTE\nFINISH\n"));
      TEST_ASSERT(readCommandsThrows(system,
                  "READ_W_BASIS " + wFile + "\nCOMPUTE\nFINISH\n"));
      TEST_ASSERT(eq(system.domain().unitCell().parameter(2), 1.2));
   }

   /*
   * Sweep of gamma that starts at pi/2.
   */
   void testSweep()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmCellChangeSweep.log");
      double pi = Constants::Pi;

      Rp::System<2,CPT> system;
      initSystem(system, 0.5*pi, true);
      TEST_ASSERT(system.hasSweep());

      // The first state (gamma = pi/2) is solved, the next one rejected
      bool thrown = false;
      try {
         system.sweep();
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(thrown);
      double gamma = system.domain().unitCell().parameter(2);
      TEST_ASSERT(gamma < 0.5*pi - 0.01);
      TEST_ASSERT(gamma > 0.5*pi - 0.3);
   }

};

TEST_BEGIN(FilmCellChangeTest)
TEST_ADD(FilmCellChangeTest, testSetUnitCell)
TEST_ADD(FilmCellChangeTest, testReadWHeader)
TEST_ADD(FilmCellChangeTest, testSweep)
TEST_END(FilmCellChangeTest)

#endif
