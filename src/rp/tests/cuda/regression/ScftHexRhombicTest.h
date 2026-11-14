#ifndef RPG_SCFT_HEX_RHOMBIC_TEST_H
#define RPG_SCFT_HEX_RHOMBIC_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/scft/ScftThermo.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>

#include <util/containers/FSArray.h>
#include <util/format/Dbl.h>
#include <util/math/Constants.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Hexagonal state in a rhombic cell without a space group, 12x12 mesh.
*
* kSq of wavevectors with a Nyquist index was averaged over the signs of
* the Nyquist indices for every cell parameter of a variable-angle
* lattice. This average is not invariant under the 6-fold rotation of a
* rhombic cell at 120 degrees, so that the converged hexagonal state of
* the hexagonal lattice (diblock, f = 0.3, chi = 20) had a different
* free energy in the rhombic cell (by 1.4e-4) and a flexible rhombic
* cell relaxed to 7e-4 rad from 120 degrees, with a free energy 8.7e-5
* above the hexagonal value. The average is now used only near right
* angles: the free energy at 120 degrees equals the hexagonal value, and
* the flexible cell relaxes to within 1e-4 rad of 120 degrees (where
* the minimum image itself has a kink), with the hexagonal free energy.
*/
class ScftHexRhombicTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   void setupSystem(Rp::System<2,CUT>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

   /*
   * Free energy of the hexagonal state in the hexagonal lattice.
   */
   double hexagonalFreeEnergy()
   {
      Rp::System<2,CUT> system;
      setupSystem(system, "in/Scft/hexrh_hex12_param");
      system.w().readRGrid("in/Scft/hexrh_hex12_w.rf");
      system.compute();
      system.scft().compute();
      return system.scft().fHelmholtz();
   }

   /*
   * Rigid rhombic cell at (a, 120 degrees): same free energy.
   */
   void testRigid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftHexRhombicRigid.log");
      double fHex = hexagonalFreeEnergy();

      Rp::System<2,CUT> system;
      setupSystem(system, "in/Scft/hexrh_rhombic12_param");
      system.w().readRGrid("in/Scft/hexrh_rhombic12_w.rf");
      system.compute();
      system.scft().compute();
      double f = system.scft().fHelmholtz();
      Log::file() << "f hexagonal " << Dbl(fHex, 20, 13)
                  << "  f rhombic " << Dbl(f, 20, 13) << std::endl;
      TEST_ASSERT(std::abs(f - fHex) < 1.0E-10);
   }

   /*
   * Flexible rhombic cell started at (1.01 a, 120 degrees + 0.02).
   */
   void testFlexible()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftHexRhombicFlexible.log");
      double fHex = hexagonalFreeEnergy();

      Rp::System<2,CUT> system;
      setupSystem(system, "in/Scft/hexrh_rhombic12_param");
      system.w().readRGrid("in/Scft/hexrh_rhombic12_w.rf");
      double third = 2.0*Constants::Pi/3.0;
      FSArray<double, 6> p = system.domain().unitCell().parameters();
      p[0] *= 1.01;
      p[1] = third + 0.02;
      system.setUnitCell(p);
      TEST_ASSERT(system.iterate() == 0);
      double gamma = system.domain().unitCell().parameter(1);
      double f = system.scft().fHelmholtz();
      Log::file() << "gamma - 2pi/3 " << Dbl(gamma - third, 20, 13)
                  << "  f - f hexagonal " << Dbl(f - fHex, 20, 13)
                  << std::endl;
      TEST_ASSERT(std::abs(gamma - third) < 1.0E-4);
      TEST_ASSERT(std::abs(f - fHex) < 1.0E-8);
   }

};

TEST_BEGIN(ScftHexRhombicTest)
TEST_ADD(ScftHexRhombicTest, testRigid)
TEST_ADD(ScftHexRhombicTest, testFlexible)
TEST_END(ScftHexRhombicTest)

#endif
