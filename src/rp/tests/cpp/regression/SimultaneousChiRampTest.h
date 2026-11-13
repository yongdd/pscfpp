#ifndef RPC_SIMULTANEOUS_CHI_RAMP_TEST_H
#define RPC_SIMULTANEOUS_CHI_RAMP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <rp/system/System.h>
#include <rp/scft/ScftThermo.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/ramp/LinearRamp.h>
#include <rp/fts/ramp/RampParameter.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <pscf/field/RField.h>
#include <pscf/backend/VecOp.h>
#include <pscf/interaction/Interaction.h>
#include <util/tests/LogFileUnitTest.h>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace Util;
using namespace Pscf;

// Exact homogeneous thermodynamics and simultaneous parameter updates.
class SimultaneousChiRampTest : public LogFileUnitTest
{
public:
   void setUp() { setVerbose(0); }

   void testSimultaneousChiRamp()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SimultaneousChiRampNumericalChiRamp.log");
      Rp::System<1,CPT> system;
      std::istringstream input(
         "System{ Mixture{ nMonomer 3 monomers[ 1 1 1 ] nPolymer 3 "
         "Polymer{ type linear nBlock 1 blocks[ 0 1 ] phi .3333333333333333 } "
         "Polymer{ type linear nBlock 1 blocks[ 1 1 ] phi .3333333333333333 } "
         "Polymer{ type linear nBlock 1 blocks[ 2 1 ] phi .3333333333333333 } "
         "ds .05 } Interaction{ chi( 0 1 8 0 2 8 1 2 8 ) } "
         "Domain{ mesh 8 lattice lamellar } }");
      system.readParam(input);
      Rp::BdSimulator<1,CPT> simulator(system);
      std::istringstream simInput("BdSimulator{ seed 1234 }");
      simulator.readParam(simInput);
      simulator.analyzeChi();
      Rp::LinearRamp<1,CPT> ramp(simulator);
      std::istringstream rampInput(
         "LinearRamp{ nParameter 3 parameters[ "
         "chi 0 1 72 chi 0 2 72 chi 1 2 72 ] }");
      ramp.readParam(rampInput);
      ramp.setup(1);
      ramp.setParameters(1);
      for (int i = 0; i < 3; ++i) {
         for (int j = i+1; j < 3; ++j) {
            TEST_ASSERT(std::abs(system.interaction().chi(i,j)-80.0) < 1e-12);
         }
      }
      for (int i = 0; i < 2; ++i) {
         TEST_ASSERT(std::abs(simulator.chiEval(i)+80.0) < 1e-10);
      }
   }

};

TEST_BEGIN(SimultaneousChiRampTest)
TEST_ADD(SimultaneousChiRampTest, testSimultaneousChiRamp)
TEST_END(SimultaneousChiRampTest)

#endif
