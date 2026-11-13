#ifndef RPG_RAMP_PARAMETER_ROUND_TRIP_TEST_H
#define RPG_RAMP_PARAMETER_ROUND_TRIP_TEST_H

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
class RampParameterRoundTripTest : public LogFileUnitTest
{
public:
   void setUp() { setVerbose(0); }

   void testRampSerialization()
   {
      printMethod(TEST_FUNC);
      char const * parameters[] = {"v_monomer 0.1", "chi 0 1 0.1"};
      for (char const * parameter : parameters) {
         Rp::RampParameter<1,CUT> original, copy;
         std::istringstream in(parameter);
         in >> original;
         std::ostringstream out;
         out << original;
         std::istringstream reread(out.str());
         reread >> copy;
         TEST_ASSERT(copy.type() == original.type());
         TEST_ASSERT(std::abs(copy.change()-original.change()) < 1e-15);
      }
   }
};

TEST_BEGIN(RampParameterRoundTripTest)
TEST_ADD(RampParameterRoundTripTest, testRampSerialization)
TEST_END(RampParameterRoundTripTest)

#endif
