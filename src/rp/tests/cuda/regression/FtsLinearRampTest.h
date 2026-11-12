#ifndef RPG_FTS_LINEAR_RAMP_TEST_H
#define RPG_FTS_LINEAR_RAMP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/field/WFields.h>
#include <pscf/interaction/Interaction.h>

#include <util/misc/Exception.h>
#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): LinearRamp parameters and end values.
*/
class FtsLinearRampTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   template <int D>
   void initSystem(Rp::System<D,CUT>& system, std::string filename)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile(filename, in);
      system.readParam(in);
      in.close();
   }

};

TEST_BEGIN(FtsLinearRampTest)
TEST_ADD(FtsLinearRampTest, testIdsOutOfRange)
TEST_ADD(FtsLinearRampTest, testClosedPhiSum)
TEST_END(FtsLinearRampTest)

#endif
