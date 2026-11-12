#ifndef RPG_FTS_BD_SIMULATOR_TEST_H
#define RPG_FTS_BD_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/field/WFields.h>

#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): BdSimulator: ramps and final states.
*/
class FtsBdSimulatorTest : public LogFileUnitTest
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

   template <int D>
   void initSimulator(Rp::BdSimulator<D,CUT>& simulator,
                      std::string filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      simulator.readParam(in);
      in.close();
   }

   /*
   * A LinearRamp of the phi of a closed polymer in a mixed ensemble
   * (closed diblock, open homopolymer). The changes of phi need not
   * add up to zero.
   */
   void testPhiRampMixedEnsemble()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsPhiRampMixedEnsemble.log");

      Rp::System<1,CUT> system;
      initSystem(system, "in/Fts/param_system_mixed_1D");
      Rp::BdSimulator<1,CUT> simulator(system);
      bool success = true;
      try {
         initSimulator(simulator, "in/Fts/param_BdSimulator_phiRamp");
      } catch (Exception& e) {
         success = false;
      }
      TEST_ASSERT(success);
      TEST_ASSERT(simulator.hasRamp());

      system.w().readRGrid("in/Fts/w_mixed_1D.rf");
      simulator.simulate(4);
      TEST_ASSERT(simulator.iStep() == 4);
      double phi = system.mixture().polymer(0).phi();
      TEST_ASSERT(std::abs(phi - 0.7) < 1.0E-10);
   }

};

TEST_BEGIN(FtsBdSimulatorTest)
TEST_ADD(FtsBdSimulatorTest, testPhiRampMixedEnsemble)
TEST_END(FtsBdSimulatorTest)

#endif
