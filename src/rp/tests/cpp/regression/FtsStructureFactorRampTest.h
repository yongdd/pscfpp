#ifndef RPC_FTS_STRUCTURE_FACTOR_RAMP_TEST_H
#define RPC_FTS_STRUCTURE_FACTOR_RAMP_TEST_H

#include <test/UnitTestRunner.h>
#include <rp/system/System.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>
#include <fstream>
#include <sstream>
#include <string>

using namespace Util;
using namespace Pscf;

/*
* Regression test (pscfpp fixes): SIMULATE with BinaryStructureFactor
* and a ramp that changes the unit cell.
*/
class FtsStructureFactorRampTest : public LogFileUnitTest
{

public:

   /*
   * SIMULATE 4 steps (cubic cell 3.0) with simulator file paramFile,
   * in which @PARAMETER@ is replaced by the ramp parameter. Return the
   * Exception message, or "" if none. A rejected run must not have
   * taken a step or changed the cell.
   */
   template <class SimulatorT>
   std::string simulate(std::string paramFile, std::string parameter)
   {
      Rp::System<3,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile("in/Fts/param_system_small", in);
      system.readParam(in);
      in.close();
      SimulatorT simulator(system);
      openInputFile(paramFile, in);
      std::stringstream text;
      text << in.rdbuf();
      std::string s = text.str();
      s.replace(s.find("@PARAMETER@"), 11, parameter);
      std::istringstream param(s);
      simulator.readParam(param);
      system.w().readRGrid("in/Fts/w_small.rf");
      std::string message;
      try {
         simulator.simulate(4);
      } catch (Exception& e) {
         message = e.message();
         TEST_ASSERT(simulator.iStep() == 0);
         TEST_ASSERT(system.domain().unitCell().parameter(0) == 3.0);
      }
      return message;
   }

   /*
   * A cell_param ramp is rejected before the first step, with BD and
   * MC, naming the analyzer (the run threw at the second sample).
   * Control: a chi ramp runs.
   */
   void testCellRampRejected()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsStructureFactorRamp.log");
      std::string bd, mc, chi;
      bd = simulate< Rp::BdSimulator<3,CPT> >("in/Fts/sfRamp_bd",
                                              "cell_param 0 0.4");
      mc = simulate< Rp::McSimulator<3,CPT> >("in/Fts/sfRamp_mc",
                                              "cell_param 0 0.4");
      chi = simulate< Rp::BdSimulator<3,CPT> >("in/Fts/sfRamp_bd",
                                               "chi 0 1 0.5");
      std::string const name = "Analyzer BinaryStructureFactor requires";
      TEST_ASSERT(bd.find(name) != std::string::npos);
      TEST_ASSERT(mc.find(name) != std::string::npos);
      TEST_ASSERT(chi == "");
   }

};

TEST_BEGIN(FtsStructureFactorRampTest)
TEST_ADD(FtsStructureFactorRampTest, testCellRampRejected)
TEST_END(FtsStructureFactorRampTest)

#endif
