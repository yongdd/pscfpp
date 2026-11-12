#ifndef RPG_FTS_BD_MOVE_N_STEP_TEST_H
#define RPG_FTS_BD_MOVE_N_STEP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/fts/montecarlo/McMove.h>
#include <rp/fts/montecarlo/McMoveManager.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>

#include <util/misc/Exception.h>
#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): the nStep parameter of BdMove is a
* positive integer. It was read as a double, so that nStep 2.5 gave
* moves of 3 BD steps and nStep 0 moves without BD steps.
*/
class FtsBdMoveNStepTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   void initSystem(Rp::System<3,CUT>& system)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile("in/Fts/param_system_disordered", in);
      system.readParam(in);
      in.close();
   }

   /*
   * McSimulator parameter block with one BdMove of the given nStep.
   */
   std::string simulatorParam(std::string const & nStep)
   {
      return "McSimulator{\n"
             "   seed              1784106389\n"
             "   McMoveManager{\n"
             "      BdMove{\n"
             "         probability   1.0\n"
             "         mobility      1.0e-4\n"
             "         nStep         " + nStep + "\n"
             "      }\n"
             "   }\n"
             "   LrAmCompressor{\n"
             "      epsilon      1.0e-6\n"
             "      maxItr       200\n"
             "      maxHist      30\n"
             "      verbose      0\n"
             "      errorType    rms\n"
             "   }\n"
             "}\n";
   }

   /*
   * Does reading a McSimulator with this nStep throw an Exception?
   */
   bool readThrows(std::string const & nStep)
   {
      Rp::System<3,CUT> system;
      initSystem(system);
      Rp::McSimulator<3,CUT> simulator(system);
      std::istringstream in(simulatorParam(nStep));
      bool thrown = false;
      try {
         simulator.readParam(in);
      } catch (Exception& e) {
         thrown = true;
      }
      if (verbose() > 0) {
         std::cout << "\nnStep " << nStep << ": thrown = " << thrown;
      }
      return thrown;
   }

   /*
   * Values that are not positive integers are rejected.
   */
   void testInvalidNStep()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsBdMoveInvalidNStep.log");

      TEST_ASSERT(readThrows("2.5"));
      TEST_ASSERT(readThrows("1e2"));
      TEST_ASSERT(readThrows("0"));
      TEST_ASSERT(readThrows("-3"));
      TEST_ASSERT(!readThrows("3"));
   }

   /*
   * A move with nStep 3 takes 3 BD steps.
   */
   void testNStep()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsBdMoveNStep.log");

      Rp::System<3,CUT> system;
      initSystem(system);
      Rp::McSimulator<3,CUT> simulator(system);
      std::istringstream in(simulatorParam("3"));
      simulator.readParam(in);
      system.w().readRGrid("in/Fts/w_dis.rf");
      simulator.simulate(2);

      Rp::McMove<3,CUT>& move = simulator.mcMoveManager()[0];
      TEST_ASSERT(move.nAttempt() == 2);
      TEST_ASSERT(move.nAccept() == 2);

      // Line written by the move at the end of the simulation
      logFile_.flush();
      std::ifstream file;
      openInputFile("out/FtsBdMoveNStep.log", file);
      std::stringstream buffer;
      buffer << file.rdbuf();
      bool found = (buffer.str().find("BdMove: 6 BD steps, 0 failed")
                    != std::string::npos);
      TEST_ASSERT(found);
   }

};

TEST_BEGIN(FtsBdMoveNStepTest)
TEST_ADD(FtsBdMoveNStepTest, testInvalidNStep)
TEST_ADD(FtsBdMoveNStepTest, testNStep)
TEST_END(FtsBdMoveNStepTest)

#endif
