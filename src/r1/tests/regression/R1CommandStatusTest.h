#ifndef R1_COMMAND_STATUS_TEST_H
#define R1_COMMAND_STATUS_TEST_H

#include "R1TestBase.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

/*
* Regression tests (pscfpp fixes): the status returned by the command
* loop, which is the exit status of pscf_1d.
*/
class R1CommandStatusTest : public R1TestBase
{

public:

   /*
   * Status 0 if every command was executed.
   */
   void testSuccess()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1CommandStatusTestSuccess.log");

      System sys;
      readSystem(sys, "in/R1/am_lam.prm");
      TEST_ASSERT(readCommands(sys,
                  "READ_W in/R1/am_lam.w\nITERATE\nFINISH\n") == 0);

      // End of the script without FINISH
      TEST_ASSERT(readCommands(sys,
                  "READ_W in/R1/am_lam.w\nCOMPUTE\n") == 0);
   }

   /*
   * An ITERATE command that does not converge stops the command loop
   * with status 1. The exit status of pscf_1d was 0.
   */
   void testIterateFailure()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1CommandStatusTestIterateFailure.log");

      // The iterator stops after 1 iteration, far from the solution
      System sys;
      readSystem(sys, "in/R1/am_lam.prm");
      std::string wFile = "out/R1CommandStatusTestIterateFailure.w";
      std::remove((filePrefix() + wFile).c_str());
      TEST_ASSERT(readCommands(sys, "READ_W in/R1/planar.w\nITERATE\n"
                               "WRITE_W " + wFile + "\nFINISH\n") == 1);

      // The command after ITERATE was not executed
      std::ifstream file((filePrefix() + wFile).c_str());
      TEST_ASSERT(!file.is_open());
   }

   /*
   * An unknown command stops the command loop with status 1.
   */
   void testUnknownCommand()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1CommandStatusTestUnknownCommand.log");

      System sys;
      readSystem(sys, "in/R1/am_lam.prm");
      TEST_ASSERT(readCommands(sys,
                  "READ_W in/R1/am_lam.w\nNOT_A_COMMAND\nFINISH\n") == 1);
   }

private:

   /*
   * Execute the commands in a string, return the status.
   */
   int readCommands(System& sys, std::string const & commands)
   {
      std::istringstream in(commands);
      return sys.readCommands(in);
   }

};

TEST_BEGIN(R1CommandStatusTest)
TEST_ADD(R1CommandStatusTest, testSuccess)
TEST_ADD(R1CommandStatusTest, testIterateFailure)
TEST_ADD(R1CommandStatusTest, testUnknownCommand)
TEST_END(R1CommandStatusTest)

#endif
