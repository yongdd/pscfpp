#ifndef RPG_COMMAND_STATUS_TEST_H
#define RPG_COMMAND_STATUS_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/field/WFields.h>

#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): the status returned by the command
* loop, which is the exit status of pscf_rpc and pscf_rpg.
*/
class CommandStatusTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Status 0 if every command was executed.
   */
   void testSuccess()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testCommandStatusSuccess.log");

      Rp::System<1,CUT> system;
      initSystem(system, "in/Command/lam_param");
      TEST_ASSERT(readCommands(system,
                  "READ_W_BASIS in/Command/lam_w.bf\nITERATE\nFINISH\n")
                  == 0);

      // End of the script without FINISH
      TEST_ASSERT(readCommands(system,
                  "READ_W_BASIS in/Command/lam_w.bf\nCOMPUTE\n") == 0);
   }

   /*
   * An ITERATE command that does not converge stops the command loop
   * with status 1. The exit status of the program was 0.
   */
   void testIterateFailure()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testCommandStatusIterateFailure.log");

      // The iterator stops after 1 iteration
      Rp::System<1,CUT> system;
      initSystem(system, "in/Command/status_lam_param");
      std::string wFile = "out/testCommandStatusIterateFailure_w.bf";
      removeFile(wFile);
      TEST_ASSERT(readCommands(system,
                  "READ_W_BASIS in/Command/lam_w.bf\nITERATE\n"
                  "WRITE_W_BASIS " + wFile + "\nFINISH\n") == 1);
      TEST_ASSERT(!fileExists(wFile));
   }

   /*
   * A COMPRESS command that does not converge stops the command loop
   * with status 1.
   */
   void testCompressFailure()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testCommandStatusCompressFailure.log");

      // The compressor of this system stops after 2 iterations
      Rp::System<1,CUT> system;
      initSystem(system, "in/Command/fts_param_fail");
      std::string wFile = "out/testCommandStatusCompressFailure_w.rf";
      removeFile(wFile);
      TEST_ASSERT(readCommands(system,
                  "READ_W_RGRID in/Command/fts_w.rf\nCOMPRESS\n"
                  "WRITE_W_RGRID " + wFile + "\nFINISH\n") == 1);
      TEST_ASSERT(!fileExists(wFile));
   }

   /*
   * An unknown command stops the command loop with status 1.
   */
   void testUnknownCommand()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testCommandStatusUnknownCommand.log");

      Rp::System<1,CUT> system;
      initSystem(system, "in/Command/lam_param");
      TEST_ASSERT(readCommands(system,
                  "READ_W_BASIS in/Command/lam_w.bf\nNOT_A_COMMAND\n"
                  "FINISH\n") == 1);
   }

private:

   /*
   * Read a parameter file.
   */
   template <int D>
   void initSystem(Rp::System<D,CUT>& system,
                   std::string const & paramFileName)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile(paramFileName, in);
      system.readParam(in);
      in.close();
   }

   /*
   * Execute a command script given as a string, return the status.
   */
   template <int D>
   int readCommands(Rp::System<D,CUT>& system,
                    std::string const & commands)
   {
      std::istringstream in(commands);
      return system.readCommands(in);
   }

   /*
   * Does a file exist? (name relative to the test directory)
   */
   bool fileExists(std::string const & filename)
   {
      std::ifstream file((filePrefix() + filename).c_str());
      return file.is_open();
   }

   /*
   * Remove an output file left by an earlier run.
   */
   void removeFile(std::string const & filename)
   {  std::remove((filePrefix() + filename).c_str()); }

};

TEST_BEGIN(CommandStatusTest)
TEST_ADD(CommandStatusTest, testSuccess)
TEST_ADD(CommandStatusTest, testIterateFailure)
TEST_ADD(CommandStatusTest, testCompressFailure)
TEST_ADD(CommandStatusTest, testUnknownCommand)
TEST_END(CommandStatusTest)

#endif
