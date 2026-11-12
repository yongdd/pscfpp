#ifndef RPG_COMMAND_TEST_H
#define RPG_COMMAND_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <rp/field/Mask.h>
#include <rp/solvers/Mixture.h>
#include <rp/environment/Environment.h>

#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>

using namespace Util;
using namespace Pscf;

/*
* Regression tests of the command script loop, System::readCommands.
*/
class CommandTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * ANALYZE in a system without a Simulator throws an Exception.
   */
   void testAnalyzeNoSimulator()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testAnalyzeNoSimulator.log");

      Rp::System<1,CUT> system;
      initSystem(system, "in/Command/lam_param", "in/Command/lam_w.bf");
      TEST_ASSERT(!system.hasSimulator());

      TEST_ASSERT(readCommandsThrows(system,
                  "ANALYZE 0 10 HamiltonianAnalyzer out/none\nFINISH\n"));
   }

   /*
   * COMPUTE in a film system without an Iterator.
   */
   void testFilmNoIterator()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testFilmNoIterator.log");

      Rp::System<1,CUT> system;
      initSystem(system, "in/Command/film_param", "in/Command/film_w.bf");
      TEST_ASSERT(system.hasEnvironment());
      TEST_ASSERT(!system.hasIterator());

      // The FilmEnvironment generates the mask and the wall fields
      TEST_ASSERT(!readCommandsThrows(system, "COMPUTE\nFINISH\n"));
      TEST_ASSERT(system.mask().hasData());
      TEST_ASSERT(system.h().hasData());
      TEST_ASSERT(system.c().hasData());
   }

   /*
   * WRITE_STRESS in a film system without an Iterator.
   */
   void testFilmNoIteratorStress()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testFilmNoIteratorStress.log");

      Rp::System<1,CUT> system;
      initSystem(system, "in/Command/film_param", "in/Command/film_w.bf");
      TEST_ASSERT(!system.hasIterator());

      std::string stressFile = "out/testFilmNoIteratorStress_stress";
      removeFile(stressFile);
      TEST_ASSERT(!readCommandsThrows(system,
                  "COMPUTE\nWRITE_STRESS " + stressFile + "\nFINISH\n"));
      TEST_ASSERT(system.mixture().hasStress());
      TEST_ASSERT(system.environment().hasStress());
      TEST_ASSERT(fileExists(stressFile));
   }

   // Utility functions

   /*
   * Read a parameter file and, if a file name is given, w fields in
   * basis format.
   */
   template <int D>
   void initSystem(Rp::System<D,CUT>& system,
                   std::string const & paramFileName,
                   std::string const & wFileName = "")
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile(paramFileName, in);
      system.readParam(in);
      in.close();

      if (!wFileName.empty()) {
         system.w().readBasis(wFileName);
      }
   }

   /*
   * Execute a command script given as a string, as from a command file.
   *
   * Returns true if a Util::Exception was thrown.
   */
   template <int D>
   bool readCommandsThrows(Rp::System<D,CUT>& system,
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

TEST_BEGIN(CommandTest)
TEST_ADD(CommandTest, testAnalyzeNoSimulator)
TEST_ADD(CommandTest, testFilmNoIterator)
TEST_ADD(CommandTest, testFilmNoIteratorStress)
TEST_END(CommandTest)

#endif
