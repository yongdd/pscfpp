#ifndef RPC_SCFT_PARAM_FILE_END_TEST_H
#define RPC_SCFT_PARAM_FILE_END_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>

#include <util/misc/Exception.h>
#include <util/misc/FileMaster.h>
#include <util/param/Label.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <string>

using namespace Util;
using namespace Pscf;

/*
* SCFT parameter files whose closing bracket is the last character of
* the file, without a final newline. System reads the optional Sweep
* and Simulator blocks after the Iterator: the attempt to read them
* left the closing '}' of System in the label buffer at end-of-file,
* and reading that required '}' then stopped with 'EOF before reading
* required label', with or without a Sweep block.
*/
class ScftParamFileEndTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * AmIterator, no Sweep.
   */
   void testNoFinalNewline()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftParamFileEndNoSweep.log");

      Rp::System<1,CPT> system;
      TEST_ASSERT(lastCharacter("in/Scft/eof_param") == '}');
      TEST_ASSERT(!readParamThrows(system, "in/Scft/eof_param"));
      TEST_ASSERT(system.hasIterator());
      TEST_ASSERT(!system.hasSweep());
      TEST_ASSERT(!system.hasSimulator());
      TEST_ASSERT(system.mixture().nPolymer() == 1);
   }

   /*
   * AmIterator and LinearSweep.
   */
   void testNoFinalNewlineSweep()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftParamFileEndSweep.log");

      Rp::System<1,CPT> system;
      TEST_ASSERT(lastCharacter("in/Scft/eof_sweep_param") == '}');
      TEST_ASSERT(!readParamThrows(system, "in/Scft/eof_sweep_param"));
      TEST_ASSERT(system.hasIterator());
      TEST_ASSERT(system.hasSweep());
      TEST_ASSERT(!system.hasSimulator());
   }

private:

   /*
   * Read the parameter file, return true if this throws an Exception.
   */
   bool readParamThrows(Rp::System<1,CPT>& system, std::string fileName)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      bool thrown = false;
      try {
         std::ifstream in;
         openInputFile(fileName, in);
         system.readParam(in);
         in.close();
      } catch (Exception& e) {
         thrown = true;
         Label::clear(); // Discard any label left by the failed read
      }
      return thrown;
   }

   /*
   * Return the last character of a file, to check that the input
   * still ends without a newline.
   */
   char lastCharacter(std::string fileName)
   {
      std::ifstream in;
      openInputFile(fileName, in);
      in.seekg(-1, std::ios::end);
      char c = '\0';
      in.get(c);
      in.close();
      return c;
   }

};

TEST_BEGIN(ScftParamFileEndTest)
TEST_ADD(ScftParamFileEndTest, testNoFinalNewline)
TEST_ADD(ScftParamFileEndTest, testNoFinalNewlineSweep)
TEST_END(ScftParamFileEndTest)

#endif
