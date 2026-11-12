#ifndef RPC_FTS_AM_COMPRESSOR_LOG_TEST_H
#define RPC_FTS_AM_COMPRESSOR_LOG_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/compressor/AmCompressor.h>
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
* Regression tests (pscfpp fixes): AmCompressor with verbose 0 writes
* nothing to the log. It wrote " Entering setup", " Exiting setup",
* " Entering compress" and " Exiting compress" on every call.
*/
class FtsAmCompressorLogTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Contents of the log file written so far.
   */
   std::string logContents(std::string const & filename)
   {
      logFile_.flush();
      std::ifstream file;
      openInputFile(filename, file);
      std::stringstream buffer;
      buffer << file.rdbuf();
      return buffer.str();
   }

   void testCompressLog()
   {
      printMethod(TEST_FUNC);
      std::string const logName = "out/FtsAmCompressorLog.log";
      openLogFile(logName.c_str());

      Rp::System<3,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile("in/Fts/param_system_small", in);
      system.readParam(in);
      in.close();

      // AmCompressor with verbose 0
      Rp::AmCompressor<3,CPT> compressor(system);
      openInputFile("in/Fts/param_AmCompressor_relNorm", in);
      compressor.readParam(in);
      in.close();
      system.w().readRGrid("in/Fts/w_small.rf");

      // Two calls of compress (each calls setup)
      std::string before = logContents(logName);
      TEST_ASSERT(compressor.compress() == 0);
      TEST_ASSERT(compressor.compress() == 0);
      std::string after = logContents(logName);
      TEST_ASSERT(after.find("Entering") == std::string::npos);
      TEST_ASSERT(after.find("Exiting") == std::string::npos);
      TEST_ASSERT(after == before);
   }

};

TEST_BEGIN(FtsAmCompressorLogTest)
TEST_ADD(FtsAmCompressorLogTest, testCompressLog)
TEST_END(FtsAmCompressorLogTest)

#endif
