#ifndef RPC_FTS_ANALYZE_PREFIX_LOG_TEST_H
#define RPC_FTS_ANALYZE_PREFIX_LOG_TEST_H

#include <test/UnitTestRunner.h>
#include <rp/system/System.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>
#include <fstream>
#include <sstream>
#include <string>

using namespace Util;
using namespace Pscf;

/*
* Regression test (pscfpp fixes): ANALYZE with an input prefix writes
* no error message to the log when the trajectory file is found.
*/
class FtsAnalyzePrefixLogTest : public LogFileUnitTest
{

public:

   /*
   * ANALYZE trajectory fileName with the input prefix out/, with the
   * log written to logName. Return true if the trajectory is read and
   * the log holds no echoed Exception message.
   */
   bool analyzeQuietly(std::string fileName, std::string logName)
   {
      openLogFile(logName);
      bool success = true;
      {
         Rp::System<3,CPT> system;
         system.fileMaster().setInputPrefix(filePrefix());
         system.fileMaster().setOutputPrefix(filePrefix());
         std::ifstream in;
         openInputFile("in/Fts/param_system_small", in);
         system.readParam(in);
         in.close();
         Rp::BdSimulator<3,CPT> simulator(system);
         openInputFile("in/Fts/param_BdSimulator_analyze", in);
         simulator.readParam(in);
         in.close();
         system.fileMaster().setInputPrefix(filePrefix() + "out/");
         try {
            simulator.analyze(0, 3, "RGridTrajectoryReader", fileName);
         } catch (Exception& e) {
            success = false;
         }
      }
      closeLogFile();
      std::ifstream log;
      openInputFile(logName, log);
      std::stringstream text;
      text << log.rdbuf();
      return success
             && text.str().find("Read Header") != std::string::npos
             && text.str().find("Message :") == std::string::npos;
   }

   /*
   * A name relative to the prefix and a name that already includes it
   * (the earlier usage) both read the trajectory. The fallback to the
   * second form logged the failed attempt at the doubled prefix.
   */
   void testAnalyzePrefixLog()
   {
      printMethod(TEST_FUNC);
      std::string name = "FtsAnalyzePrefixLog_traj.rf";
      {
         std::ifstream in;
         openInputFile("in/Fts/traj_small.rf", in);
         std::ofstream out;
         openOutputFile("out/" + name, out);
         out << in.rdbuf();
      }
      TEST_ASSERT(analyzeQuietly(name, "out/FtsAnalyzePrefixLog_1.log"));
      TEST_ASSERT(analyzeQuietly(filePrefix() + "out/" + name,
                                 "out/FtsAnalyzePrefixLog_2.log"));
   }

};

TEST_BEGIN(FtsAnalyzePrefixLogTest)
TEST_ADD(FtsAnalyzePrefixLogTest, testAnalyzePrefixLog)
TEST_END(FtsAnalyzePrefixLogTest)

#endif
