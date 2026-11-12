#ifndef RPC_FTS_ANALYZER_TEST_H
#define RPC_FTS_ANALYZER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/analyzer/AnalyzerManager.h>
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
* Regression tests (pscfpp fixes): Analyzers, trajectory files and ANALYZE.
*/
class FtsAnalyzerTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   template <int D>
   void initSystem(Rp::System<D,CPT>& system, std::string filename)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile(filename, in);
      system.readParam(in);
      in.close();
   }

   template <int D>
   void initSimulator(Rp::BdSimulator<D,CPT>& simulator,
                      std::string filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      simulator.readParam(in);
      in.close();
   }

   /*
   * Return the number of nonempty lines in a file.
   */
   int countLines(std::string filename)
   {
      std::ifstream file;
      openInputFile(filename, file);
      std::string line;
      int n = 0;
      while (std::getline(file, line)) {
         if (line.find_first_not_of(" \t\r") != std::string::npos) {
            ++n;
         }
      }
      file.close();
      return n;
   }

   /*
   * Two SIMULATE runs with averaging analyzers that write .dat files.
   *
   * The .dat files of an AverageListAnalyzer (HamiltonianAnalyzer) and
   * an AverageAnalyzer (MaxOrderParameter) must hold the samples of
   * both runs.
   */
   void testSimulateTwice()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsSimulateTwice.log");

      Rp::System<3,CPT> system;
      initSystem(system, "in/Fts/param_system_small");
      Rp::BdSimulator<3,CPT> simulator(system);
      initSimulator(simulator, "in/Fts/param_BdSimulator_twice");
      system.w().readRGrid("in/Fts/w_small.rf");

      bool success = true;
      try {
         simulator.simulate(3);
         simulator.simulate(3);
      } catch (Exception& e) {
         success = false;
      }
      TEST_ASSERT(success);

      // Steps 0 to 3 of each run, one line per sample
      TEST_ASSERT(countLines("out/FtsTwice_hamiltonian.dat") == 8);
      TEST_ASSERT(countLines("out/FtsTwice_maxOrder.dat") == 8);
   }

};

TEST_BEGIN(FtsAnalyzerTest)
TEST_ADD(FtsAnalyzerTest, testSimulateTwice)
TEST_END(FtsAnalyzerTest)

#endif
