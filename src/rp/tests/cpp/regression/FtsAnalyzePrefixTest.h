#ifndef RPC_FTS_ANALYZE_PREFIX_TEST_H
#define RPC_FTS_ANALYZE_PREFIX_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/field/Domain.h>
#include <pscf/crystal/UnitCell.h>

#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): ANALYZE with an input prefix.
*/
class FtsAnalyzePrefixTest : public LogFileUnitTest
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
   * Copy a trajectory file, replacing the line after the "cell_param"
   * label of the header by the given cell parameter line.
   */
   void copyTrajectory(std::istream& in, std::ostream& out,
                       std::string cellLine)
   {
      std::string line;
      bool isCellLine = false;
      while (std::getline(in, line)) {
         if (isCellLine) {
            out << cellLine << "\n";
            isCellLine = false;
         } else {
            out << line << "\n";
            std::istringstream s(line);
            std::string label;
            s >> label;
            if (label == "cell_param") isCellLine = true;
         }
      }
   }

   /*
   * ANALYZE of a trajectory given by a name relative to the input
   * prefix. Return false if it throws an Exception, and set cell to the
   * unit cell parameter of the system after the analysis.
   */
   bool analyzeWithPrefix(std::string fileName, double& cell)
   {
      Rp::System<3,CPT> system;
      initSystem(system, "in/Fts/param_system_small");
      Rp::BdSimulator<3,CPT> simulator(system);
      initSimulator(simulator, "in/Fts/param_BdSimulator_analyze");
      system.fileMaster().setInputPrefix(filePrefix() + "out/");
      bool success = true;
      try {
         simulator.analyze(0, 3, "RGridTrajectoryReader", fileName);
      } catch (Exception& e) {
         success = false;
      }
      cell = system.domain().unitCell().parameter(0);
      return success;
   }

   /*
   * ANALYZE with an input prefix (as set by pscf_rpc -i), out/ here:
   *
   *  (1) a trajectory name without the prefix is found under the prefix;
   *  (2) with a decoy file of the same name in the working directory
   *      (other unit cell), the prefixed file is still the one read;
   *  (3) a name that already includes the prefix is still accepted.
   */
   void testAnalyzeInputPrefix()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsAnalyzeInputPrefix.log");
      const std::string name = "FtsAnalyzePrefix_traj.rf";

      // The trajectory (cell 3.0) at the prefixed path
      {
         std::ifstream in;
         openInputFile("in/Fts/traj_small.rf", in);
         std::ofstream out;
         openOutputFile("out/" + name, out);
         copyTrajectory(in, out, "      3.0000000000e+00");
      }

      // (1) No file of that name in the working directory
      double cell = 0.0;
      bool success = analyzeWithPrefix(name, cell);
      TEST_ASSERT(success);
      TEST_ASSERT(std::abs(cell - 3.0) < 1.0E-10);

      // (2) A decoy (cell 3.3) in the working directory, removed afterwards
      {
         std::ifstream in;
         openInputFile("in/Fts/traj_small.rf", in);
         std::ofstream out(name.c_str());
         TEST_ASSERT(out.is_open());
         copyTrajectory(in, out, "      3.3000000000e+00");
      }
      cell = 0.0;
      success = analyzeWithPrefix(name, cell);
      std::remove(name.c_str());
      TEST_ASSERT(success);
      TEST_ASSERT(std::abs(cell - 3.0) < 1.0E-10);

      // (3) The name with the prefix
      cell = 0.0;
      success = analyzeWithPrefix(filePrefix() + "out/" + name, cell);
      TEST_ASSERT(success);
      TEST_ASSERT(std::abs(cell - 3.0) < 1.0E-10);
   }

};

TEST_BEGIN(FtsAnalyzePrefixTest)
TEST_ADD(FtsAnalyzePrefixTest, testAnalyzeInputPrefix)
TEST_END(FtsAnalyzePrefixTest)

#endif
