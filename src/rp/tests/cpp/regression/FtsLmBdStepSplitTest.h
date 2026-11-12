#ifndef RPC_FTS_LM_BD_STEP_SPLIT_TEST_H
#define RPC_FTS_LM_BD_STEP_SPLIT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/field/WFields.h>
#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldComparison.h>

#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <string>
#include <vector>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): LMBdStep keeps its noise chain
* across SIMULATE commands.
*/
class FtsLmBdStepSplitTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Run BD simulations of nSteps[0], nSteps[1], ... steps in turn
   * (one SIMULATE command each), from the same initial state and
   * seed, and return the final w fields.
   */
   void runBd(std::vector<int> const & nSteps,
              DArray< RField<3,CPT> >& w)
   {
      Rp::System<3,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile("in/Fts/param_system_small", in);
      system.readParam(in);
      in.close();

      Rp::BdSimulator<3,CPT> simulator(system);
      openInputFile("in/Fts/param_BdSimulator_split", in);
      simulator.readParam(in);
      in.close();

      system.w().readRGrid("in/Fts/w_small.rf");
      for (unsigned int i = 0; i < nSteps.size(); ++i) {
         simulator.simulate(nSteps[i]);
      }

      int nMonomer = system.w().rgrid().capacity();
      w.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w[i].allocate(system.w().rgrid(i).meshDimensions());
         w[i] = system.w().rgrid(i);
      }
   }

   /*
   * A run split into SIMULATE 2, SIMULATE 1 and SIMULATE 3 must give
   * the same trajectory as SIMULATE 6, bitwise: the LM step of each
   * command must continue with the old noise of the previous step.
   */
   void testSplitRun()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsLmBdStepSplitRun.log");

      DArray< RField<3,CPT> > w1, w2;
      std::vector<int> nSteps(1, 6);
      runBd(nSteps, w1);
      nSteps.clear();
      nSteps.push_back(2);
      nSteps.push_back(1);
      nSteps.push_back(3);
      runBd(nSteps, w2);

      RFieldComparison<3,CPT> comparison;
      comparison.compare(w1, w2);
      Log::file() << "Max difference, split vs one run = "
                  << comparison.maxDiff() << std::endl;
      TEST_ASSERT(comparison.maxDiff() == 0.0);
   }

};

TEST_BEGIN(FtsLmBdStepSplitTest)
TEST_ADD(FtsLmBdStepSplitTest, testSplitRun)
TEST_END(FtsLmBdStepSplitTest)

#endif
