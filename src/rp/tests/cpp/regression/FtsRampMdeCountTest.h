#ifndef RPC_FTS_RAMP_MDE_COUNT_TEST_H
#define RPC_FTS_RAMP_MDE_COUNT_TEST_H

#include <test/UnitTestRunner.h>
#include <rp/system/System.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/fts/compressor/Compressor.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldComparison.h>
#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): MDE solutions of BD and MC runs with
* a LinearRamp. The steps leave the w fields unchanged (zero mobility or
* move size), so that the compressor of a step needs one MDE solution.
*/
class FtsRampMdeCountTest : public LogFileUnitTest
{
public:

   /*
   * Run 4 BD or MC steps with the given perturbation and ramp blocks.
   * Check that c and H are those of the final w and parameters, and
   * return the number of MDE solutions after the initial compression.
   */
   template <class SimulatorT>
   int run(bool bd, std::string const & blocks)
   {
      if (!logFile_.is_open()) openLogFile("out/FtsRampMdeCount.log");
      Rp::System<3,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      std::ifstream file;
      openInputFile("in/Fts/param_system_small", file);
      system.readParam(file);
      SimulatorT simulator(system);
      std::istringstream in((bd ? "BdSimulator{ seed 5 "
                  "ExplicitBdStep{ mobility 0 } " : "McSimulator{ seed 5 "
                  "McMoveManager{ RealMove{ probability 1 sigma 0 } } ")
                  + std::string("LrAmCompressor{ epsilon 1e-8 maxItr 200 "
                  "maxHist 20 verbose 0 errorType rms } ") + blocks + " }");
      simulator.readParam(in);
      system.w().readRGrid("in/Fts/w_small.rf");
      simulator.simulate(4);
      TEST_ASSERT(simulator.iStep() == 4);
      int count = simulator.compressor().mdeCounter();
      Log::file() << "MDE solutions in 4 steps: " << count << std::endl;

      RField<3,CPT> c;
      c.allocate(system.c().rgrid(0).meshDimensions());
      c = system.c().rgrid(0);
      double hamiltonian = simulator.hamiltonian();
      system.compute();
      simulator.computeWc();
      simulator.computeHamiltonian();
      RFieldComparison<3,CPT> comparison;
      comparison.compare(c, system.c().rgrid(0));
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
      double diff = std::abs(simulator.hamiltonian() - hamiltonian);
      TEST_ASSERT(diff < 1.0E-10*std::abs(hamiltonian));
      return count;
   }

   std::string ramp(std::string const & parameter)
   {  return "LinearRamp{ nParameter 1 parameters[ " + parameter + " ] }"; }

   // A ramp of chi (10 to 12) or of the lambda of a perturbation (0 to
   // 1) changes neither c(w) nor the compression: one MDE solution per
   // step. Each ramp step solved the MDE twice more.
   void testChiRamp()
   {
      printMethod(TEST_FUNC);
      std::string chi = ramp("chi 0 1 2.0");
      TEST_ASSERT((run< Rp::BdSimulator<3,CPT> >(true, chi) == 4));
      TEST_ASSERT((run< Rp::McSimulator<3,CPT> >(false, chi) == 4));
   }

   void testLambdaRampBd()
   {
      printMethod(TEST_FUNC);
      std::string ec = "EinsteinCrystalPerturbation{ lambda 0 "
                       "fieldFileName in/Fts/w_small_ref.rf } ";
      int n = run< Rp::BdSimulator<3,CPT> >(true, ec + ramp("lambda_pert 1"));
      TEST_ASSERT(n == 4);
   }

   // A ramp of a segment length (1.0 to 1.1) changes c(w): the fields
   // are compressed again after each ramp step.
   void testKuhnRampBd()
   {
      printMethod(TEST_FUNC);
      int n = run< Rp::BdSimulator<3,CPT> >(true, ramp("kuhn 0 0.1"));
      TEST_ASSERT(n >= 8);
   }
};

TEST_BEGIN(FtsRampMdeCountTest)
TEST_ADD(FtsRampMdeCountTest, testChiRamp)
TEST_ADD(FtsRampMdeCountTest, testLambdaRampBd)
TEST_ADD(FtsRampMdeCountTest, testKuhnRampBd)
TEST_END(FtsRampMdeCountTest)

#endif
