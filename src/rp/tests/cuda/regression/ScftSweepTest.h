#ifndef RPG_SCFT_SWEEP_TEST_H
#define RPG_SCFT_SWEEP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/scft/sweep/Sweep.h>
#include <rp/scft/sweep/LinearSweep.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Solvent.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>

#include <util/misc/Exception.h>
#include <util/param/Label.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests for SCFT parameter sweeps.
*/
class ScftSweepTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Sweep of phi_polymer in a mixed ensemble (closed polymer, open
   * solvent). The open solvent takes up the change in volume fraction,
   * so the swept phi values need not add up to zero.
   */
   void testLinearSweepPhiMixed()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftSweepPhiMixed.log");

      Rp::System<1,CUT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      bool thrown = false;
      try {
         std::ifstream in;
         openInputFile("in/Scft/mixed_param", in);
         system.readParam(in);
         in.close();
      } catch (Exception& e) {
         thrown = true;
         Label::clear(); // Discard any label left by the failed read
      }
      TEST_ASSERT(!thrown);
      TEST_ASSERT(!system.mixture().isCanonical());

      // Sweep from phi_polymer = 0.50 to 0.55 in 2 steps
      system.w().readBasis("../sweep/in/solvent/w.bf");
      system.sweep();

      // Check that the sweep reached the end without backtracking
      std::ifstream f;
      f.open(std::string(filePrefix()
                         + "out/ScftSweepPhiMixed_2_w.bf").c_str());
      TEST_ASSERT(f.good());
      f.close();
      f.open(std::string(filePrefix()
                         + "out/ScftSweepPhiMixed_3_w.bf").c_str());
      TEST_ASSERT(!f.good());

      // Check the final volume fractions
      double phiP = system.mixture().polymer(0).phi();
      double phiS = system.mixture().solvent(0).phi();
      if (verbose() > 0) {
         std::cout << "\nphi polymer = " << phiP
                   << "\nphi solvent = " << phiS << "\n";
      }
      TEST_ASSERT(std::abs(phiP - 0.55) < 1.0E-10);
      TEST_ASSERT(std::abs(phiP + phiS - 1.0) < 1.0E-6);
   }

};

TEST_BEGIN(ScftSweepTest)
TEST_ADD(ScftSweepTest, testLinearSweepPhiMixed)
TEST_END(ScftSweepTest)

#endif
