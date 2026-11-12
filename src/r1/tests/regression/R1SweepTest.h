#ifndef R1_SWEEP_TEST_H
#define R1_SWEEP_TEST_H

#include "R1TestBase.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <r1/solvers/Solvent.h>

#include <cmath>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

class R1SweepTest : public R1TestBase
{

public:

   /*
   * LinearSweep of a polymer volume fraction in a mixed ensemble.
   */
   void testPhiMixedEnsemble()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1SweepTestPhiMixedEnsemble.log");

      // Closed diblock and open solvent: the sweep changes the diblock
      // volume fraction from 0.95 to 0.92, the solvent takes up the rest
      System sys;
      bool thrown = false;
      try {
         readSystem(sys, "in/R1/planar_mixed_sweep.prm");
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
      runCommands(sys, "READ_W in/R1/planar.w\nSWEEP\nFINISH\n");
      TEST_ASSERT(eq(sys.mixture().polymer(0).phi(), 0.92));

      // Direct solution for a diblock volume fraction of 0.92
      System ref;
      readSystem(ref, "in/R1/planar_mixed.prm");
      runCommands(ref, "READ_W in/R1/planar.w\nITERATE\nFINISH\n");

      double phiS = sys.mixture().solvent(0).phi();
      double phiSRef = ref.mixture().solvent(0).phi();
      Log::file() << "Sweep:  fHelmholtz = " << sys.fHelmholtz()
                  << ", solvent phi = " << phiS << "\n"
                  << "Direct: fHelmholtz = " << ref.fHelmholtz()
                  << ", solvent phi = " << phiSRef << std::endl;
      TEST_ASSERT(std::abs(sys.fHelmholtz() - ref.fHelmholtz()) < 1.0E-6);
      TEST_ASSERT(std::abs(sys.pressure() - ref.pressure()) < 1.0E-6);
      TEST_ASSERT(std::abs(phiS - phiSRef) < 1.0E-6);
   }

};

TEST_BEGIN(R1SweepTest)
TEST_ADD(R1SweepTest, testPhiMixedEnsemble)
TEST_END(R1SweepTest)

#endif
