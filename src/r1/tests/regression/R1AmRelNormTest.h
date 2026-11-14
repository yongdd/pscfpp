#ifndef R1_AM_REL_NORM_TEST_H
#define R1_AM_REL_NORM_TEST_H

#include "R1TestBase.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/field/Domain.h>
#include <r1/iterator/Iterator.h>

#include <cmath>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

/*
* Regression tests (pscfpp fixes): the relNorm error of the r1
* AmIterator for w fields with a diverged pressure field.
*/
class R1AmRelNormTest : public R1TestBase
{

public:

   /*
   * Converged lamellar solution (diblock, chi N = 30), to which a smooth
   * field xi(x) = A cos(pi x/L) is added in both w fields. The residual
   * contains this pressure-like part of w only through the volume
   * fractions, which the Crank-Nicolson solver keeps bounded for large
   * fields. A diverged AM iteration of this system (maxHist 50) blew it
   * up to |xi| ~ 1e8 and stopped with "Converged" and fHelmholtz 3.9e4
   * (correct 3.79), because the relNorm error |R|/|w| is about 2e-9 for
   * A = 1e10. With maxItr = 1, solve() returns 0 only if the error of
   * the given state is below epsilon (1e-7).
   */
   void testPressureField()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1AmRelNormTestPressureField.log");

      System sys;
      readSystem(sys, "in/R1/am_lam.prm");

      // The solution itself is converged (ITERATE computes the free
      // energy only after convergence)
      runCommands(sys, "READ_W in/R1/am_lam.w\nITERATE\nFINISH\n");
      TEST_ASSERT(std::abs(sys.fHelmholtz() - 3.7910820) < 1.0E-6);

      const int nx = sys.domain().nx();
      const double pi = std::acos(-1.0);
      const int nAmp = 4;
      double amplitudes[nAmp] = {1.0E2, 1.0E5, 1.0E8, 1.0E10};
      for (int i = 0; i < nAmp; ++i) {
         double amp = amplitudes[i];
         runCommands(sys, "READ_W in/R1/am_lam.w\nFINISH\n");
         for (int k = 0; k < nx; ++k) {
            double xi = amp*std::cos(pi*double(k)/double(nx - 1));
            sys.wField(0)[k] += xi;
            sys.wField(1)[k] += xi;
         }

         // These fields are far from a solution: the total volume
         // fraction differs from 1 by more than 1e-2 somewhere
         runCommands(sys, "COMPUTE\nFINISH\n");
         double maxIncompress = 0.0;
         for (int k = 0; k < nx; ++k) {
            double r = sys.cField(0)[k] + sys.cField(1)[k] - 1.0;
            if (std::abs(r) > maxIncompress) maxIncompress = std::abs(r);
         }
         Log::file() << "A = " << amp << ": max |sum c - 1| = "
                     << maxIncompress << std::endl;
         TEST_ASSERT(maxIncompress > 1.0E-2);

         // The iterator must not report convergence
         TEST_ASSERT(sys.iterator().solve() != 0);
      }
   }

};

TEST_BEGIN(R1AmRelNormTest)
TEST_ADD(R1AmRelNormTest, testPressureField)
TEST_END(R1AmRelNormTest)

#endif
