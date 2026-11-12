#ifndef R1_NR_SHELL_TEST_H
#define R1_NR_SHELL_TEST_H

#include "R1TestBase.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/field/Domain.h>
#include <r1/iterator/Iterator.h>
#include <util/math/Constants.h>

#include <cmath>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

class R1NrShellTest : public R1TestBase
{

public:

   /*
   * Set initial w fields for a micelle with an A core (chi = 80).
   */
   void setGuess(System& sys)
   {
      int nx = sys.domain().nx();
      double chi = 80.0;
      double cs;
      for (int i = 0; i < nx; ++i) {
         cs = std::cos(Constants::Pi*double(i)/double(nx-1));
         sys.wField(0)[i] = -0.5*chi*cs;
         sys.wField(1)[i] = +0.5*chi*cs;
      }
      double shift = sys.wField(1)[nx-1];
      for (int i = 0; i < nx; ++i) {
         sys.wField(0)[i] -= shift;
         sys.wField(1)[i] -= shift;
      }
   }

   /*
   * NrIterator for a spherical shell with xMin = 1e-9 and dx = 0.006.
   * The inner grid point has the weight xMin^2/2 and is tied to the
   * next point by the factor (1 + dx/(2 xMin))^2 of the Crank-Nicolson
   * matrix, so a change of w there changes no concentration, and the
   * jacobian has two dependent columns. The LU decomposition of GSL
   * 2.7 and later then read pivot indices it had never set, and
   * pscf_1d crashed (segmentation fault). The iteration must converge
   * to the solution for xMin = 1e-6, whose jacobian is regular.
   */
   void testTinyInnerRadius()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1NrShellTestTinyInnerRadius.log");

      System ref;
      readSystem(ref, "in/R1/shell_nr_6.prm");
      setGuess(ref);
      TEST_ASSERT(ref.iterate() == 0);

      System sys;
      readSystem(sys, "in/R1/shell_nr_9.prm");
      setGuess(sys);
      TEST_ASSERT(sys.iterate() == 0);

      Log::file() << "fHelmholtz (xMin = 1e-9, 1e-6) = " 
                  << sys.fHelmholtz() << "  " << ref.fHelmholtz() 
                  << std::endl;
      TEST_ASSERT(std::abs(sys.fHelmholtz() - ref.fHelmholtz()) < 1.0E-7);
      int nx = sys.domain().nx();
      double dc, dcMax = 0.0;
      for (int im = 0; im < 2; ++im) {
         for (int i = 0; i < nx; ++i) {
            dc = std::abs(sys.cField(im)[i] - ref.cField(im)[i]);
            if (!(dc <= dcMax)) dcMax = dc;
         }
      }
      Log::file() << "max |c(xMin = 1e-9) - c(xMin = 1e-6)| = " << dcMax
                  << std::endl;
      TEST_ASSERT(dcMax < 1.0E-4);
   }

};

TEST_BEGIN(R1NrShellTest)
TEST_ADD(R1NrShellTest, testTinyInnerRadius)
TEST_END(R1NrShellTest)

#endif
