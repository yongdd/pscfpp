#ifndef RPG_SCFT_AM_ITERATOR_TEST_H
#define RPG_SCFT_AM_ITERATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/scft/ScftThermo.h>
#include <rp/scft/iterator/Iterator.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <rp/field/Mask.h>
#include <rp/field/FieldIo.h>
#include <pscf/interaction/Interaction.h>

#include <util/misc/Exception.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;

/*
* Regression tests for the Anderson-mixing SCFT iterators.
*/
class ScftAmIteratorTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Read a parameter file to initialize a System.
   */
   template <int D>
   void setupSystem(Rp::System<D,CUT>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

   /*
   * Solves that follow a solve that failed with a NaN.
   *
   * A large jump in chi drives the AM iteration to a NaN, which ends
   * the solve with failure. As in a Sweep that backtracks, the state
   * before the jump is then restored, and a smaller step is solved
   * with isContinuation = true. This solve must be able to start (the
   * AM timers must have been stopped after the NaN), and must not
   * reuse the AM history of the failed solve (which would apply the
   * displacement that led to the NaN at its first step).
   */
   void testBasisNanRecovery()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftAmIteratorBasisNanRecovery.log");

      // Converged solution at chi = 12
      Rp::System<1,CUT> system;
      setupSystem(system, "in/Scft/lam_nan");
      system.w().readBasis("in/Scft/lam_w.bf");
      int error = system.iterate();
      TEST_ASSERT(!error);

      // Jump to chi = 412: the solve ends with a NaN
      system.interaction().setChi(0, 1, 412.0);
      bool thrown = false;
      error = 0;
      try {
         error = system.iterate(true);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
      TEST_ASSERT(error);

      // Restore the solution at chi = 12, and solve for chi = 112
      system.w().readBasis("in/Scft/lam_w.bf");
      system.interaction().setChi(0, 1, 112.0);
      thrown = false;
      error = 1;
      try {
         error = system.iterate(true);
      } catch (Exception& e) {
         thrown = true;
      }
      if (verbose() > 0) {
         std::cout << "\nRetry: exception = " << thrown
                   << ", error = " << error;
      }
      TEST_ASSERT(!thrown);
      TEST_ASSERT(!error);
   }

};

TEST_BEGIN(ScftAmIteratorTest)
TEST_ADD(ScftAmIteratorTest, testBasisNanRecovery)
TEST_END(ScftAmIteratorTest)

#endif
