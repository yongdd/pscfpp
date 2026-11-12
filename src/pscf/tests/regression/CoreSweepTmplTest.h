#ifndef PSCF_CORE_SWEEP_TMPL_TEST_H
#define PSCF_CORE_SWEEP_TMPL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/sweep/SweepTmpl.h>
#include <util/misc/Log.h>

#include <fstream>
#include <sstream>
#include <vector>

using namespace Pscf;
using namespace Util;

/*
* Sweep with trivial states (for CoreSweepTmplTest).
*
* Records the value of s of every accepted state. The solve function
* fails at chosen attempts (counted from 1, the initial state being
* attempt 0), so that the sweep halves its step.
*/
class CoreSweepTmplMock : public SweepTmpl<double>
{

public:

   CoreSweepTmplMock()
    : SweepTmpl<double>(2),
      nAttempt_(-1)
   {}

   /*
   * Read the parameters, with ns steps.
   */
   void readNs(int ns)
   {
      std::ostringstream text;
      text << "ns  " << ns << "\n"
           << "baseFileName  out/\n"
           << "historyCapacity  2\n"
           << "reuseState  1\n";
      std::istringstream in(text.str());
      readParameters(in);
   }

   /// Attempts (1, 2, ...) at which solve fails.
   std::vector<int> failures;

   /// Values of s of the accepted states.
   std::vector<double> sAccepted;

protected:

   void checkAllocation(double& state)
   {}

   void setup()
   {
      initialize();
      nAttempt_ = -1;
      sAccepted.clear();
   }

   void setParameters(double sNew)
   {}

   void extrapolate(double sNew)
   {}

   int solve(bool isContinuation)
   {
      ++nAttempt_;
      for (size_t i = 0; i < failures.size(); ++i) {
         if (failures[i] == nAttempt_) return 1;
      }
      return 0;
   }

   void reset()
   {}

   void getSolution()
   {  sAccepted.push_back(s(0)); }

private:

   int nAttempt_;

};

/*
* Regression tests (pscfpp fixes): the path variable s of SweepTmpl.
*/
class CoreSweepTmplTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {
      if (logFile_.is_open()) {
         Log::close();
      }
   }

   /*
   * Run a sweep with ns steps, failing at the given attempts.
   */
   void runSweep(CoreSweepTmplMock& sweep, int ns,
                 std::vector<int> const & failures)
   {
      if (!logFile_.is_open()) {
         openOutputFile("out/CoreSweepTmpl.log", logFile_);
         Log::setFile(logFile_);
      }
      sweep.readNs(ns);
      sweep.failures = failures;
      sweep.sweep();
   }

   /*
   * Without failures, the states are at s = i/ns exactly, for
   * i = 0, ..., ns. Repeated addition of ds = 1/ns ended at
   * s = 1 + 2.2e-16 for ns = 9 and at 1 - 1.1e-16 for ns = 10.
   */
   void testExactS()
   {
      printMethod(TEST_FUNC);

      int nsList[4] = {9, 10, 3, 49};
      for (int j = 0; j < 4; ++j) {
         int ns = nsList[j];
         CoreSweepTmplMock sweep;
         runSweep(sweep, ns, std::vector<int>());
         TEST_ASSERT((int) sweep.sAccepted.size() == ns + 1);
         for (int i = 0; i <= ns; ++i) {
            TEST_ASSERT(sweep.sAccepted[i] == double(i)/double(ns));
         }
         TEST_ASSERT(sweep.sAccepted[ns] == 1.0);
      }
   }

   /*
   * After a failure, the step is halved, and the states are at exact
   * multiples of the new step, the last one at s = 1.
   */
   void testExactSHalved()
   {
      printMethod(TEST_FUNC);

      // ns = 3, attempt 2 (s = 2/3) fails: s = 0, 1/3, 1/2, 2/3, 5/6, 1
      {
         CoreSweepTmplMock sweep;
         std::vector<int> failures(1, 2);
         runSweep(sweep, 3, failures);
         TEST_ASSERT(sweep.sAccepted.size() == 6);
         int k[6] = {0, 2, 3, 4, 5, 6};
         for (int i = 0; i < 6; ++i) {
            TEST_ASSERT(sweep.sAccepted[i] == double(k[i])/6.0);
         }
         TEST_ASSERT(sweep.sAccepted[5] == 1.0);
      }

      // ns = 9, attempts 3 and 4 fail (step 1/9, then 1/18, then
      // 1/36): s = 0, 1/9, 2/9, 9/36, 10/36, ..., 36/36
      {
         CoreSweepTmplMock sweep;
         std::vector<int> failures;
         failures.push_back(3);
         failures.push_back(4);
         runSweep(sweep, 9, failures);
         TEST_ASSERT(sweep.sAccepted.size() == 31);
         TEST_ASSERT(sweep.sAccepted[1] == 1.0/9.0);
         TEST_ASSERT(sweep.sAccepted[2] == 2.0/9.0);
         for (int i = 3; i < 31; ++i) {
            TEST_ASSERT(sweep.sAccepted[i] == double(i + 6)/36.0);
         }
         TEST_ASSERT(sweep.sAccepted[30] == 1.0);
      }
   }

private:

   std::ofstream logFile_;

};

TEST_BEGIN(CoreSweepTmplTest)
TEST_ADD(CoreSweepTmplTest, testExactS)
TEST_ADD(CoreSweepTmplTest, testExactSHalved)
TEST_END(CoreSweepTmplTest)

#endif
