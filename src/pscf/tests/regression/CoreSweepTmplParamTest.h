#ifndef PSCF_CORE_SWEEP_TMPL_PARAM_TEST_H
#define PSCF_CORE_SWEEP_TMPL_PARAM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/sweep/SweepTmpl.h>

#include <sstream>
#include <string>

using namespace Pscf;
using namespace Util;

/*
* Sweep with trivial states (for CoreSweepTmplParamTest), which only
* reads its parameters.
*/
class CoreSweepTmplParamMock : public SweepTmpl<double>
{

public:

   CoreSweepTmplParamMock()
    : SweepTmpl<double>(2)
   {}

   /*
   * Read the parameters from a string.
   */
   void readText(std::string const & text)
   {
      std::istringstream in(text);
      readParameters(in);
   }

protected:

   void checkAllocation(double& state)
   {}

   void setup()
   {  initialize(); }

   void setParameters(double sNew)
   {}

   void extrapolate(double sNew)
   {}

   int solve(bool isContinuation)
   {  return 0; }

   void reset()
   {}

   void getSolution()
   {}

};

/*
* Regression tests (pscfpp fixes): parameters read by SweepTmpl.
*/
class CoreSweepTmplParamTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   /*
   * Does reading these parameters throw an Exception?
   */
   bool readThrows(std::string const & text)
   {
      CoreSweepTmplParamMock sweep;
      try {
         sweep.readText(text);
      } catch (Exception& e) {
         return true;
      }
      return false;
   }

   /*
   * The number of steps ns must be positive. With ns = 0 the step was
   * 1/0, and a negative ns ended the sweep after a step backwards.
   */
   void testNs()
   {
      printMethod(TEST_FUNC);
      TEST_ASSERT(readThrows("ns  0\n"));
      TEST_ASSERT(readThrows("ns  -2\n"));
      TEST_ASSERT(!readThrows("ns  1\n"));
      TEST_ASSERT(!readThrows("ns  10\nbaseFileName out/\n"));
   }

   /*
   * The history capacity must be positive.
   */
   void testHistoryCapacity()
   {
      printMethod(TEST_FUNC);
      TEST_ASSERT(readThrows("ns  2\nhistoryCapacity  0\n"));
      TEST_ASSERT(readThrows("ns  2\nhistoryCapacity  -1\n"));
      TEST_ASSERT(!readThrows("ns  2\nhistoryCapacity  1\n"));
   }

};

TEST_BEGIN(CoreSweepTmplParamTest)
TEST_ADD(CoreSweepTmplParamTest, testNs)
TEST_ADD(CoreSweepTmplParamTest, testHistoryCapacity)
TEST_END(CoreSweepTmplParamTest)

#endif
