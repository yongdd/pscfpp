#ifndef RP_CPT_SCFT_RANGE_FAILURE_TEST_H
#define RP_CPT_SCFT_RANGE_FAILURE_TEST_H
#include <rp/tests/cpp/regression/ScftRangeFailureChecks.h>
#include <test/UnitTestRunner.h>
class ScftRangeFailureTest
   : public Pscf::Rp::ScftRangeFailureChecks<Pscf::CPT>
{
public:
   void testIterateBasis()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftRangeFailureBasis.log");
      checkIterate("AmIteratorBasis");
   }
   void testIterateGrid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftRangeFailureGrid.log");
      checkIterate("AmIteratorGrid");
   }
   void testCommand()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftRangeFailureCommand.log");
      checkCommand();
   }
   void testSweep()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftRangeFailureSweep.log");
      checkSweep();
   }
};
TEST_BEGIN(ScftRangeFailureTest)
TEST_ADD(ScftRangeFailureTest, testIterateBasis)
TEST_ADD(ScftRangeFailureTest, testIterateGrid)
TEST_ADD(ScftRangeFailureTest, testCommand)
TEST_ADD(ScftRangeFailureTest, testSweep)
TEST_END(ScftRangeFailureTest)
#endif
