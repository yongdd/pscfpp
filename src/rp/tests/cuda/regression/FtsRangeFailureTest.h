#ifndef RP_CUT_FTS_RANGE_FAILURE_TEST_H
#define RP_CUT_FTS_RANGE_FAILURE_TEST_H
#include <rp/tests/cpp/regression/FtsRangeFailureChecks.h>
#include <test/UnitTestRunner.h>
class FtsRangeFailureTest
   : public Pscf::Rp::FtsRangeFailureChecks<Pscf::CUT>
{
public:
   void testCompressors()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsRangeFailureCompressors.log");
      checkCompressors();
   }
   void testMcMoves()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsRangeFailureMcMoves.log");
      checkMcMoves();
   }
   void testBdStep()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsRangeFailureBdStep.log");
      checkBdStep();
   }
};
TEST_BEGIN(FtsRangeFailureTest)
TEST_ADD(FtsRangeFailureTest, testCompressors)
TEST_ADD(FtsRangeFailureTest, testMcMoves)
TEST_ADD(FtsRangeFailureTest, testBdStep)
TEST_END(FtsRangeFailureTest)
#endif
