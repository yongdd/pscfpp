#ifndef RP_CUT_SCFT_CANONICAL_REL_NORM_TEST_H
#define RP_CUT_SCFT_CANONICAL_REL_NORM_TEST_H
#include <rp/tests/cpp/regression/ScftCanonicalRelNormChecks.h>
#include <test/UnitTestRunner.h>
class ScftCanonicalRelNormTest
   : public Pscf::Rp::ScftCanonicalRelNormChecks<Pscf::CUT>
{
public:
   void testPressureField()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftCanonicalRelNormPressure.log");
      checkPressureField("AmIteratorBasis");
      checkPressureField("AmIteratorGrid");
   }
   void testExternalField()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftCanonicalRelNormExternal.log");
      checkExternalField("AmIteratorBasis");
      checkExternalField("AmIteratorGrid");
   }
   void testDivergedIteration()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftCanonicalRelNormDiverged.log");
      checkDivergedIteration("AmIteratorBasis");
      checkDivergedIteration("AmIteratorGrid");
   }
};
TEST_BEGIN(ScftCanonicalRelNormTest)
TEST_ADD(ScftCanonicalRelNormTest, testPressureField)
TEST_ADD(ScftCanonicalRelNormTest, testExternalField)
TEST_ADD(ScftCanonicalRelNormTest, testDivergedIteration)
TEST_END(ScftCanonicalRelNormTest)
#endif
