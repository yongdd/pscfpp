#ifndef RP_CPT_FLEXIBLE_PARAMS_REUSE_TEST_H
#define RP_CPT_FLEXIBLE_PARAMS_REUSE_TEST_H
#include <rp/tests/FlexibleParamsReuseChecks.h>
#include <test/UnitTestRunner.h>

class FlexibleParamsReuseTest
   : public Pscf::Rp::FlexibleParamsReuseChecks<Pscf::CPT>
{
public:
   void testBasis() { printMethod(TEST_FUNC); check(true); }
   void testGrid() { printMethod(TEST_FUNC); check(false); }
};
TEST_BEGIN(FlexibleParamsReuseTest)
TEST_ADD(FlexibleParamsReuseTest, testBasis)
TEST_ADD(FlexibleParamsReuseTest, testGrid)
TEST_END(FlexibleParamsReuseTest)
#endif
