#ifndef RP_CPT_BEAD_FIELD_NORMALIZATION_TEST_H
#define RP_CPT_BEAD_FIELD_NORMALIZATION_TEST_H
#include <rp/tests/cpp/regression/BeadFieldNormalizationChecks.h>
#include <test/UnitTestRunner.h>
class BeadFieldNormalizationTest
   : public Pscf::Rp::BeadFieldNormalizationChecks<Pscf::CPT>
{
public:
   void testLinear() { printMethod(TEST_FUNC); checkLinear(); }
   void testBranched() { printMethod(TEST_FUNC); checkBranched(); }
   void testRange() { printMethod(TEST_FUNC); checkRangeAndRecovery(); }
};
TEST_BEGIN(BeadFieldNormalizationTest)
TEST_ADD(BeadFieldNormalizationTest, testLinear)
TEST_ADD(BeadFieldNormalizationTest, testBranched)
TEST_ADD(BeadFieldNormalizationTest, testRange)
TEST_END(BeadFieldNormalizationTest)
#endif
