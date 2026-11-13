#ifndef RP_CUT_SOLVENT_NORMALIZATION_TEST_H
#define RP_CUT_SOLVENT_NORMALIZATION_TEST_H

#include <rp/tests/SolventNormalizationChecks.h>
#include <test/UnitTestRunner.h>

class SolventNormalizationTest
   : public Pscf::Rp::SolventNormalizationChecks<Pscf::CUT>
{
public:
   void testClosed()
   { printMethod(TEST_FUNC); checkClosed(); }

   void testOpen()
   { printMethod(TEST_FUNC); checkOpen(); }
};

TEST_BEGIN(SolventNormalizationTest)
TEST_ADD(SolventNormalizationTest, testClosed)
TEST_ADD(SolventNormalizationTest, testOpen)
TEST_END(SolventNormalizationTest)

#endif
