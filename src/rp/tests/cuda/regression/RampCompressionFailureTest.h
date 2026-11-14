#ifndef RP_CUT_RAMP_COMPRESSION_FAILURE_TEST_H
#define RP_CUT_RAMP_COMPRESSION_FAILURE_TEST_H

#include <rp/tests/cpp/regression/RampCompressionChecks.h>
#include <test/UnitTestRunner.h>

class RampCompressionFailureTest
   : public Pscf::Rp::RampCompressionChecks<Pscf::CUT>
{
public:
   void testBd()
   { printMethod(TEST_FUNC); checkSimulator(true); }

   void testMc()
   { printMethod(TEST_FUNC); checkSimulator(false); }
};

TEST_BEGIN(RampCompressionFailureTest)
TEST_ADD(RampCompressionFailureTest, testBd)
TEST_ADD(RampCompressionFailureTest, testMc)
TEST_END(RampCompressionFailureTest)

#endif
