#ifndef RP_CPT_SCFT_FIELD_SHIFT_TEST_H
#define RP_CPT_SCFT_FIELD_SHIFT_TEST_H
#include <rp/tests/cpp/regression/ScftFieldShiftChecks.h>
#include <test/UnitTestRunner.h>
class ScftFieldShiftTest
   : public Pscf::Rp::ScftFieldShiftChecks<Pscf::CPT>
{
public:
   void testThreadStrongFields()
   {  printMethod(TEST_FUNC); checkThreadStrongFields(); }
   void testBeadSlices()
   {  printMethod(TEST_FUNC); checkBeadSlices(); }
   void testRaisedShift()
   {  printMethod(TEST_FUNC); checkRaisedShift(); }
};
TEST_BEGIN(ScftFieldShiftTest)
TEST_ADD(ScftFieldShiftTest, testThreadStrongFields)
TEST_ADD(ScftFieldShiftTest, testBeadSlices)
TEST_ADD(ScftFieldShiftTest, testRaisedShift)
TEST_END(ScftFieldShiftTest)
#endif
