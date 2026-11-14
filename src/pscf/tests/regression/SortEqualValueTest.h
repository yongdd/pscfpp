#ifndef PSCF_SORT_EQUAL_VALUE_TEST_H
#define PSCF_SORT_EQUAL_VALUE_TEST_H
#include <pscf/math/Sort.h>
#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <vector>

class SortEqualValueTest : public UnitTest
{
public:
   void testFloatGroups()
   {
      printMethod(TEST_FUNC);
      check<float>(1.0e-8f);
      check<float>(1.0e-4f);
   }
   void testDoubleGroups()
   {
      printMethod(TEST_FUNC);
      check<double>(1.0e-17);
      check<double>(1.0e-8);
   }
private:
   template <typename T>
   void check(T epsilon)
   {
      std::vector<Pscf::Sort::Item<T> > values =
         {{T(2),0}, {T(1),1}, {T(1),2}, {T(3),3}, {T(2),4}};
      Pscf::Sort::sort(values);
      Util::GArray<Pscf::Sort::Bunch> bunches;
      Pscf::Sort::findBunches(values,bunches,epsilon);
      TEST_ASSERT(bunches.size() == 3);
      TEST_ASSERT(bunches[0][0] == 0 && bunches[0][1] == 2);
      TEST_ASSERT(bunches[1][0] == 2 && bunches[1][1] == 4);
      TEST_ASSERT(bunches[2][0] == 4 && bunches[2][1] == 5);
      // Retain grouping relative to the first value, not chained neighbors.
      values = {{T(0),0}, {epsilon*T(0.75),1}, {epsilon*T(1.5),2}};
      Pscf::Sort::findBunches(values,bunches,epsilon);
      TEST_ASSERT(bunches.size() == 2);
      TEST_ASSERT(bunches[0][0] == 0 && bunches[0][1] == 2);
      TEST_ASSERT(bunches[1][0] == 2 && bunches[1][1] == 3);
   }
};
TEST_BEGIN(SortEqualValueTest)
TEST_ADD(SortEqualValueTest, testFloatGroups)
TEST_ADD(SortEqualValueTest, testDoubleGroups)
TEST_END(SortEqualValueTest)
#endif
