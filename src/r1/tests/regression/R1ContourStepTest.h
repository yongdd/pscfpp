#ifndef R1_CONTOUR_STEP_TEST_H
#define R1_CONTOUR_STEP_TEST_H

#include "R1TestBase.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

/*
* A contour step so small that the number of steps of a block does
* not fit an int.
*/
class R1ContourStepTest : public R1TestBase
{

public:

   /*
   * Blocks of length 0.5 with ds = 1E-10 (5E9 steps): reading the
   * parameter file must stop with an Exception.
   */
   void testTinyDs()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1ContourStepTestTinyDs.log");

      System sys;
      bool thrown = false;
      try {
         readSystem(sys, "in/R1/tiny_ds.prm");
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(thrown);
   }

};

TEST_BEGIN(R1ContourStepTest)
TEST_ADD(R1ContourStepTest, testTinyDs)
TEST_END(R1ContourStepTest)

#endif
