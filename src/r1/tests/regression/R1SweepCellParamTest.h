#ifndef R1_SWEEP_CELL_PARAM_TEST_H
#define R1_SWEEP_CELL_PARAM_TEST_H

#include "R1TestBase.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

/*
* Regression tests (pscfpp fixes): sweep parameter types of pscf_1d.
*/
class R1SweepCellParamTest : public R1TestBase
{

public:

   /*
   * A sweep of type cell_param (a parameter of a periodic unit cell)
   * must be rejected when the parameter file is read. The diblock of
   * this system is open (mu given), and the parameters are not echoed:
   * the file was then read without an error (else the "Invalid type_"
   * Exception of type() was thrown), and the sweep failed with "Invalid
   * type_" in get_.
   */
   void testCellParam()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1SweepCellParamTestCellParam.log");
      ParamComponent::setEcho(false);

      System sys;
      bool thrown = false;
      try {
         readSystem(sys, "in/R1/sweep_cell_param.prm");
      } catch (Exception& e) {
         thrown = e.message().find("cell_param is not available")
                  != std::string::npos;
      }
      TEST_ASSERT(thrown);
   }

};

TEST_BEGIN(R1SweepCellParamTest)
TEST_ADD(R1SweepCellParamTest, testCellParam)
TEST_END(R1SweepCellParamTest)

#endif
