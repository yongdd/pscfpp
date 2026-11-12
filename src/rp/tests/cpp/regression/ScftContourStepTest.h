#ifndef RPC_SCFT_CONTOUR_STEP_TEST_H
#define RPC_SCFT_CONTOUR_STEP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>

#include <util/misc/Exception.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;

/*
* A contour step so small that the number of contour steps of a block
* does not fit an int.
*/
class ScftContourStepTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Contour step ds = 1E-10 for a block of length 1800: reading the
   * parameter file must stop with an Exception.
   */
   void testTinyDs()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftContourStepTestTinyDs.log");
      Rp::System<1,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      bool thrown = false;
      try {
         std::ifstream in;
         openInputFile("in/Scft/tiny_ds", in);
         system.readParam(in);
         in.close();
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(thrown);
   }

};

TEST_BEGIN(ScftContourStepTest)
TEST_ADD(ScftContourStepTest, testTinyDs)
TEST_END(ScftContourStepTest)

#endif
