/*
* This program runs all unit tests in the r1/tests directory.
*/ 

#include <util/global.h>
#include <util/param/BracketPolicy.h>

#include "R1dTestComposite.h"
#include "regression/RegressionTestComposite.h"

#include <test/CompositeTestRunner.h>

using namespace R1;
using namespace Util;

TEST_COMPOSITE_BEGIN(R1NsTestComposite)
addChild(new R1dTestComposite, "");
addChild(new RegressionTestComposite, "regression/");
TEST_COMPOSITE_END

int main(int argc, char* argv[])
{

   BracketPolicy::set(BracketPolicy::Optional);

   try {

      R1NsTestComposite runner;
  
      // Optionally add file prefix given as command line argument 
      if (argc > 2) {
         UTIL_THROW("Too many arguments");
      }
      if (argc == 2) {
         runner.addFilePrefix(argv[1]);
      }
   
      // Run all unit test functions
      int failures = runner.run();

      return (failures != 0);

   } catch (...) {

      std::cerr << "Uncaught exception in r1/tests/Test.cc" << std::endl;
      return 1;

   }
}
