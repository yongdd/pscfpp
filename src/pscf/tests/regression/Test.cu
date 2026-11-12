/*
* This program runs the regression tests in pscf/tests/regression.
*/

#include <util/global.h>
#include "RegressionTestComposite.h"

#include <util/param/BracketPolicy.h>
#include <test/TestRunner.h>
#include <test/CompositeTestRunner.h>

using namespace Util;

int main(int argc, char* argv[])
{
   BracketPolicy::set(BracketPolicy::Optional);
   RegressionTestComposite runner;

   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
   }

   int failures = runner.run();
   if (failures != 0) {
      failures = 1;
   }
   return failures;
}
