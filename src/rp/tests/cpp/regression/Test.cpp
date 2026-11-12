/*
* This program runs all unit tests in the rp/tests/cpp/regression directory.
*/

#include <util/global.h>
#include "RegressionTestComposite.h"

#include <util/param/BracketPolicy.h>
#include <test/TestRunner.h>
#include <test/CompositeTestRunner.h>

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
   return (failures != 0) ? 1 : 0;
}
