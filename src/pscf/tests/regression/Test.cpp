/*
* This program runs the regression tests in the pscf/tests/regression
* directory.
*/ 

#include <util/global.h>
#include <util/param/BracketPolicy.h>
#include "RegressionTestComposite.h"

#include <test/TestRunner.h>
#include <test/CompositeTestRunner.h>

using namespace Pscf;
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
   return (failures != 0);
}
