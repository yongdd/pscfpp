#ifndef PSCF_REGRESSION_TEST_COMPOSITE_H
#define PSCF_REGRESSION_TEST_COMPOSITE_H

/*
* Regression tests for fixed bugs, kept apart from the unit tests.
*/

#include <test/CompositeTestRunner.h>

#include "CoreBasisTest.h"

TEST_COMPOSITE_BEGIN(RegressionTestComposite)
TEST_COMPOSITE_ADD_UNIT(CoreBasisTest);
TEST_COMPOSITE_END

#endif
