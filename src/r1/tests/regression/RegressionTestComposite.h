#ifndef R1_REGRESSION_TEST_COMPOSITE_H
#define R1_REGRESSION_TEST_COMPOSITE_H

/*
* Regression tests for fixed bugs, kept apart from the unit tests.
*/

#include <test/CompositeTestRunner.h>

#include "CoreMixtureTest.h"
#include "R1IteratorTest.h"
#include "R1BlockTest.h"
#include "R1CommandTest.h"

TEST_COMPOSITE_BEGIN(RegressionTestComposite)
TEST_COMPOSITE_ADD_UNIT(CoreMixtureTest);
TEST_COMPOSITE_ADD_UNIT(R1IteratorTest);
TEST_COMPOSITE_ADD_UNIT(R1BlockTest);
TEST_COMPOSITE_ADD_UNIT(R1CommandTest);
TEST_COMPOSITE_END

#endif
