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
#include "R1SweepTest.h"
#include "R1HomogeneousTest.h"
#include "R1HomogeneousLimitTest.h"
#include "R1LinearSweepTest.h"
#include "R1HomogeneousUnderflowTest.h"

TEST_COMPOSITE_BEGIN(RegressionTestComposite)
TEST_COMPOSITE_ADD_UNIT(CoreMixtureTest);
TEST_COMPOSITE_ADD_UNIT(R1IteratorTest);
TEST_COMPOSITE_ADD_UNIT(R1BlockTest);
TEST_COMPOSITE_ADD_UNIT(R1CommandTest);
TEST_COMPOSITE_ADD_UNIT(R1SweepTest);
TEST_COMPOSITE_ADD_UNIT(R1HomogeneousTest);
TEST_COMPOSITE_ADD_UNIT(R1HomogeneousLimitTest);
TEST_COMPOSITE_ADD_UNIT(R1LinearSweepTest);
TEST_COMPOSITE_ADD_UNIT(R1HomogeneousUnderflowTest);
TEST_COMPOSITE_END

#endif
