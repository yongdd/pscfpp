#ifndef PSCF_REGRESSION_TEST_COMPOSITE_H
#define PSCF_REGRESSION_TEST_COMPOSITE_H

/*
* Regression tests for fixed bugs, kept apart from the unit tests.
*/

#include <test/CompositeTestRunner.h>

#include "CoreBasisTest.h"
#include "CoreLuSolverTest.h"
#include "CoreFhInteractionTest.h"
#include "CoreAmbdInteractionTest.h"

TEST_COMPOSITE_BEGIN(RegressionTestComposite)
TEST_COMPOSITE_ADD_UNIT(CoreBasisTest);
TEST_COMPOSITE_ADD_UNIT(CoreLuSolverTest);
TEST_COMPOSITE_ADD_UNIT(CoreFhInteractionTest);
TEST_COMPOSITE_ADD_UNIT(CoreAmbdInteractionTest);
TEST_COMPOSITE_END

#endif
