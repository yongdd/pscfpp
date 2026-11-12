#ifndef RPG_REGRESSION_TEST_COMPOSITE_H
#define RPG_REGRESSION_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "ScftAmIteratorTest.h"
#include "FieldIoRegressionTest.h"

TEST_COMPOSITE_BEGIN(RegressionTestComposite)
TEST_COMPOSITE_ADD_UNIT(ScftAmIteratorTest);
TEST_COMPOSITE_ADD_UNIT(FieldIoRegressionTest);
TEST_COMPOSITE_END

#endif
