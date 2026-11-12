#ifndef RPC_REGRESSION_TEST_COMPOSITE_H
#define RPC_REGRESSION_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "ScftAmIteratorTest.h"

TEST_COMPOSITE_BEGIN(RegressionTestComposite)
TEST_COMPOSITE_ADD_UNIT(ScftAmIteratorTest);
TEST_COMPOSITE_END

#endif
