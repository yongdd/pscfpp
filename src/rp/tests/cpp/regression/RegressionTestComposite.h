#ifndef RPC_REGRESSION_TEST_COMPOSITE_H
#define RPC_REGRESSION_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "ScftAmIteratorTest.h"
#include "FieldIoRegressionTest.h"
#include "FilmStressTest.h"
#include "FilmMaskCheckTest.h"
#include "FilmSweepTest.h"
#include "ScftSweepTest.h"
#include "CommandTest.h"

TEST_COMPOSITE_BEGIN(RegressionTestComposite)
TEST_COMPOSITE_ADD_UNIT(ScftAmIteratorTest);
TEST_COMPOSITE_ADD_UNIT(FieldIoRegressionTest);
TEST_COMPOSITE_ADD_UNIT(FilmStressTest);
TEST_COMPOSITE_ADD_UNIT(FilmMaskCheckTest);
TEST_COMPOSITE_ADD_UNIT(FilmSweepTest);
TEST_COMPOSITE_ADD_UNIT(ScftSweepTest);
TEST_COMPOSITE_ADD_UNIT(CommandTest);
TEST_COMPOSITE_END

#endif
