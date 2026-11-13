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
#include "R1NrShellTest.h"
#include "R1GaugeTest.h"
#include "R1ContourStepTest.h"
#include "R1AmRelNormTest.h"
#include "R1CommandStatusTest.h"
#include "R1SweepCellParamTest.h"

#include "FieldReadTransactionTest.h"

#include "DerivedOutputValidityTest.h"

#include "MeshTransformBoundsTest.h"

#include "DomainResizeTest.h"

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
TEST_COMPOSITE_ADD_UNIT(R1NrShellTest);
TEST_COMPOSITE_ADD_UNIT(R1GaugeTest);
TEST_COMPOSITE_ADD_UNIT(R1ContourStepTest);
TEST_COMPOSITE_ADD_UNIT(R1AmRelNormTest);
TEST_COMPOSITE_ADD_UNIT(R1CommandStatusTest);
TEST_COMPOSITE_ADD_UNIT(R1SweepCellParamTest);
TEST_COMPOSITE_ADD_UNIT(FieldReadTransactionTest);
TEST_COMPOSITE_ADD_UNIT(DerivedOutputValidityTest);
TEST_COMPOSITE_ADD_UNIT(MeshTransformBoundsTest);
TEST_COMPOSITE_ADD_UNIT(DomainResizeTest);
TEST_COMPOSITE_END

#endif
