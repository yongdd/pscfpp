#ifndef PSCF_REGRESSION_TEST_COMPOSITE_H
#define PSCF_REGRESSION_TEST_COMPOSITE_H

/*
* Regression tests for fixed bugs, kept apart from the unit tests.
* Classes with the prefix Gpu test CUDA code (only with PSCF_CUDA).
*/

#include <test/CompositeTestRunner.h>

#include "CoreBasisTest.h"
#include "CoreLuSolverTest.h"
#include "CoreFhInteractionTest.h"
#include "CoreAmbdInteractionTest.h"
#include "CorePolymerTopologyTest.h"
#include "CoreEdgeTest.h"
#include "CoreFhMixtureTest.h"

#ifdef PSCF_CUDA
#include "GpuWaveListTest.h"
#include "GpuVecRandomTest.h"
#include "GpuFftTest.h"
#endif

TEST_COMPOSITE_BEGIN(RegressionTestComposite)
TEST_COMPOSITE_ADD_UNIT(CoreBasisTest);
TEST_COMPOSITE_ADD_UNIT(CoreLuSolverTest);
TEST_COMPOSITE_ADD_UNIT(CoreFhInteractionTest);
TEST_COMPOSITE_ADD_UNIT(CoreAmbdInteractionTest);
TEST_COMPOSITE_ADD_UNIT(CorePolymerTopologyTest);
TEST_COMPOSITE_ADD_UNIT(CoreEdgeTest);
TEST_COMPOSITE_ADD_UNIT(CoreFhMixtureTest);
#ifdef PSCF_CUDA
TEST_COMPOSITE_ADD_UNIT(GpuWaveListTest);
TEST_COMPOSITE_ADD_UNIT(GpuVecRandomTest);
TEST_COMPOSITE_ADD_UNIT(GpuFftTest);
#endif
TEST_COMPOSITE_END

#endif
