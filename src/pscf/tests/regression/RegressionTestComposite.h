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
#include "FhMixtureDependentSpeciesTest.h"
#include "CoreGroupInputTest.h"
#include "CoreKSqOverflowTest.h"
#include "CoreBasisOutputTest.h"
#include "CoreSymmetryCheckTest.h"
#include "CoreGroupFileTest.h"
#include "CoreSweepTmplTest.h"
#include "CoreAmIteratorTest.h"
#include "CoreFhMixtureUnderflowTest.h"
#include "CoreLuSolverSingularTest.h"
#include "CoreNyquistTest.h"
#include "CoreFhMixtureSingularTest.h"
#include "CoreSweepTmplParamTest.h"
#include "CoreMixAndMatchEnvTest.h"
#include "CoreSingularCellTest.h"
#include "CoreBasisDataTest.h"
#include "CoreNyquistGroupTest.h"
#include "CoreNyquistPlaneCacheTest.h"

#ifdef PSCF_CUDA
#include "GpuWaveListTest.h"
#include "GpuVecRandomTest.h"
#include "GpuFftTest.h"
#include "GpuKSqOverflowTest.h"
#include "GpuNyquistTest.h"
#include "GpuNyquistGroupTest.h"
#endif

#include "MinimumImageScaleTest.h"

#ifdef PSCF_CUDA
#include "MinimumImageScaleGpuTest.h"
#endif

#include "BeadCountValidationTest.h"

#include "FftSliceAlignmentTest.h"

#include "FieldReadEndOfFileTest.h"

#include "FourierColumnSpacingTest.h"

#include "BasisScaleTest.h"

#include "WaveShellScaleTest.h"

#ifdef PSCF_CUDA
#include "WaveShellScaleGpuTest.h"
#endif

#include "IdentityGroupReloadTest.h"

#include "AnisotropicMirrorMeshTest.h"

#ifdef PSCF_CUDA
#include "GpuComplexExpAliasTest.h"
#endif

#ifdef PSCF_CUDA
#include "GpuThreadMeshDimensionTest.h"
#endif

#include "EdgeArchiveStateTest.h"
#include "MonomerArchiveTest.h"

#include "UnitCellArchiveStateTest.h"

#include "FieldArchiveRoundTripTest.h"

#ifdef PSCF_CPP
#include "FftwComplexArrayCopyTest.h"
#endif

#include "EdgeOutputModeTest.h"

TEST_COMPOSITE_BEGIN(RegressionTestComposite)
TEST_COMPOSITE_ADD_UNIT(CoreBasisTest);
TEST_COMPOSITE_ADD_UNIT(CoreLuSolverTest);
TEST_COMPOSITE_ADD_UNIT(CoreFhInteractionTest);
TEST_COMPOSITE_ADD_UNIT(CoreAmbdInteractionTest);
TEST_COMPOSITE_ADD_UNIT(CorePolymerTopologyTest);
TEST_COMPOSITE_ADD_UNIT(CoreEdgeTest);
TEST_COMPOSITE_ADD_UNIT(CoreFhMixtureTest);
TEST_COMPOSITE_ADD_UNIT(FhMixtureDependentSpeciesTest);
TEST_COMPOSITE_ADD_UNIT(CoreGroupInputTest);
TEST_COMPOSITE_ADD_UNIT(CoreKSqOverflowTest);
TEST_COMPOSITE_ADD_UNIT(CoreBasisOutputTest);
TEST_COMPOSITE_ADD_UNIT(CoreSymmetryCheckTest);
TEST_COMPOSITE_ADD_UNIT(CoreGroupFileTest);
TEST_COMPOSITE_ADD_UNIT(CoreSweepTmplTest);
TEST_COMPOSITE_ADD_UNIT(CoreAmIteratorTest);
TEST_COMPOSITE_ADD_UNIT(CoreFhMixtureUnderflowTest);
TEST_COMPOSITE_ADD_UNIT(CoreLuSolverSingularTest);
TEST_COMPOSITE_ADD_UNIT(CoreNyquistTest);
TEST_COMPOSITE_ADD_UNIT(CoreFhMixtureSingularTest);
TEST_COMPOSITE_ADD_UNIT(CoreSweepTmplParamTest);
TEST_COMPOSITE_ADD_UNIT(CoreMixAndMatchEnvTest);
TEST_COMPOSITE_ADD_UNIT(CoreSingularCellTest);
TEST_COMPOSITE_ADD_UNIT(CoreBasisDataTest);
TEST_COMPOSITE_ADD_UNIT(CoreNyquistGroupTest);
TEST_COMPOSITE_ADD_UNIT(CoreNyquistPlaneCacheTest);
#ifdef PSCF_CUDA
TEST_COMPOSITE_ADD_UNIT(GpuWaveListTest);
TEST_COMPOSITE_ADD_UNIT(GpuVecRandomTest);
TEST_COMPOSITE_ADD_UNIT(GpuFftTest);
TEST_COMPOSITE_ADD_UNIT(GpuKSqOverflowTest);
TEST_COMPOSITE_ADD_UNIT(GpuNyquistTest);
TEST_COMPOSITE_ADD_UNIT(GpuNyquistGroupTest);
#endif
TEST_COMPOSITE_ADD_UNIT(MinimumImageScaleTest);
#ifdef PSCF_CUDA
TEST_COMPOSITE_ADD_UNIT(MinimumImageScaleGpuTest);
#endif
TEST_COMPOSITE_ADD_UNIT(BeadCountValidationTest);
TEST_COMPOSITE_ADD_UNIT(FftSliceAlignmentTest);
TEST_COMPOSITE_ADD_UNIT(FieldReadEndOfFileTest);
TEST_COMPOSITE_ADD_UNIT(FourierColumnSpacingTest);
TEST_COMPOSITE_ADD_UNIT(BasisScaleTest);
TEST_COMPOSITE_ADD_UNIT(WaveShellScaleTest);
#ifdef PSCF_CUDA
TEST_COMPOSITE_ADD_UNIT(WaveShellScaleGpuTest);
#endif
TEST_COMPOSITE_ADD_UNIT(IdentityGroupReloadTest);
TEST_COMPOSITE_ADD_UNIT(AnisotropicMirrorMeshTest);
#ifdef PSCF_CUDA
TEST_COMPOSITE_ADD_UNIT(GpuComplexExpAliasTest);
#endif
#ifdef PSCF_CUDA
TEST_COMPOSITE_ADD_UNIT(GpuThreadMeshDimensionTest);
#endif
TEST_COMPOSITE_ADD_UNIT(EdgeArchiveStateTest);
TEST_COMPOSITE_ADD_UNIT(MonomerArchiveTest);
TEST_COMPOSITE_ADD_UNIT(UnitCellArchiveStateTest);
TEST_COMPOSITE_ADD_UNIT(FieldArchiveRoundTripTest);
#ifdef PSCF_CPP
TEST_COMPOSITE_ADD_UNIT(FftwComplexArrayCopyTest);
#endif
TEST_COMPOSITE_ADD_UNIT(EdgeOutputModeTest);
TEST_COMPOSITE_END

#endif
