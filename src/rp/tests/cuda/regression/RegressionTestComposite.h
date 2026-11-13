#ifndef RPG_REGRESSION_TEST_COMPOSITE_H
#define RPG_REGRESSION_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "ScftAmIteratorTest.h"
#include "FieldIoRegressionTest.h"
#include "FilmStressTest.h"
#include "FilmMaskCheckTest.h"
#include "FilmCpuGpuTest.h"
#include "FilmSweepTest.h"
#include "ScftSweepTest.h"
#include "FieldIoHeaderTest.h"
#include "FilmSemiGrandTest.h"
#include "CommandTest.h"
#include "FtsAnalyzerTest.h"
#include "FtsBdSimulatorTest.h"
#include "FtsCompressorTest.h"
#include "FtsMcSimulatorTest.h"
#include "FtsPerturbationTest.h"
#include "FtsStructureFactorTest.h"
#include "FtsAnalyzePrefixTest.h"
#include "ScftThermoRefTest.h"
#include "CommandReadHTest.h"
#include "ScftScaleStressTest.h"
#include "FieldIoStarsHeaderTest.h"
#include "CommandExpandTest.h"
#include "FilmOrthogonalityTest.h"
#include "FilmSpaceGroupTest.h"
#include "FilmCellChangeTest.h"
#include "FilmWriteFieldsTest.h"
#include "FilmParamStreamTest.h"
#include "ScftLinearSweepTest.h"
#include "FtsLinearRampTest.h"
#include "FtsBdMoveFailTest.h"
#include "FtsLmBdStepSplitTest.h"
#include "FtsRelNormTest.h"
#include "ScftLongChainTest.h"
#include "ScftContourStepTest.h"
#include "ScftNyquistStressTest.h"
#include "CommandStatusTest.h"
#include "FtsBdMoveNStepTest.h"
#include "FtsAmCompressorLogTest.h"
#include "CommandReplicateTest.h"
#include "ScftBeadDsTest.h"
#include "ScftNyquistGroupTest.h"
#include "FtsLrNyquistTest.h"
#include "FtsLrResponseTest.h"

#include "DiluteThermodynamicsTest.h"

#include "SimultaneousChiRampTest.h"

#include "RampParameterRoundTripTest.h"

#include "InteractionEnergyScaleTest.h"

#include "HomogeneousClumpScaleTest.h"

TEST_COMPOSITE_BEGIN(RegressionTestComposite)
TEST_COMPOSITE_ADD_UNIT(ScftAmIteratorTest);
TEST_COMPOSITE_ADD_UNIT(FieldIoRegressionTest);
TEST_COMPOSITE_ADD_UNIT(FilmStressTest);
TEST_COMPOSITE_ADD_UNIT(FilmMaskCheckTest);
TEST_COMPOSITE_ADD_UNIT(FilmCpuGpuTest);
TEST_COMPOSITE_ADD_UNIT(FilmSweepTest);
TEST_COMPOSITE_ADD_UNIT(ScftSweepTest);
TEST_COMPOSITE_ADD_UNIT(FieldIoHeaderTest);
TEST_COMPOSITE_ADD_UNIT(FilmSemiGrandTest);
TEST_COMPOSITE_ADD_UNIT(CommandTest);
TEST_COMPOSITE_ADD_UNIT(FtsAnalyzerTest);
TEST_COMPOSITE_ADD_UNIT(FtsBdSimulatorTest);
TEST_COMPOSITE_ADD_UNIT(FtsCompressorTest);
TEST_COMPOSITE_ADD_UNIT(FtsMcSimulatorTest);
TEST_COMPOSITE_ADD_UNIT(FtsPerturbationTest);
TEST_COMPOSITE_ADD_UNIT(FtsStructureFactorTest);
TEST_COMPOSITE_ADD_UNIT(FtsAnalyzePrefixTest);
TEST_COMPOSITE_ADD_UNIT(ScftThermoRefTest);
TEST_COMPOSITE_ADD_UNIT(CommandReadHTest);
TEST_COMPOSITE_ADD_UNIT(ScftScaleStressTest);
TEST_COMPOSITE_ADD_UNIT(FieldIoStarsHeaderTest);
TEST_COMPOSITE_ADD_UNIT(CommandExpandTest);
TEST_COMPOSITE_ADD_UNIT(FilmOrthogonalityTest);
TEST_COMPOSITE_ADD_UNIT(FilmSpaceGroupTest);
TEST_COMPOSITE_ADD_UNIT(FilmCellChangeTest);
TEST_COMPOSITE_ADD_UNIT(FilmWriteFieldsTest);
TEST_COMPOSITE_ADD_UNIT(FilmParamStreamTest);
TEST_COMPOSITE_ADD_UNIT(ScftLinearSweepTest);
TEST_COMPOSITE_ADD_UNIT(FtsLinearRampTest);
TEST_COMPOSITE_ADD_UNIT(FtsBdMoveFailTest);
TEST_COMPOSITE_ADD_UNIT(FtsLmBdStepSplitTest);
TEST_COMPOSITE_ADD_UNIT(FtsRelNormTest);
TEST_COMPOSITE_ADD_UNIT(ScftLongChainTest);
TEST_COMPOSITE_ADD_UNIT(ScftContourStepTest);
TEST_COMPOSITE_ADD_UNIT(ScftNyquistStressTest);
TEST_COMPOSITE_ADD_UNIT(CommandStatusTest);
TEST_COMPOSITE_ADD_UNIT(FtsBdMoveNStepTest);
TEST_COMPOSITE_ADD_UNIT(FtsAmCompressorLogTest);
TEST_COMPOSITE_ADD_UNIT(CommandReplicateTest);
TEST_COMPOSITE_ADD_UNIT(ScftBeadDsTest);
TEST_COMPOSITE_ADD_UNIT(ScftNyquistGroupTest);
TEST_COMPOSITE_ADD_UNIT(FtsLrNyquistTest);
TEST_COMPOSITE_ADD_UNIT(FtsLrResponseTest);
TEST_COMPOSITE_ADD_UNIT(DiluteThermodynamicsTest);
TEST_COMPOSITE_ADD_UNIT(SimultaneousChiRampTest);
TEST_COMPOSITE_ADD_UNIT(RampParameterRoundTripTest);
TEST_COMPOSITE_ADD_UNIT(InteractionEnergyScaleTest);
TEST_COMPOSITE_ADD_UNIT(HomogeneousClumpScaleTest);
TEST_COMPOSITE_END

#endif
