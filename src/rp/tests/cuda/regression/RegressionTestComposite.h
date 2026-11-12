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
TEST_COMPOSITE_END

#endif
