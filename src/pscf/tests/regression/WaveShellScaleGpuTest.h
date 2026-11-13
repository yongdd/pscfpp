#ifndef PSCF_WAVE_SHELL_SCALE_GPU_TEST_H
#define PSCF_WAVE_SHELL_SCALE_GPU_TEST_H

#include "WaveShellScaleChecks.h"
#include <test/UnitTestRunner.h>

using namespace Pscf;

class WaveShellScaleGpuTest : public UnitTest
{
public:
   void testCubicShells()
   {
      printMethod(TEST_FUNC);
      checkWaveShellScaleCubicShells<CUT>();
   }

   void testDistinctShells()
   {
      printMethod(TEST_FUNC);
      checkWaveShellScaleDistinctShells<CUT>();
   }
};

TEST_BEGIN(WaveShellScaleGpuTest)
TEST_ADD(WaveShellScaleGpuTest, testCubicShells)
TEST_ADD(WaveShellScaleGpuTest, testDistinctShells)
TEST_END(WaveShellScaleGpuTest)

#endif
