#ifndef PSCF_WAVE_SHELL_SCALE_TEST_H
#define PSCF_WAVE_SHELL_SCALE_TEST_H

#include "WaveShellScaleChecks.h"
#include <test/UnitTestRunner.h>
#include <pscf/crystal/Basis.h>
#include <pscf/crystal/SpaceGroup.h>
#include <pscf/crystal/unitCellHeader.h>
#include <util/math/Constants.h>
#include <util/misc/Exception.h>
#include <set>
#include <sstream>

using namespace Util;
using namespace Pscf;

class WaveShellScaleTest : public UnitTest
{
public:
   void testCubicShells()
   {
      printMethod(TEST_FUNC);
      checkWaveShellScaleCubicShells<CPT>();
   }

   void testDistinctShells()
   {
      printMethod(TEST_FUNC);
      checkWaveShellScaleDistinctShells<CPT>();
   }

};

TEST_BEGIN(WaveShellScaleTest)
TEST_ADD(WaveShellScaleTest, testCubicShells)
TEST_ADD(WaveShellScaleTest, testDistinctShells)
TEST_END(WaveShellScaleTest)

#endif
