#ifndef RPC_SCFT_LINEAR_SWEEP_TEST_H
#define RPC_SCFT_LINEAR_SWEEP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <pscf/crystal/UnitCell.h>

#include <util/misc/Exception.h>
#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): LinearSweep parameters and states.
*/
class ScftLinearSweepTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

};

TEST_BEGIN(ScftLinearSweepTest)
TEST_ADD(ScftLinearSweepTest, testIdsOutOfRange)
TEST_ADD(ScftLinearSweepTest, testEndsAtOne)
TEST_END(ScftLinearSweepTest)

#endif
