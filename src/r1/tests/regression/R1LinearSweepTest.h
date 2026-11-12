#ifndef R1_LINEAR_SWEEP_TEST_H
#define R1_LINEAR_SWEEP_TEST_H

#include "R1TestBase.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <util/misc/Exception.h>
#include <util/misc/FileMaster.h>

#include <cmath>
#include <sstream>
#include <string>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

/*
* Regression tests (pscfpp fixes): LinearSweep parameters and states.
*/
class R1LinearSweepTest : public R1TestBase
{

public:

};

TEST_BEGIN(R1LinearSweepTest)
TEST_ADD(R1LinearSweepTest, testIdsOutOfRange)
TEST_ADD(R1LinearSweepTest, testEndsAtOne)
TEST_END(R1LinearSweepTest)

#endif
