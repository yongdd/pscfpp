#ifndef R1_BLOCK_TEST_H
#define R1_BLOCK_TEST_H

#include "R1TestBase.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/field/Domain.h>
#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <r1/solvers/Block.h>
#include <util/math/Constants.h>

#include <cmath>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

class R1BlockTest : public R1TestBase
{

public:

   /*
   * Solve the MDE for a block shorter than ds/2.
   */
   void testShortBlock()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1BlockTestShortBlock.log");

      // ABA triblock with blocks 0.5, 0.004 and 0.5, ds = 0.01
      System sys;
      readSystem(sys, "in/R1/planar_short.prm");

      Mixture& mix = sys.mixture();
      Domain& domain = sys.domain();
      Polymer& polymer = mix.polymer(0);

      int nx = domain.nx();
      double cs;
      for (int i = 0; i < nx; ++i) {
         cs = cos(2.0*Constants::Pi*double(i)/double(nx-1));
         sys.wField(0)[i] = 0.5 + cs;
         sys.wField(1)[i] = 0.7 - cs;
      }

      // Middle block of length 0.004 from the parameter file, then of
      // length 0.002 set by setLength (as in a sweep)
      double lengths[2] = {0.004, 0.002};
      for (int k = 0; k < 2; ++k) {
         if (k > 0) {
            polymer.block(1).setLength(lengths[k]);
         }
         TEST_ASSERT(polymer.block(1).ns() >= 3);

         mix.compute(sys.wFields(), sys.cFields());
         for (int i = 0; i < nx; ++i) {
            TEST_ASSERT(std::isfinite(sys.cField(0)[i]));
            TEST_ASSERT(std::isfinite(sys.cField(1)[i]));
         }

         // Monomer type 1 (the middle block) takes up its fraction
         // of the chain length
         double phi0 = domain.spatialAverage(sys.cField(0));
         double phi1 = domain.spatialAverage(sys.cField(1));
         double fraction = lengths[k]/polymer.length();
         Log::file() << "length = " << lengths[k]
                     << "  ns = " << polymer.block(1).ns()
                     << "  phi0 = " << phi0
                     << "  phi1 = " << phi1
                     << "  length fraction = " << fraction << std::endl;
         TEST_ASSERT(std::abs(phi1/fraction - 1.0) < 1.0E-6);
         TEST_ASSERT(std::abs(phi0 + phi1 - 1.0) < 1.0E-6);
      }
   }

};

TEST_BEGIN(R1BlockTest)
TEST_ADD(R1BlockTest, testShortBlock)
TEST_END(R1BlockTest)

#endif
