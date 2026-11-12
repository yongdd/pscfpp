#ifndef R1_HOMOGENEOUS_LIMIT_TEST_H
#define R1_HOMOGENEOUS_LIMIT_TEST_H

#include "R1HomogeneousTest.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/field/Domain.h>
#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <r1/solvers/Solvent.h>
#include <pscf/interaction/Interaction.h>

#include <string>
#include <cmath>
#include <cstdio>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

/*
* COMPARE_HOMOGENEOUS for extremely dilute species and for volume
* fractions whose sum differs from 1 by rounding. Uses the helpers of
* R1HomogeneousTest.
*/
class R1HomogeneousLimitTest : public R1HomogeneousTest
{

public:

   /*
   * Check a homogeneous comparison in mode 1 or 2 for a blend of A
   * and B homopolymers of length 1 (vMonomer = 1), species iA and iB,
   * against the phase rich in A (isARich) or B with the chemical
   * potentials muA and muB. Return the minority volume fraction of
   * the output.
   */
   double checkBlendPhase(std::string const & text, double chi,
                          int iA, int iB, double muA, double muB,
                          bool isARich)
   {
      double e;    // minority volume fraction
      if (isARich) {
         e = minorityFraction(chi, muA - muB);
      } else {
         e = minorityFraction(chi, muB - muA);
      }
      double phiA = isARich ? 1.0 - e : e;
      double phiB = isARich ? e : 1.0 - e;

      // f = sum_i phi_i (log(phi_i) - 1) + chi phiA phiB
      double f = (1.0 - e)*(std::log1p(-e) - 1.0)
               + e*(std::log(e) - 1.0) + chi*e*(1.0 - e);
      double p = -f + phiA*muA + phiB*muB;
      double eHomo = speciesValue(text, isARich ? iB : iA, 2);
      Log::file() << "Expected: f (homo) = " << f
                  << ", p (homo) = " << p
                  << ", minority phi = " << e << std::endl;
      TEST_ASSERT(std::abs(labeledValue(text, "f (homo)") - f) < 1.0E-10);
      TEST_ASSERT(std::abs(labeledValue(text, "p (homo)") - p) < 1.0E-8);
      TEST_ASSERT(std::abs(eHomo/e - 1.0) < 1.0E-8);
      return eHomo;
   }

   /*
   * Modes 1 and 2 for an A/B interface at chi = 44, in both orders of
   * the species.
   */
   void testVeryDilute()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1HomogeneousLimitTestVeryDilute.log");

      // Minority fraction of the bulk phases about 8E-20. Reference
      // values: pscfpp devel (no guess floor) for the same fields,
      // in mode 1 (B-rich phase at xMax) and mode 2 (A-rich at xMin)
      double ref[2] = {7.78113380406E-20, 7.78113067821E-20};

      char const * paramFiles[2] = {"in/R1/blend_chi44_ab.prm",
                                    "in/R1/blend_chi44_ba.prm"};
      for (int order = 0; order < 2; ++order) {

         // Species of the A and B homopolymers
         int iA = (order == 0) ? 0 : 1;
         int iB = 1 - iA;

         System sys;
         readSystem(sys, paramFiles[order]);
         runCommands(sys, "READ_W in/R1/blend_chi44.w\nCOMPUTE\nFINISH\n");
         double chi = sys.interaction().chi(0,1);
         double muA = sys.mixture().polymer(iA).mu();
         double muB = sys.mixture().polymer(iB).mu();
         int nx = sys.domain().nx();

         std::string text;
         int ix;
         double eHomo;
         for (int mode = 1; mode <= 2; ++mode) {
            text = compareHomogeneous(sys, mode);
            TEST_ASSERT(!text.empty());
            ix = (mode == 1) ? nx - 1 : 0;
            eHomo = checkBlendPhase(text, chi, iA, iB, muA, muB,
                                    sys.cField(0)[ix] > 0.5);
            TEST_ASSERT(std::abs(eHomo/ref[mode-1] - 1.0) < 1.0E-8);
         }
      }
   }

   /*
   * Sweep of chi from 44 to 50 with homogeneousMode 2.
   */
   void testVeryDiluteSweep()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1HomogeneousLimitTestVeryDiluteSweep.log");

      // B homopolymer is species 0, A is species 1
      std::remove((filePrefix() + "out/blendChi44Sweep3.stt").c_str());
      System sys;
      readSystem(sys, "in/R1/blend_chi44_sweep.prm");
      bool thrown = false;
      try {
         runCommands(sys, "READ_W in/R1/blend_chi44.w\nSWEEP\nFINISH\n");
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
      double chi = sys.interaction().chi(0,1);
      TEST_ASSERT(eq(chi, 50.0));

      // Final state: phase at xMin with the mu of the final solution
      std::string text = fileContents("out/blendChi44Sweep3.stt");
      checkBlendPhase(text, chi, 1, 0,
                      sys.mixture().polymer(1).mu(),
                      sys.mixture().polymer(0).mu(),
                      sys.cField(0)[0] > 0.5);
   }

   /*
   * Mode 0 when the fractions of the species other than the last add
   * up to 1 + 1 ulp in floating point, and the last species is absent.
   */
   void testRoundedSum()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1HomogeneousLimitTestRoundedSum.log");

      // AB diblock 0.56, A 0.33, B 0.11 (sum 1 + 2.2E-16), solvent 0
      System sys;
      readSystem(sys, "in/R1/blend_phi_sum.prm");
      runCommands(sys, "READ_W in/R1/planar.w\nCOMPUTE\nFINISH\n");
      std::string text = compareHomogeneous(sys, 0);
      TEST_ASSERT(!text.empty());

      // f = sum_i phi_i (log(phi_i) - 1) + chi cA cB, cA = 0.61
      double chi = sys.interaction().chi(0,1);
      double cA = 0.5*0.56 + 0.33;
      double f = 0.56*(std::log(0.56) - 1.0) + 0.33*(std::log(0.33) - 1.0)
               + 0.11*(std::log(0.11) - 1.0) + chi*cA*(1.0 - cA);
      TEST_ASSERT(std::abs(labeledValue(text, "f (homo)") - f) < 1.0E-10);
      TEST_ASSERT(std::isfinite(labeledValue(text, "delta f")));
      TEST_ASSERT(speciesValue(text, 3, 2) == 0.0);

      // Sweep with homogeneousMode 0 from A 0.01, B 0.56, AB 0.43 to
      // A 0.1, B 0.9, AB 0 (A + B = 1 + 1 ulp), in 4 steps
      std::remove((filePrefix() + "out/blendPhiSumSweep4.stt").c_str());
      System sweepSys;
      readSystem(sweepSys, "in/R1/blend_phi_sum_sweep.prm");
      int nx = sweepSys.domain().nx();
      for (int i = 0; i < nx; ++i) {
         sweepSys.wField(0)[i] = 0.0;
         sweepSys.wField(1)[i] = 0.0;
      }
      bool thrown = false;
      try {
         runCommands(sweepSys, "SWEEP\nFINISH\n");
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
      TEST_ASSERT(eq(sweepSys.mixture().polymer(2).phi(), 0.0));

      // Final state: binary blend A 0.1, B 0.9
      double phiA = sweepSys.mixture().polymer(0).phi();
      double phiB = sweepSys.mixture().polymer(1).phi();
      chi = sweepSys.interaction().chi(0,1);
      f = phiA*(std::log(phiA) - 1.0) + phiB*(std::log(phiB) - 1.0)
        + chi*phiA*phiB;
      text = fileContents("out/blendPhiSumSweep4.stt");
      TEST_ASSERT(std::abs(labeledValue(text, "f (homo)") - f) < 1.0E-10);
      TEST_ASSERT(std::isfinite(labeledValue(text, "delta f")));
   }

};

TEST_BEGIN(R1HomogeneousLimitTest)
TEST_ADD(R1HomogeneousLimitTest, testVeryDilute)
TEST_ADD(R1HomogeneousLimitTest, testVeryDiluteSweep)
TEST_ADD(R1HomogeneousLimitTest, testRoundedSum)
TEST_END(R1HomogeneousLimitTest)

#endif
