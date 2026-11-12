#ifndef R1_HOMOGENEOUS_TEST_H
#define R1_HOMOGENEOUS_TEST_H

#include "R1TestBase.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/system/HomogeneousComparison.h>
#include <r1/field/Domain.h>
#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <pscf/interaction/Interaction.h>

#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstdio>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

class R1HomogeneousTest : public R1TestBase
{

public:

   /*
   * Compare the current state with a homogeneous mixture. Return the
   * output, or an empty string if an Exception was thrown.
   */
   std::string compareHomogeneous(System& sys, int mode)
   {
      HomogeneousComparison comparison(sys);
      std::ostringstream out;
      try {
         comparison.compute(mode);
         comparison.output(mode, out);
      } catch (Exception& e) {
         return std::string();
      }
      Log::file() << "\nHomogeneous comparison, mode " << mode << "\n"
                  << out.str();
      return out.str();
   }

   /*
   * Return the number after "label =" in text (nan and inf allowed).
   */
   double labeledValue(std::string const & text, std::string const & label)
   {
      std::size_t pos = text.find(label);
      TEST_ASSERT(pos != std::string::npos);
      pos = text.find('=', pos);
      TEST_ASSERT(pos != std::string::npos);
      std::istringstream in(text.substr(pos + 1));
      std::string word;
      in >> word;
      return std::strtod(word.c_str(), 0);
   }

   /*
   * Return column 1 (mu) or 2 (phi) of the row for species i in the
   * species table of a homogeneous comparison.
   */
   double speciesValue(std::string const & text, int i, int column)
   {
      std::size_t pos = text.find("Species:");
      TEST_ASSERT(pos != std::string::npos);
      std::istringstream in(text.substr(pos));
      std::string line;
      std::getline(in, line);   // "Species:"
      std::getline(in, line);   // column labels
      for (int j = 0; j <= i; ++j) {
         std::getline(in, line);
      }
      std::istringstream row(line);
      int id;
      std::string word;
      row >> id;
      TEST_ASSERT(id == i);
      for (int j = 0; j < column; ++j) {
         row >> word;
      }
      return std::strtod(word.c_str(), 0);
   }

   /*
   * Minority volume fraction e of a homogeneous binary blend of
   * homopolymers of length 1 (vMonomer = 1) on the minority branch
   * e < e(spinodal), for chi > 2 and the difference dMu of the
   * chemical potentials of the majority and minority species:
   *
   *    log((1 - e)/e) - chi*(1 - 2e) = dMu
   *
   * The left side decreases with e on this branch (bisection in log e).
   */
   double minorityFraction(double chi, double dMu)
   {
      double lo = std::log(1.0E-300);
      double hi = std::log(0.5*(1.0 - std::sqrt(1.0 - 2.0/chi)));
      double mid, e;
      for (int k = 0; k < 200; ++k) {
         mid = 0.5*(lo + hi);
         e = std::exp(mid);
         if (std::log((1.0 - e)/e) - chi*(1.0 - 2.0*e) > dMu) {
            lo = mid;
         } else {
            hi = mid;
         }
      }
      return std::exp(0.5*(lo + hi));
   }

   /*
   * Check a homogeneous comparison in mode 1 or 2 for a blend of A
   * and B homopolymers of length 1 (species 0 and 1, vMonomer = 1)
   * against the homogeneous phase rich in A (isARich) or B with the
   * chemical potentials muA and muB.
   */
   void checkBinaryBlend(std::string const & text, double chi,
                         double muA, double muB, bool isARich)
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
      double eHomo = speciesValue(text, isARich ? 1 : 0, 2);
      Log::file() << "Expected: f (homo) = " << f
                  << ", p (homo) = " << p
                  << ", minority phi = " << e << std::endl;
      TEST_ASSERT(std::abs(labeledValue(text, "f (homo)") - f) < 1.0E-10);
      TEST_ASSERT(std::abs(labeledValue(text, "p (homo)") - p) < 1.0E-8);
      TEST_ASSERT(std::abs(eHomo/e - 1.0) < 1.0E-6);
   }

   /*
   * COMPARE_HOMOGENEOUS in modes 1 and 2 for a strongly segregated
   * blend.
   */
   void testSegregated()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1HomogeneousTestSegregated.log");

      // A and B homopolymers at chi = 20, fields of a converged
      // solution: nearly pure A at both ends (B fraction about 2E-9),
      // B in the middle
      System sys;
      readSystem(sys, "in/R1/blend_chi20.prm");
      runCommands(sys, "READ_W in/R1/blend_chi20.w\nCOMPUTE\nFINISH\n");
      int nx = sys.domain().nx();
      Log::file() << "c(A) - 1 at xMin, xMax: "
                  << sys.cField(0)[0] - 1.0 << ", "
                  << sys.cField(0)[nx-1] - 1.0 << std::endl;

      // Mode 1 starts from the composition at xMax, mode 2 from that
      // at xMin. Both must give the A-rich phase.
      double chi = sys.interaction().chi(0,1);
      double muA = sys.mixture().polymer(0).mu();
      double muB = sys.mixture().polymer(1).mu();
      std::string text;
      for (int mode = 1; mode <= 2; ++mode) {
         text = compareHomogeneous(sys, mode);
         TEST_ASSERT(!text.empty());
         checkBinaryBlend(text, chi, muA, muB, true);
      }
   }

   /*
   * COMPARE_HOMOGENEOUS with a closed species of volume fraction 0.
   */
   void testPhiZero()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1HomogeneousTestPhiZero.log");

      // A and B homopolymers with phi = 0.5 and an AB diblock with
      // phi = 0, at chi = 6
      System sys;
      readSystem(sys, "in/R1/blend_phi0.prm");
      double chi = sys.interaction().chi(0,1);

      // Mode 0: the binary blend with phi = 0.5, in which
      // f = log(0.5) - 1 + chi/4 and mu = log(0.5) + chi/2
      std::string text = compareHomogeneous(sys, 0);
      TEST_ASSERT(!text.empty());
      double f = std::log(0.5) - 1.0 + 0.25*chi;
      double mu = std::log(0.5) + 0.5*chi;
      TEST_ASSERT(std::abs(labeledValue(text, "f (homo)") - f) < 1.0E-10);
      TEST_ASSERT(std::abs(speciesValue(text, 0, 1) - mu) < 1.0E-10);
      TEST_ASSERT(std::abs(speciesValue(text, 1, 1) - mu) < 1.0E-10);
      TEST_ASSERT(speciesValue(text, 2, 2) == 0.0);

      // SCFT solution, from the fields of the chi = 20 blend
      runCommands(sys, "READ_W in/R1/blend_chi20.w\nFINISH\n");
      TEST_ASSERT(sys.iterate() == 0);

      // Modes 1 and 2: the binary blend phase with the same mu of A
      // and B as the SCFT solution, without the diblock
      double muA = sys.mixture().polymer(0).mu();
      double muB = sys.mixture().polymer(1).mu();
      int nx = sys.domain().nx();
      int ix;
      for (int mode = 1; mode <= 2; ++mode) {
         text = compareHomogeneous(sys, mode);
         TEST_ASSERT(!text.empty());
         ix = (mode == 1) ? nx - 1 : 0;
         checkBinaryBlend(text, chi, muA, muB, sys.cField(0)[ix] > 0.5);
         TEST_ASSERT(speciesValue(text, 2, 2) == 0.0);
      }

      // Sweep with homogeneousMode 0 in which the diblock volume
      // fraction goes from 0.1 to 0 (2 steps), from this solution
      runCommands(sys, "WRITE_W out/blendPhi0.w\nFINISH\n");
      std::remove((filePrefix() + "out/blendPhi0Sweep2.stt").c_str());
      System sweepSys;
      readSystem(sweepSys, "in/R1/blend_phi0_sweep.prm");
      runCommands(sweepSys, "READ_W out/blendPhi0.w\nSWEEP\nFINISH\n");
      TEST_ASSERT(eq(sweepSys.mixture().polymer(2).phi(), 0.0));

      // The state file of the final state gives the binary blend
      text = fileContents("out/blendPhi0Sweep2.stt");
      TEST_ASSERT(std::abs(labeledValue(text, "f (homo)") - f) < 1.0E-10);
      TEST_ASSERT(std::isfinite(labeledValue(text, "delta f")));
   }

};

TEST_BEGIN(R1HomogeneousTest)
TEST_ADD(R1HomogeneousTest, testSegregated)
TEST_ADD(R1HomogeneousTest, testPhiZero)
TEST_END(R1HomogeneousTest)

#endif
