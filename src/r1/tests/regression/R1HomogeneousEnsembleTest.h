#ifndef R1_HOMOGENEOUS_ENSEMBLE_TEST_H
#define R1_HOMOGENEOUS_ENSEMBLE_TEST_H

#include "R1HomogeneousTest.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/field/Domain.h>
#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <r1/solvers/Solvent.h>
#include <pscf/interaction/Interaction.h>
#include <util/misc/ioUtil.h>

#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

/*
* COMPARE_HOMOGENEOUS in closed, mixed and grand-canonical ensembles,
* in modes 0, 1 and 2, checked against independent Flory-Huggins
* computations. In a mixed or grand-canonical ensemble, the species
* volume fractions of an SCFT state add up to 1 only within the SCFT
* error tolerance. Uses the helpers of R1HomogeneousTest.
*/
class R1HomogeneousEnsembleTest : public R1HomogeneousTest
{

public:

   /*
   * Read a System from parameter file text, with the file prefixes of
   * the unit test.
   */
   void readSystemText(System& sys, std::string const & text)
   {
      std::istringstream in(text);
      sys.readParam(in);
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());
   }

   /*
   * Species volume fractions of the current state of sys (polymers,
   * then solvents), and their sum.
   */
   double speciesPhi(System const & sys, std::vector<double>& phi)
   {
      int np = sys.mixture().nPolymer();
      int ns = sys.mixture().nSolvent();
      phi.resize(np + ns);
      double sum = 0.0;
      for (int i = 0; i < np; ++i) {
         phi[i] = sys.mixture().polymer(i).phi();
         sum += phi[i];
      }
      for (int i = 0; i < ns; ++i) {
         phi[np + i] = sys.mixture().solvent(i).phi();
         sum += phi[np + i];
      }
      return sum;
   }

   /*
   * Check the output of mode 0 against the Flory-Huggins mixture of two
   * monomer types with parameter chi whose species volume fractions
   * are those of the state of sys rescaled to a sum of 1. Species i
   * consists of size[i] monomers of type monomer[i]; a species with
   * monomer[i] < 0 (several clumps) must be absent. If checkDeltaF,
   * also check delta f against the SCFT free energy of sys.
   */
   void checkModeZero(std::string const & text, System const & sys,
                      int const * monomer, double const * size,
                      double chi, bool checkDeltaF)
   {
      std::vector<double> phi;
      double sum = speciesPhi(sys, phi);
      int nt = phi.size();
      double c[2] = {0.0, 0.0};
      for (int i = 0; i < nt; ++i) {
         phi[i] /= sum;
         if (monomer[i] < 0) {
            TEST_ASSERT(phi[i] == 0.0);
         } else {
            c[monomer[i]] += phi[i];
         }
      }

      // f = sum_i phi_i (log(phi_i) - 1)/size_i + chi cA cB,
      // mu_i = log(phi_i) + size_i w(monomer_i), wA = chi cB, wB = chi cA
      double w[2] = {chi*c[1], chi*c[0]};
      double f = chi*c[0]*c[1];
      for (int i = 0; i < nt; ++i) {
         if (phi[i] > 0.0) {
            f += phi[i]*(std::log(phi[i]) - 1.0)/size[i];
         }
      }
      Log::file() << std::setprecision(17)
                  << "Sum of SCFT species phi - 1 = " << sum - 1.0
                  << "\nExpected: f (homo) = " << f << std::endl;
      TEST_ASSERT(std::abs(labeledValue(text, "f (homo)") - f) < 1.0E-10);
      if (checkDeltaF) {
         double df = sys.fHelmholtz() - f;
         TEST_ASSERT(std::abs(labeledValue(text, "delta f") - df)
                     < 1.0E-10);
      }
      double phiHomo, muHomo, mu;
      for (int i = 0; i < nt; ++i) {
         phiHomo = speciesValue(text, i, 2);
         muHomo = speciesValue(text, i, 1);
         if (phi[i] > 0.0) {
            mu = std::log(phi[i]) + size[i]*w[monomer[i]];
            Log::file() << "Expected: species " << i << ", phi = "
                        << phi[i] << ", mu = " << mu << std::endl;
            TEST_ASSERT(std::abs(phiHomo/phi[i] - 1.0) < 1.0E-10);
            TEST_ASSERT(std::abs(muHomo - mu) < 1.0E-10);
         } else {
            TEST_ASSERT(phiHomo == 0.0);
            TEST_ASSERT(std::isinf(muHomo) && muHomo < 0.0);
         }
      }
   }

   /*
   * Volume fraction phiA of A in the homogeneous blend of A and B
   * homopolymers of length 1 with chi < 2, the unique root of
   *
   *    log(phiA/(1 - phiA)) + chi*(1 - 2 phiA) = dMu
   *
   * (the left side increases with phiA for chi < 2). Bisection in
   * t = log(phiA/(1 - phiA)).
   */
   double blendFraction(double chi, double dMu)
   {
      UTIL_CHECK(chi < 2.0);
      double lo = -700.0;
      double hi = 700.0;
      double mid, phiA;
      for (int k = 0; k < 200; ++k) {
         mid = 0.5*(lo + hi);
         phiA = 1.0/(1.0 + std::exp(-mid));
         if (mid + chi*(1.0 - 2.0*phiA) < dMu) {
            lo = mid;
         } else {
            hi = mid;
         }
      }
      return 1.0/(1.0 + std::exp(-0.5*(lo + hi)));
   }

   /*
   * Check the output of mode 1 or 2 for a one-phase blend (chi < 2) of
   * A and B homopolymers of length 1 (species 0 and 1, chemical
   * potentials muA and muB). If iS >= 0, species iS is a solvent of
   * size 1 of monomer type B with chemical potential muS: it is the
   * same molecule as B, so phiS = phiB exp(muS - muB) and the blend is
   * binary in A and cB = phiB + phiS. Species iAbsent (if >= 0) is
   * absent.
   */
   void checkUniformBlend(std::string const & text, double chi,
                          double muA, double muB, int iS, double muS,
                          int iAbsent)
   {
      double r = (iS >= 0) ? std::exp(muS - muB) : 0.0;
      double phiA = blendFraction(chi, muA - muB - std::log1p(r));
      double cB = 1.0 - phiA;
      double phiB = cB/(1.0 + r);
      double phiS = r*phiB;
      double f = phiA*(std::log(phiA) - 1.0) + phiB*(std::log(phiB) - 1.0)
               + chi*phiA*cB;
      double p = phiA*muA + phiB*muB;
      if (iS >= 0) {
         f += phiS*(std::log(phiS) - 1.0);
         p += phiS*muS;
      }
      p -= f;
      Log::file() << std::setprecision(17)
                  << "Expected: f (homo) = " << f << ", p (homo) = " << p
                  << ", phi(A) = " << phiA << std::endl;
      TEST_ASSERT(std::abs(labeledValue(text, "f (homo)") - f) < 1.0E-10);
      TEST_ASSERT(std::abs(labeledValue(text, "p (homo)") - p) < 1.0E-9);
      TEST_ASSERT(std::abs(speciesValue(text, 0, 2)/phiA - 1.0) < 1.0E-9);
      TEST_ASSERT(std::abs(speciesValue(text, 1, 2)/phiB - 1.0) < 1.0E-9);
      if (iS >= 0) {
         Log::file() << "Expected: phi(S) = " << phiS << std::endl;
         TEST_ASSERT(std::abs(speciesValue(text, iS, 2)/phiS - 1.0)
                     < 1.0E-9);
      }
      if (iAbsent >= 0) {
         TEST_ASSERT(speciesValue(text, iAbsent, 2) == 0.0);
      }
   }

   /*
   * Compare the state of sys in mode 0. Return the output, and set
   * logText to the text that the comparison wrote to logName, the log
   * file of the test (opened by openLogFile).
   */
   std::string compareModeZero(System& sys, char const * logName,
                               std::string& logText)
   {
      Log::file() << std::flush;
      std::size_t start = fileContents(logName).size();
      std::string text = compareHomogeneous(sys, 0);
      Log::file() << std::flush;
      logText = fileContents(logName).substr(start);
      return text;
   }

   /*
   * Check the log text of a mode 0 comparison: it gives the sum of the
   * SCFT volume fractions, and the warning about the rescaling appears
   * if and only if isWarned (|sum - 1| above 1e-6).
   */
   void checkModeZeroLog(std::string const & logText, double sum,
                         bool isWarned)
   {
      double logged = labeledValue(logText,
                                   "Sum of species volume fractions");
      TEST_ASSERT(std::abs(logged/sum - 1.0) < 1.0E-12);
      std::size_t pos = logText.find("Warning: Species volume fractions");
      TEST_ASSERT((pos != std::string::npos) == isWarned);
   }

   /*
   * Modes 0, 1 and 2 for an A/B interface (A at both ends, B in the
   * middle) at chi = 6 in the closed, mixed and grand-canonical
   * ensembles.
   */
   void testInterface()
   {
      printMethod(TEST_FUNC);
      char const * logName = "out/R1HomogeneousEnsembleTestInterface.log";
      openLogFile(logName);

      int monomer[2] = {0, 1};
      double size[2] = {1.0, 1.0};

      // Closed: A and B with phi = 0.5
      System closed;
      readSystem(closed, "in/R1/blend_ens_closed.prm");
      runCommands(closed, "READ_W in/R1/blend_chi20.w\nFINISH\n");
      TEST_ASSERT(closed.iterate() == 0);
      runCommands(closed, "WRITE_W out/blendEnsClosed.w\nFINISH\n");

      // Mixed: A closed with phi = 0.5, B open with mu = -6. Open: A and
      // B with mu = -6 (at coexistence of the two phases).
      System mixed;
      readSystem(mixed, "in/R1/blend_ens_mixed.prm");
      runCommands(mixed, "READ_W out/blendEnsClosed.w\nFINISH\n");
      TEST_ASSERT(mixed.iterate() == 0);
      System open;
      readSystem(open, "in/R1/blend_ens_open.prm");
      runCommands(open, "READ_W out/blendEnsClosed.w\nFINISH\n");
      TEST_ASSERT(open.iterate() == 0);

      System* systems[3] = {&closed, &mixed, &open};
      char const * names[3] = {"closed", "mixed", "grand-canonical"};
      std::string text, logText;
      std::vector<double> phi;
      double chi, muA, muB;
      int nx, ix;
      for (int k = 0; k < 3; ++k) {
         System& sys = *systems[k];
         Log::file() << "\nEnsemble: " << names[k] << std::endl;
         chi = sys.interaction().chi(0,1);

         // Mode 0: the blend with the average composition, without a
         // warning for these converged states
         text = compareModeZero(sys, logName, logText);
         TEST_ASSERT(!text.empty());
         checkModeZero(text, sys, monomer, size, chi, true);
         checkModeZeroLog(logText, speciesPhi(sys, phi), false);

         // Modes 1 and 2: the A-rich phase with the same mu
         muA = sys.mixture().polymer(0).mu();
         muB = sys.mixture().polymer(1).mu();
         nx = sys.domain().nx();
         for (int mode = 1; mode <= 2; ++mode) {
            text = compareHomogeneous(sys, mode);
            TEST_ASSERT(!text.empty());
            ix = (mode == 1) ? nx - 1 : 0;
            TEST_ASSERT(sys.cField(0)[ix] > 0.5);
            checkBinaryBlend(text, chi, muA, muB, true);
         }
      }
   }

   /*
   * Sweeps with homogeneousMode 0 in the mixed and grand-canonical
   * ensembles, which change mu of B (mixed) or of A and B (open).
   */
   void testInterfaceSweep()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1HomogeneousEnsembleTestInterfaceSweep.log");

      // Initial fields: closed solution
      System closed;
      readSystem(closed, "in/R1/blend_ens_closed.prm");
      runCommands(closed, "READ_W in/R1/blend_chi20.w\nFINISH\n");
      TEST_ASSERT(closed.iterate() == 0);
      runCommands(closed, "WRITE_W out/blendEnsSweep.w\nFINISH\n");

      int monomer[2] = {0, 1};
      double size[2] = {1.0, 1.0};
      char const * params[2] = {"in/R1/blend_ens_mixed_sweep.prm",
                                "in/R1/blend_ens_open_sweep.prm"};
      char const * prefixes[2] = {"out/blendEnsMixedSweep",
                                  "out/blendEnsOpenSweep"};
      std::string text, file;
      double chi, f, phi0, phi1;
      bool thrown;
      for (int k = 0; k < 2; ++k) {
         for (int i = 0; i <= 2; ++i) {
            file = std::string(prefixes[k]) + toString(i) + ".stt";
            std::remove((filePrefix() + file).c_str());
         }
         System sys;
         readSystem(sys, params[k]);
         thrown = false;
         try {
            runCommands(sys, "READ_W out/blendEnsSweep.w\nSWEEP\nFINISH\n");
         } catch (Exception& e) {
            thrown = true;
         }
         TEST_ASSERT(!thrown);
         chi = sys.interaction().chi(0,1);

         // Every state: phi(homo) adds up to 1 and gives f (homo)
         for (int i = 0; i <= 2; ++i) {
            file = std::string(prefixes[k]) + toString(i) + ".stt";
            text = fileContents(file.c_str());
            phi0 = speciesValue(text, 0, 2);
            phi1 = speciesValue(text, 1, 2);
            TEST_ASSERT(std::abs(phi0 + phi1 - 1.0) < 1.0E-10);
            f = phi0*(std::log(phi0) - 1.0) + phi1*(std::log(phi1) - 1.0)
              + chi*phi0*phi1;
            TEST_ASSERT(std::abs(labeledValue(text, "f (homo)") - f)
                        < 1.0E-9);
            TEST_ASSERT(std::isfinite(labeledValue(text, "delta f")));
         }

         // Final state, at full precision
         checkModeZero(text, sys, monomer, size, chi, true);
      }
   }

   /*
   * Modes 0, 1 and 2 for uniform states of a one-phase A/B blend
   * (chi = 1.5) in the closed, mixed and grand-canonical ensembles,
   * with chemical potentials that put the sum of the volume fractions
   * slightly off 1, as an SCFT tolerance does. Also with a dilute open
   * solvent and an absent closed diblock, and with a single present
   * species whose volume fraction exceeds 1 slightly. Mode 0 logs the
   * sum, and warns about the rescaling above |sum - 1| = 1e-6.
   */
   void testUniform()
   {
      printMethod(TEST_FUNC);
      char const * logName = "out/R1HomogeneousEnsembleTestUniform.log";
      openLogFile(logName);

      // Uniform fields of the Flory-Huggins solution with phiA = 0.3,
      // phiB = 0.7: wA = chi*phiB, wB = chi*phiA. For a homopolymer of
      // length 1, or a solvent of size 1, phi = exp(mu - w).
      double const chi = 1.5;
      double const wA = chi*0.7;
      double const wB = chi*0.3;
      std::string const head =
         "System{\nMixture{\nnMonomer 2\nmonomers 1.0 1.0\n";
      std::ostringstream tail;
      tail << "ds 0.01\n}\nInteraction{\nchi( 0 1 " << chi << " )\n}\n"
           << "Domain{\nmode Planar\nxMin 0.0\nxMax 1.0\nnx 5\n}\n"
           << "NrIterator{\nepsilon 1.0E-8\n}\n}\n";
      std::string const polyA =
         "Polymer{\ntype linear\nnBlock 1\nblocks 0 1.0\n";
      std::string const polyB =
         "Polymer{\ntype linear\nnBlock 1\nblocks 1 1.0\n";
      std::string const diblock =
         "Polymer{\ntype linear\nnBlock 2\nblocks 0 0.5\n1 0.5\nphi 0.0\n}\n";

      // Species: A, B; closed, mixed (A closed), open
      std::ostringstream specA[3], specB[3];
      specA[0] << "phi 0.3\n}\n";
      specB[0] << "phi 0.7\n}\n";
      specA[1] << "phi 0.3\n}\n";
      specB[1] << std::setprecision(17)
               << "mu " << std::log(0.7) + wB + 1.0E-8 << "\n}\n";
      specA[2] << std::setprecision(17)
               << "mu " << std::log(0.3) + wA + 3.0E-8 << "\n}\n";
      specB[2] << std::setprecision(17)
               << "mu " << std::log(0.7) + wB - 2.0E-8 << "\n}\n";
      char const * names[3] = {"closed", "mixed", "grand-canonical"};
      int monomer[4] = {0, 1, -1, 1};
      double size[4] = {1.0, 1.0, 1.0, 1.0};
      std::string text, logText;
      double sum, muA, muB, muS;
      std::vector<double> phi;
      int nx;
      for (int k = 0; k < 3; ++k) {
         Log::file() << "\nEnsemble: " << names[k] << std::endl;
         System sys;
         readSystemText(sys, head + "nPolymer 2\nnSolvent 0\n"
                        + polyA + specA[k].str() + polyB + specB[k].str()
                        + tail.str());
         nx = sys.domain().nx();
         for (int i = 0; i < nx; ++i) {
            sys.wField(0)[i] = wA;
            sys.wField(1)[i] = wB;
         }
         sys.compute();
         sum = speciesPhi(sys, phi);
         if (k > 0) {
            TEST_ASSERT(std::abs(sum - 1.0) > 1.0E-9);
            TEST_ASSERT(std::abs(sum - 1.0) < 1.0E-7);
         }
         text = compareModeZero(sys, logName, logText);
         TEST_ASSERT(!text.empty());
         checkModeZero(text, sys, monomer, size, chi, true);
         checkModeZeroLog(logText, sum, false);
         muA = sys.mixture().polymer(0).mu();
         muB = sys.mixture().polymer(1).mu();
         for (int mode = 1; mode <= 2; ++mode) {
            text = compareHomogeneous(sys, mode);
            TEST_ASSERT(!text.empty());
            checkUniformBlend(text, chi, muA, muB, -1, 0.0, -1);
         }
      }

      // Mixed: A and B open, a closed AB diblock with phi = 0, and an
      // open solvent of size 1 of B monomers with phi = 1E-20
      {
         Log::file() << "\nDilute open solvent, absent diblock"
                     << std::endl;
         std::ostringstream text0;
         text0 << std::setprecision(17)
               << head << "nPolymer 3\nnSolvent 1\n"
               << polyA << "mu " << std::log(0.3) + wA + 1.0E-8 << "\n}\n"
               << polyB << "mu " << std::log(0.7) + wB << "\n}\n"
               << diblock
               << "Solvent{\nmonomerId 1\nsize 1.0\nmu "
               << std::log(1.0E-20) + wB << "\n}\n"
               << tail.str();
         System sys;
         readSystemText(sys, text0.str());
         nx = sys.domain().nx();
         for (int i = 0; i < nx; ++i) {
            sys.wField(0)[i] = wA;
            sys.wField(1)[i] = wB;
         }
         sys.compute();
         sum = speciesPhi(sys, phi);
         TEST_ASSERT(std::abs(sum - 1.0) > 1.0E-9);
         TEST_ASSERT(std::abs(phi[3]/1.0E-20 - 1.0) < 1.0E-6);
         TEST_ASSERT(phi[2] == 0.0);
         text = compareModeZero(sys, logName, logText);
         TEST_ASSERT(!text.empty());
         checkModeZero(text, sys, monomer, size, chi, true);
         checkModeZeroLog(logText, sum, false);
         muA = sys.mixture().polymer(0).mu();
         muB = sys.mixture().polymer(1).mu();
         muS = sys.mixture().solvent(0).mu();
         for (int mode = 1; mode <= 2; ++mode) {
            text = compareHomogeneous(sys, mode);
            TEST_ASSERT(!text.empty());
            checkUniformBlend(text, chi, muA, muB, 3, muS, 2);
         }
      }

      // Mixed: A open with phi = 1 + 1E-9 (wA = 0), B closed with
      // phi = 0. The reference is pure A: f = -1, mu(A) = 0 in mode 0,
      // and p = 1 + mu(A) in modes 1 and 2.
      {
         Log::file() << "\nSingle present species" << std::endl;
         std::ostringstream text0;
         text0 << std::setprecision(17)
               << head << "nPolymer 2\nnSolvent 0\n"
               << polyA << "mu " << 1.0E-9 << "\n}\n"
               << polyB << "phi 0.0\n}\n"
               << tail.str();
         System sys;
         readSystemText(sys, text0.str());
         nx = sys.domain().nx();
         for (int i = 0; i < nx; ++i) {
            sys.wField(0)[i] = 0.0;
            sys.wField(1)[i] = chi;
         }
         sys.compute();
         sum = speciesPhi(sys, phi);
         TEST_ASSERT(phi[0] > 1.0);
         TEST_ASSERT(std::abs(sum - 1.0) < 1.0E-7);
         TEST_ASSERT(phi[1] == 0.0);
         text = compareModeZero(sys, logName, logText);
         TEST_ASSERT(!text.empty());
         checkModeZero(text, sys, monomer, size, chi, true);
         checkModeZeroLog(logText, sum, false);
         TEST_ASSERT(std::abs(labeledValue(text, "f (homo)") + 1.0)
                     < 1.0E-12);
         TEST_ASSERT(std::abs(speciesValue(text, 0, 1)) < 1.0E-12);
         muA = sys.mixture().polymer(0).mu();
         for (int mode = 1; mode <= 2; ++mode) {
            text = compareHomogeneous(sys, mode);
            TEST_ASSERT(!text.empty());
            TEST_ASSERT(std::abs(labeledValue(text, "f (homo)") + 1.0)
                        < 1.0E-12);
            TEST_ASSERT(std::abs(labeledValue(text, "p (homo)")
                                 - (1.0 + muA)) < 1.0E-12);
            TEST_ASSERT(std::abs(speciesValue(text, 0, 2) - 1.0)
                        < 1.0E-12);
            TEST_ASSERT(speciesValue(text, 1, 2) == 0.0);
         }
      }

      // Grand canonical, mu of A and B shifted by log(s): the volume
      // fractions are 0.3 s and 0.7 s. The reference composition is
      // (0.3, 0.7) for every s, and the log warns about the rescaling
      // only above |s - 1| = 1e-6. The sum s = 0.8 is that of a state
      // far from an SCFT solution, as after COMPUTE of a poor guess.
      double const sums[3] = {1.0 + 5.0E-7, 1.0 - 2.0E-6, 0.8};
      for (int k = 0; k < 3; ++k) {
         Log::file() << "\nGrand canonical, sum of phi = " << sums[k]
                     << std::endl;
         double shift = std::log(sums[k]);
         std::ostringstream text0;
         text0 << std::setprecision(17)
               << head << "nPolymer 2\nnSolvent 0\n"
               << polyA << "mu " << std::log(0.3) + wA + shift << "\n}\n"
               << polyB << "mu " << std::log(0.7) + wB + shift << "\n}\n"
               << tail.str();
         System sys;
         readSystemText(sys, text0.str());
         nx = sys.domain().nx();
         for (int i = 0; i < nx; ++i) {
            sys.wField(0)[i] = wA;
            sys.wField(1)[i] = wB;
         }
         sys.compute();
         sum = speciesPhi(sys, phi);
         TEST_ASSERT(std::abs(sum/sums[k] - 1.0) < 1.0E-12);
         text = compareModeZero(sys, logName, logText);
         TEST_ASSERT(!text.empty());
         checkModeZero(text, sys, monomer, size, chi, true);
         TEST_ASSERT(std::abs(speciesValue(text, 0, 2)/0.3 - 1.0) < 1.0E-10);
         TEST_ASSERT(std::abs(speciesValue(text, 1, 2)/0.7 - 1.0) < 1.0E-10);
         checkModeZeroLog(logText, sum, k > 0);
      }
   }

};

TEST_BEGIN(R1HomogeneousEnsembleTest)
TEST_ADD(R1HomogeneousEnsembleTest, testInterface)
TEST_ADD(R1HomogeneousEnsembleTest, testInterfaceSweep)
TEST_ADD(R1HomogeneousEnsembleTest, testUniform)
TEST_END(R1HomogeneousEnsembleTest)

#endif
