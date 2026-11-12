#ifndef R1_GAUGE_TEST_H
#define R1_GAUGE_TEST_H

#include "R1TestBase.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/field/Domain.h>
#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <r1/solvers/Solvent.h>

#include <cmath>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

/*
* Independence of r1 results from the arbitrary constant of w, and long
* chains whose partition function is not representable as a double.
*/
class R1GaugeTest : public R1TestBase
{

public:

   /*
   * Add a constant to the w fields of all monomer types.
   */
   void shiftFields(System& sys, double shift)
   {
      int nm = sys.mixture().nMonomer();
      int nx = sys.domain().nx();
      for (int m = 0; m < nm; ++m) {
         for (int i = 0; i < nx; ++i) {
            sys.wField(m)[i] += shift;
         }
      }
   }

   /*
   * Set the w field of each monomer type to a constant.
   */
   void setUniformFields(System& sys, double w0, double w1)
   {
      int nx = sys.domain().nx();
      for (int i = 0; i < nx; ++i) {
         sys.wField(0)[i] = w0;
         sys.wField(1)[i] = w1;
      }
   }

   /*
   * A constant added to all w fields of a closed system changes neither
   * the free energy nor the concentrations, and shifts mu by K times
   * the chain length (examples/1d/planar/1, ds 0.01).
   */
   void testShiftInvariance()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1GaugeTestShiftInvariance.log");

      System sys;
      readSystem(sys, "in/R1/gauge_diblock.prm");
      runCommands(sys, "READ_W in/R1/gauge_diblock.w\nFINISH\n");
      int nx = sys.domain().nx();
      double length = sys.mixture().polymer(0).length();

      sys.compute();
      double f0 = sys.fHelmholtz();
      double p0 = sys.pressure();
      double mu0 = sys.mixture().polymer(0).mu();
      DArray<double> c0;
      c0.allocate(nx);
      for (int i = 0; i < nx; ++i) {
         c0[i] = sys.cField(0)[i];
      }

      double shifts[3] = {-20.0, 10.0, 30.0};
      double total = 0.0;
      double f, mu, p, cDiff;
      for (int k = 0; k < 3; ++k) {
         shiftFields(sys, shifts[k]);
         total += shifts[k];
         sys.compute();
         f = sys.fHelmholtz();
         p = sys.pressure();
         mu = sys.mixture().polymer(0).mu();
         cDiff = 0.0;
         for (int i = 0; i < nx; ++i) {
            cDiff = std::max(cDiff, std::abs(sys.cField(0)[i] - c0[i]));
         }
         Log::file() << "shift = " << total
                     << "  f - f0 = " << f - f0
                     << "  mu - mu0 - shift*N = "
                     << mu - mu0 - total*length
                     << "  p - p0 - shift = " << p - p0 - total
                     << "  max |c - c0| = " << cDiff << std::endl;
         TEST_ASSERT(std::abs(f - f0) < 1.0E-10);
         TEST_ASSERT(std::abs(mu - mu0 - total*length) < 1.0E-10);
         TEST_ASSERT(std::abs(p - p0 - total) < 1.0E-10);
         TEST_ASSERT(cDiff < 1.0E-10);
      }
   }

   /*
   * The same diblock with monomer labels swapped converges to the same
   * free energy, although the iterator pins the w field of a different
   * monomer type at the boundary.
   */
   void testLabelSwap()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1GaugeTestLabelSwap.log");

      System sys;
      readSystem(sys, "in/R1/gauge_diblock.prm");
      runCommands(sys, "READ_W in/R1/gauge_diblock.w\nITERATE\nFINISH\n");
      double f = sys.fHelmholtz();

      System swap;
      readSystem(swap, "in/R1/gauge_diblock_swap.prm");
      runCommands(swap,
                  "READ_W in/R1/gauge_diblock_swap.w\nITERATE\nFINISH\n");
      double fSwap = swap.fHelmholtz();

      Log::file() << "f = " << f << "  f (swapped labels) = " << fSwap
                  << "  difference = " << f - fSwap << std::endl;
      TEST_ASSERT(std::abs(f - fSwap) < 1.0E-7);
   }

   /*
   * Exactly homogeneous states of a diblock (chi N = 10) in uniform
   * fields with an arbitrary constant: f = chi fA fB - 1 exactly.
   */
   void testHomogeneous()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1GaugeTestHomogeneous.log");

      char const * paramFiles[2] = {"in/R1/gauge_homo.prm",
                                    "in/R1/gauge_homo_asym.prm"};
      double fA[2] = {0.5, 0.3};
      double shifts[4] = {0.0, 5.0, 20.0, -20.0};
      for (int j = 0; j < 2; ++j) {
         System sys;
         readSystem(sys, paramFiles[j]);
         double chi = sys.interaction().chi(0,1);
         double fB = 1.0 - fA[j];
         double fExact = chi*fA[j]*fB - 1.0;
         for (int k = 0; k < 4; ++k) {
            // Saddle point fields chi*c of the other type, plus K
            setUniformFields(sys, chi*fB + shifts[k], chi*fA[j] + shifts[k]);
            sys.compute();
            double f = sys.fHelmholtz();
            double mu = sys.mixture().polymer(0).mu();
            double muExact = chi*(fA[j]*fB + fB*fA[j]) + shifts[k];
            Log::file() << "fA = " << fA[j] << "  K = " << shifts[k]
                        << "  f - exact = " << f - fExact
                        << "  mu - exact = " << mu - muExact
                        << std::endl;
            TEST_ASSERT(std::abs(f - fExact) < 1.0E-12);
            TEST_ASSERT(std::abs(mu - muExact) < 1.0E-12);
            TEST_ASSERT(std::abs(sys.cField(0)[0] - fA[j]) < 1.0E-12);
         }
      }
   }

   /*
   * Homopolymer of length 2000 in a solvent (phi 0.1, chi 0.5). The
   * partition function of the polymer, exp(-N w), is below the double
   * range; the homogeneous solution follows Flory-Huggins theory.
   */
   void testLongChain()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1GaugeTestLongChain.log");

      System sys;
      readSystem(sys, "in/R1/long_chain.prm");
      runCommands(sys, "READ_W in/R1/long_chain.w\nITERATE\nFINISH\n");

      Polymer const & polymer = sys.mixture().polymer(0);
      Solvent const & solvent = sys.mixture().solvent(0);
      double n = polymer.length();
      double phiP = polymer.phi();
      double phiS = solvent.phi();
      double chi = sys.interaction().chi(0,1);
      int nx = sys.domain().nx();

      // Flory-Huggins free energy and chemical potentials in the gauge
      // of the converged (uniform) fields
      double f = phiP*(std::log(phiP) - 1.0)/n
               + phiS*(std::log(phiS) - 1.0) + chi*phiP*phiS;
      double wP = sys.domain().spatialAverage(sys.wField(0));
      double wS = sys.domain().spatialAverage(sys.wField(1));
      double muP = std::log(phiP) + n*wP;
      double muS = std::log(phiS) + wS;
      Log::file() << "f = " << sys.fHelmholtz() << "  (FH " << f << ")"
                  << "  mu = " << polymer.mu() << "  (FH " << muP << ")"
                  << std::endl;
      TEST_ASSERT(std::isfinite(sys.fHelmholtz()));
      TEST_ASSERT(std::abs(sys.fHelmholtz() - f) < 1.0E-8);
      TEST_ASSERT(std::abs(polymer.mu() - muP) < 1.0E-6);
      TEST_ASSERT(std::abs(solvent.mu() - muS) < 1.0E-8);
      TEST_ASSERT(n*wP > 750.0);  // q = exp(-N wP) below the double range
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(std::abs(sys.cField(0)[i] - phiP) < 1.0E-8);
         TEST_ASSERT(std::abs(sys.wField(0)[i] - wP) < 1.0E-6);
      }
   }

};

TEST_BEGIN(R1GaugeTest)
TEST_ADD(R1GaugeTest, testShiftInvariance)
TEST_ADD(R1GaugeTest, testLabelSwap)
TEST_ADD(R1GaugeTest, testHomogeneous)
TEST_ADD(R1GaugeTest, testLongChain)
TEST_END(R1GaugeTest)

#endif
