#ifndef PSCF_FH_MIXTURE_DEPENDENT_SPECIES_TEST_H
#define PSCF_FH_MIXTURE_DEPENDENT_SPECIES_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/floryHuggins/FhMixture.h>
#include <pscf/floryHuggins/FhInteraction.h>
#include <util/containers/DArray.h>
#include <cmath>
#include <limits>

using namespace Pscf;
using namespace Util;

class FhMixtureDependentSpeciesTest : public UnitTest
{
public:
   void setUp() {}
   void tearDown() {}

   void initialize(FhMixture& mix, int n, double const* sizes)
   {
      mix.setNMonomer(n);
      mix.setNMolecule(n);
      for (int i = 0; i < n; ++i) {
         mix.molecule(i).setNClump(1);
         mix.molecule(i).clump(0).setMonomerId(i);
         mix.molecule(i).clump(0).setSize(sizes[i]);
         mix.molecule(i).computeSize();
      }
   }

   // Ideal mixtures have a unique solution; changing the initial majority
   // must not affect it. Includes fractions below double's underflow limit.
   void testIdealBinary()
   {
      printMethod(TEST_FUNC);
      double gaps[] = {2, 10, 20, 30, 46, 700, 718, 800};
      for (double size : {1.0, 4.0}) {
         FhMixture mix;
         double sizes[] = {size, size};
         initialize(mix, 2, sizes);
         FhInteraction interaction;
         interaction.setNMonomer(2);
         interaction.setChi(0, 1, 0.0);
         DArray<double> mu, guess;
         mu.allocate(2); guess.allocate(2);
         for (double gap : gaps) {
            long double e = std::exp(-static_cast<long double>(gap));
            double reference = double(e/(1+e));
            for (int dilute = 0; dilute < 2; ++dilute) {
               mu[dilute] = -gap;
               mu[1-dilute] = 0.0;
               for (double firstGuess : {0.4, 0.6}) {
                  guess[dilute] = firstGuess;
                  guess[1-dilute] = 1.0-firstGuess;
                  double xi = 0.0;
                  mix.computePhi(interaction, mu, guess, xi);
                  for (int m = 0; m < 2; ++m) {
                     TEST_ASSERT(std::abs(mix.mu(m)-mu[m]) < 1e-10);
                  }
                  TEST_ASSERT(std::abs(mix.phi(0)+mix.phi(1)-1) < 1e-14);
                  if (reference == 0.0) {
                     TEST_ASSERT(mix.phi(dilute) == 0.0);
                  } else {
                     TEST_ASSERT(std::abs(mix.phi(dilute)/reference-1) < 1e-9);
                     double residual = std::log(mix.phi(dilute))
                                     - std::log(mix.phi(1-dilute)) + gap;
                     TEST_ASSERT(std::abs(residual) < 1e-9);
                  }
                  TEST_ASSERT(std::abs(xi-double(std::log1p(e))/size) < 1e-10);
                  mix.computeFreeEnergy(interaction);
                  TEST_ASSERT(std::isfinite(mix.fHelmholtz()));
                  TEST_ASSERT(std::isfinite(mix.pressure()));
               }
            }
         }
      }
   }

   // Independent mu_i = ln(phi_i) + N_i*xi for unequal molecular sizes.
   // Reuse the same object with each species in turn becoming the majority.
   void testUnequalSizesAndAbsentSpecies()
   {
      printMethod(TEST_FUNC);
      double sizes[] = {1, 4, 2};
      FhMixture mix;
      initialize(mix, 3, sizes);
      FhInteraction interaction;
      interaction.setNMonomer(3);
      for (int i = 0; i < 3; ++i) for (int j = 0; j < i; ++j)
         interaction.setChi(i, j, 0.0);
      DArray<double> mu, guess;
      mu.allocate(3); guess.allocate(3);
      for (bool absent : {false, true}) {
         for (int majority = 0; majority < 3; ++majority) {
            int dilute = (majority+1)%3;
            int third = (majority+2)%3;
            double ref[3];
            ref[dilute] = std::exp(-30.0);
            ref[third] = absent ? 0.0 : 0.2;
            ref[majority] = 1.0-ref[dilute]-ref[third];
            for (int i = 0; i < 3; ++i) {
               mu[i] = ref[i] > 0 ? std::log(ref[i])+sizes[i]*0.25
                                 : -std::numeric_limits<double>::infinity();
            }
            guess[dilute] = 0.6;
            guess[third] = absent ? 0.0 : 0.1;
            guess[majority] = 0.4-guess[third];
            double xi = 0.0;
            mix.computePhi(interaction, mu, guess, xi);
            double total = 0.0;
            for (int i = 0; i < 3; ++i) {
               total += mix.phi(i);
               if (ref[i] == 0) {
                  TEST_ASSERT(mix.phi(i) == 0);
               } else {
                  TEST_ASSERT(std::abs(mix.phi(i)/ref[i]-1) < 1e-9);
                  TEST_ASSERT(std::abs(std::log(mix.phi(i))+sizes[i]*xi-mu[i]) < 1e-10);
               }
            }
            TEST_ASSERT(std::abs(total-1) < 1e-14);
            TEST_ASSERT(std::abs(xi-0.25) < 1e-10);
         }
      }
   }
};

TEST_BEGIN(FhMixtureDependentSpeciesTest)
TEST_ADD(FhMixtureDependentSpeciesTest, testIdealBinary)
TEST_ADD(FhMixtureDependentSpeciesTest, testUnequalSizesAndAbsentSpecies)
TEST_END(FhMixtureDependentSpeciesTest)
#endif
