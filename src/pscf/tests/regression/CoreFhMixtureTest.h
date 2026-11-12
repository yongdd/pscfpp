#ifndef PSCF_CORE_FH_MIXTURE_TEST_H
#define PSCF_CORE_FH_MIXTURE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/floryHuggins/FhMixture.h>
#include <pscf/floryHuggins/FhMolecule.h>
#include <pscf/floryHuggins/FhInteraction.h>
#include <util/containers/DArray.h>
#include <util/param/BracketPolicy.h>

#include <cmath>
#include <fstream>
#include <limits>

using namespace Pscf;
using namespace Util;

class CoreFhMixtureTest : public UnitTest 
{

public:

   void setUp()
   {  BracketPolicy::set(BracketPolicy::Optional); }

   void tearDown()
   {}

   /*
   * Blend of A and B homopolymers of equal size at chi N = 20, with
   * equal target chemical potentials: the A-rich solution has a B
   * fraction of about exp(-20) = 2e-9, as at the edge of a strongly
   * segregated interface. The dilute species B is the last one. The
   * last species used to take the incompressibility constraint, so
   * its fraction 1 - phi_A had no relative precision and computePhi
   * failed to converge. By symmetry, the solve from the mirror guess
   * gives the B-rich solution, in which the dilute species comes first.
   */
   void testComputePhiDiluteLast()
   {
      printMethod(TEST_FUNC);

      FhMixture mixture;
      std::ifstream in;
      openInputFile("in/Core/FhMixtureBlend", in);
      mixture.readParam(in);
      in.close();

      FhInteraction interaction;
      interaction.setNMonomer(2);
      interaction.setChi(0, 1, 20.0);

      DArray<double> mu;
      DArray<double> phi;
      mu.allocate(2);
      phi.allocate(2);
      mu[0] = -20.0;
      mu[1] = -20.0;

      // Dilute species first (B-rich solution)
      phi[0] = 1.0E-8;
      phi[1] = 1.0 - 1.0E-8;
      double xi = 0.0;
      mixture.computePhi(interaction, mu, phi, xi);
      double phiDilute = mixture.phi(0);
      mixture.computeFreeEnergy(interaction);
      double fHelmholtz = mixture.fHelmholtz();
      double pressure = mixture.pressure();

      // Dilute species last (A-rich solution)
      phi[0] = 1.0 - 1.0E-8;
      phi[1] = 1.0E-8;
      xi = 0.0;
      mixture.computePhi(interaction, mu, phi, xi);
      TEST_ASSERT(eq(mixture.mu(0), mu[0]));
      TEST_ASSERT(eq(mixture.mu(1), mu[1]));
      double phiB = mixture.phi(1);
      double phiA = mixture.phi(0);
      TEST_ASSERT(phiB > 1.0E-9);
      TEST_ASSERT(phiB < 3.0E-9);
      TEST_ASSERT(std::abs(phiA + phiB - 1.0) < 1.0E-15);
      TEST_ASSERT(std::abs(phiB/phiDilute - 1.0) < 1.0E-9);

      // Equal mu: log(phiB/phiA) + chi*(phiA - phiB) = 0
      double residual = std::log(phiB/phiA) + 20.0*(phiA - phiB);
      TEST_ASSERT(std::abs(residual) < 1.0E-9);

      mixture.computeFreeEnergy(interaction);
      TEST_ASSERT(eq(mixture.fHelmholtz(), fHelmholtz));
      TEST_ASSERT(eq(mixture.pressure(), pressure));
   }

   /*
   * A closed species with phi = 0 in the SCFT solution has mu = -inf.
   * Mixture of A and B homopolymers and an AB diblock with phi = 0, at
   * chi N = 1.5 (one phase): the results must be those of the A/B
   * blend without the diblock. The free energy was NaN (0*log(0)), and
   * computePhi gave NaN fractions and aborted.
   */
   void testComputePhiAbsentSpecies()
   {
      printMethod(TEST_FUNC);

      FhInteraction interaction;
      interaction.setNMonomer(2);
      interaction.setChi(0, 1, 1.5);

      FhMixture blend2;
      FhMixture blend3;
      readBlends(blend2, blend3);

      // Two-species reference state, phi = (0.75, 0.25)
      DArray<double> phi2;
      phi2.allocate(2);
      phi2[0] = 0.75;
      phi2[1] = 0.25;
      blend2.setComposition(phi2);
      blend2.computeMu(interaction, 0.0);
      blend2.computeFreeEnergy(interaction);
      double mu0 = blend2.mu(0);
      double mu1 = blend2.mu(1);
      double fHelmholtz = blend2.fHelmholtz();
      double pressure = blend2.pressure();

      // Same composition with the diblock at phi = 0
      DArray<double> phi3;
      phi3.allocate(3);
      phi3[0] = 0.75;
      phi3[1] = 0.25;
      phi3[2] = 0.0;
      blend3.setComposition(phi3);
      TEST_ASSERT(blend3.phi(2) == 0.0);
      blend3.computeMu(interaction, 0.0);
      blend3.computeFreeEnergy(interaction);
      TEST_ASSERT(eq(blend3.mu(0), mu0));
      TEST_ASSERT(eq(blend3.mu(1), mu1));
      TEST_ASSERT(eq(blend3.fHelmholtz(), fHelmholtz));
      TEST_ASSERT(eq(blend3.pressure(), pressure));

      // Two-species solve from phi = (0.5, 0.5)
      DArray<double> mu2;
      mu2.allocate(2);
      mu2[0] = mu0;
      mu2[1] = mu1;
      phi2[0] = 0.5;
      phi2[1] = 0.5;
      double xi = 0.0;
      blend2.computePhi(interaction, mu2, phi2, xi);
      TEST_ASSERT(eq(blend2.phi(0), 0.75));
      TEST_ASSERT(eq(blend2.phi(1), 0.25));
      blend2.computeFreeEnergy(interaction);

      // Three-species solve, diblock absent (mu = -inf, phi = 0)
      DArray<double> mu3;
      mu3.allocate(3);
      mu3[0] = mu0;
      mu3[1] = mu1;
      mu3[2] = -std::numeric_limits<double>::infinity();
      phi3[0] = 0.5;
      phi3[1] = 0.5;
      phi3[2] = 0.0;
      xi = 0.0;
      blend3.computePhi(interaction, mu3, phi3, xi);
      blend3.computeFreeEnergy(interaction);
      checkSameState(blend3, blend2);
      TEST_ASSERT(eq(blend3.fHelmholtz(), fHelmholtz));
      TEST_ASSERT(eq(blend3.pressure(), pressure));
   }

private:

   // Read the A/B blend and the A/B blend with an AB diblock
   void readBlends(FhMixture& blend2, FhMixture& blend3)
   {
      std::ifstream in;
      openInputFile("in/Core/FhMixtureBlend", in);
      blend2.readParam(in);
      in.close();
      openInputFile("in/Core/FhMixtureBlend3", in);
      blend3.readParam(in);
      in.close();
   }

   /*
   * Check that blend3, with an absent third species, is in the state
   * of the two-species blend2 (after computePhi and computeFreeEnergy).
   */
   void checkSameState(FhMixture const & blend3, FhMixture const & blend2)
   {
      TEST_ASSERT(eq(blend3.phi(0), blend2.phi(0)));
      TEST_ASSERT(eq(blend3.phi(1), blend2.phi(1)));
      TEST_ASSERT(blend3.phi(2) == 0.0);
      TEST_ASSERT(eq(blend3.mu(0), blend2.mu(0)));
      TEST_ASSERT(eq(blend3.mu(1), blend2.mu(1)));
      TEST_ASSERT(eq(blend3.c(0), blend2.c(0)));
      TEST_ASSERT(eq(blend3.c(1), blend2.c(1)));
      TEST_ASSERT(eq(blend3.fHelmholtz(), blend2.fHelmholtz()));
      TEST_ASSERT(eq(blend3.pressure(), blend2.pressure()));
   }

};

TEST_BEGIN(CoreFhMixtureTest)
TEST_ADD(CoreFhMixtureTest, testComputePhiDiluteLast)
TEST_ADD(CoreFhMixtureTest, testComputePhiAbsentSpecies)
TEST_END(CoreFhMixtureTest)

#endif
