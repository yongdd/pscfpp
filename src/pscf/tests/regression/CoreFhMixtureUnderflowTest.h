#ifndef PSCF_CORE_FH_MIXTURE_UNDERFLOW_TEST_H
#define PSCF_CORE_FH_MIXTURE_UNDERFLOW_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/floryHuggins/FhMixture.h>
#include <pscf/floryHuggins/FhMolecule.h>
#include <pscf/floryHuggins/FhInteraction.h>
#include <util/containers/DArray.h>
#include <util/param/BracketPolicy.h>
#include <util/global.h>

#include <cmath>
#include <fstream>

using namespace Pscf;
using namespace Util;

class CoreFhMixtureUnderflowTest : public UnitTest
{

public:

   void setUp()
   {  BracketPolicy::set(BracketPolicy::Optional); }

   void tearDown()
   {}

   /*
   * Blend of A and B homopolymers of size 1 at chi = 2, with target
   * chemical potentials -1 for the majority species and muDilute for
   * the dilute one. The dilute fraction is e = exp(lnPhi) with
   * lnPhi = muDilute - chi + 1 (to relative order e), and the other
   * fraction is 1. Check computePhi from the guesses 1e-310 (subnormal,
   * as COMPARE_HOMOGENEOUS takes it from an SCFT solution), 1e-300 and
   * 1e-30 of the dilute fraction, with the dilute species first and
   * last.
   */
   void checkDilute(double muDilute)
   {
      const double chi = 2.0;
      const double muMajor = -1.0;
      const double lnPhi = muDilute - chi - muMajor;

      FhInteraction interaction;
      interaction.setNMonomer(2);
      interaction.setChi(0, 1, chi);

      DArray<double> mu;
      DArray<double> phi;
      mu.allocate(2);
      phi.allocate(2);

      double guesses[3] = {1.0E-310, 1.0E-300, 1.0E-30};
      for (int ig = 0; ig < 3; ++ig) {
         for (int dilute = 0; dilute < 2; ++dilute) {
            FhMixture mixture;
            std::ifstream in;
            openInputFile("in/Core/FhMixtureBlend", in);
            mixture.readParam(in);
            in.close();

            mu[dilute] = muDilute;
            mu[1 - dilute] = muMajor;
            phi[dilute] = guesses[ig];
            phi[1 - dilute] = 1.0 - guesses[ig];
            double xi = 0.0;
            bool hasThrown = false;
            try {
               mixture.computePhi(interaction, mu, phi, xi);
            } catch (Exception& e) {
               hasThrown = true;
            }
            TEST_ASSERT(!hasThrown);

            // Chemical potentials reach their targets
            TEST_ASSERT(std::abs(mixture.mu(dilute) - muDilute) < 1.0E-9);
            TEST_ASSERT(std::abs(mixture.mu(1 - dilute) - muMajor)
                        < 1.0E-9);

            // Volume fractions (the dilute one may be subnormal or 0)
            double e = mixture.phi(dilute);
            double eRef = std::exp(lnPhi);
            TEST_ASSERT(e >= 0.0);
            TEST_ASSERT(std::abs(e - eRef) <= 1.0E-9*eRef + 1.0E-321);
            TEST_ASSERT(std::abs(mixture.phi(1 - dilute) - 1.0) < 1.0E-15);

            // Free energy and pressure of pure A: f = -1, p = 0
            mixture.computeFreeEnergy(interaction);
            TEST_ASSERT(std::abs(mixture.fHelmholtz() + 1.0) < 1.0E-12);
            TEST_ASSERT(std::abs(mixture.pressure()) < 1.0E-12);
         }
      }
   }

   /*
   * Dilute fraction of a normal double (exp(-701) = 1.5e-305). The
   * subnormal guess gave 1/phi = inf in the jacobian of computePhi.
   */
   void testComputePhiDiluteNormal()
   {
      printMethod(TEST_FUNC);
      checkDilute(-700.0);
   }

   /*
   * Subnormal dilute fraction (exp(-718) = 1.6e-312). computePhi
   * formed 1/phi = inf in its jacobian and threw 'Failed to converge'.
   */
   void testComputePhiDiluteSubnormal()
   {
      printMethod(TEST_FUNC);
      checkDilute(-717.0);
   }

   /*
   * Dilute fraction below the smallest subnormal (exp(-801) = 0):
   * phi underflows to 0 while ln(phi) and mu stay finite.
   */
   void testComputePhiDiluteZero()
   {
      printMethod(TEST_FUNC);
      checkDilute(-800.0);
   }

};

TEST_BEGIN(CoreFhMixtureUnderflowTest)
TEST_ADD(CoreFhMixtureUnderflowTest, testComputePhiDiluteNormal)
TEST_ADD(CoreFhMixtureUnderflowTest, testComputePhiDiluteSubnormal)
TEST_ADD(CoreFhMixtureUnderflowTest, testComputePhiDiluteZero)
TEST_END(CoreFhMixtureUnderflowTest)

#endif
