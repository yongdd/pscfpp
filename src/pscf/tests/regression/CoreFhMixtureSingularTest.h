#ifndef PSCF_CORE_FH_MIXTURE_SINGULAR_TEST_H
#define PSCF_CORE_FH_MIXTURE_SINGULAR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/floryHuggins/FhMixture.h>
#include <pscf/floryHuggins/FhInteraction.h>
#include <util/containers/DArray.h>
#include <util/param/BracketPolicy.h>

#include <cmath>
#include <fstream>

using namespace Pscf;
using namespace Util;

/*
* Regression tests (pscfpp fixes): FhMixture::computePhi with a
* singular Jacobian.
*/
class CoreFhMixtureSingularTest : public UnitTest
{

public:

   void setUp()
   {  BracketPolicy::set(BracketPolicy::Optional); }

   void tearDown()
   {}

   /*
   * Blend of A and B homopolymers of size 1 at chi = 2, the critical
   * point of the symmetric blend. At the guess phi = (0.5, 0.5), the
   * Jacobian of the Newton iteration of computePhi is exactly singular
   * (zero pivot). LuSolver::solve then called the GSL error handler,
   * which aborted the program. computePhi must throw an Exception.
   */
   void testSpinodalGuess()
   {
      printMethod(TEST_FUNC);

      FhMixture mixture;
      std::ifstream in;
      openInputFile("in/Core/FhMixtureBlend", in);
      mixture.readParam(in);
      in.close();

      FhInteraction interaction;
      interaction.setNMonomer(2);
      interaction.setChi(0, 1, 2.0);

      DArray<double> mu;
      DArray<double> phi;
      mu.allocate(2);
      phi.allocate(2);
      mu[0] = 0.1;
      mu[1] = 0.0;

      // Guess at the critical point: Exception
      phi[0] = 0.5;
      phi[1] = 0.5;
      double xi = 0.0;
      bool thrown = false;
      try {
         mixture.computePhi(interaction, mu, phi, xi);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(thrown);

      // Guess away from the critical point: converges to the solution,
      // ln(phiA/phiB) + chi*(phiB - phiA) = muA - muB
      phi[0] = 0.7;
      phi[1] = 0.3;
      xi = 0.0;
      mixture.computePhi(interaction, mu, phi, xi);
      TEST_ASSERT(eq(mixture.mu(0), mu[0]));
      TEST_ASSERT(eq(mixture.mu(1), mu[1]));
      double phiA = mixture.phi(0);
      double phiB = mixture.phi(1);
      TEST_ASSERT(phiA > 0.5);
      TEST_ASSERT(std::abs(phiA + phiB - 1.0) < 1.0E-14);
      double residual = std::log(phiA/phiB) + 2.0*(phiB - phiA) - 0.1;
      TEST_ASSERT(std::abs(residual) < 1.0E-9);
   }

};

TEST_BEGIN(CoreFhMixtureSingularTest)
TEST_ADD(CoreFhMixtureSingularTest, testSpinodalGuess)
TEST_END(CoreFhMixtureSingularTest)

#endif
