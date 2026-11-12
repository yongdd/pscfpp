#ifndef R1_HOMOGENEOUS_UNDERFLOW_TEST_H
#define R1_HOMOGENEOUS_UNDERFLOW_TEST_H

#include "R1HomogeneousTest.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <r1/solvers/Solvent.h>
#include <pscf/interaction/Interaction.h>

#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

class R1HomogeneousUnderflowTest : public R1HomogeneousTest
{

public:

   /*
   * COMPARE_HOMOGENEOUS in modes 1 and 2 for a closed AB diblock
   * (phi = 1, length 1, chi = 20) with an open solvent of size 0.1 and
   * chemical potential mu. In the homogeneous mixture the diblock has
   * phi = 1 to relative order phiS, w_A = w_B = chi/2 and xi =
   * muP - chi/2, where muP is the diblock chemical potential, so that
   * ln(phiS) = mu - 0.1*muP, f = chi/4 - 1 and p = -f + muP.
   */
   void checkOpenSolvent(double mu)
   {
      System sys;
      readSystem(sys, "in/R1/planar_open_solvent.prm");
      sys.mixture().solvent(0).setMu(mu);
      runCommands(sys, "READ_W in/R1/planar.w\nFINISH\n");
      TEST_ASSERT(sys.iterate() == 0);

      double chi = sys.interaction().chi(0,1);
      double muP = sys.mixture().polymer(0).mu();
      double phiS = std::exp(mu - 0.1*muP);
      double f = 0.25*chi - 1.0;
      double p = -f + muP;
      Log::file() << "Expected: f (homo) = " << f
                  << ", p (homo) = " << p
                  << ", solvent phi = " << phiS << std::endl;

      std::string text;
      double value;
      for (int mode = 1; mode <= 2; ++mode) {
         text = compareHomogeneous(sys, mode);
         TEST_ASSERT(!text.empty());
         value = labeledValue(text, "f (homo)");
         TEST_ASSERT(std::abs(value - f) < 1.0E-10);
         value = labeledValue(text, "p (homo)");
         TEST_ASSERT(std::abs(value - p) < 1.0E-9);
         TEST_ASSERT(std::abs(speciesValue(text, 0, 1) - muP) < 1.0E-9);
         TEST_ASSERT(std::abs(speciesValue(text, 1, 1) - mu) < 1.0E-9);
         TEST_ASSERT(std::abs(speciesValue(text, 0, 2) - 1.0) < 1.0E-10);
         value = speciesValue(text, 1, 2);
         TEST_ASSERT(value >= 0.0);
         TEST_ASSERT(std::abs(value - phiS) <= 1.0E-9*phiS + 1.0E-321);
      }
   }

   /*
   * mu = -720: the solvent fraction (3.1e-313) is subnormal. The
   * Newton iteration of FhMixture::computePhi formed 1/phi = inf, and
   * COMPARE_HOMOGENEOUS 1 and 2 threw 'Failed to converge'. The case
   * mu = -700, with a solvent fraction (1.5e-304) of a normal double,
   * worked before.
   */
   void testOpenSolventSubnormal()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1HomogeneousUnderflowTestSubnormal.log");
      checkOpenSolvent(-700.0);
      checkOpenSolvent(-720.0);
   }

   /*
   * mu = -800: the solvent fraction underflows to 0.
   */
   void testOpenSolventZero()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1HomogeneousUnderflowTestZero.log");
      checkOpenSolvent(-800.0);
   }

};

TEST_BEGIN(R1HomogeneousUnderflowTest)
TEST_ADD(R1HomogeneousUnderflowTest, testOpenSolventSubnormal)
TEST_ADD(R1HomogeneousUnderflowTest, testOpenSolventZero)
TEST_END(R1HomogeneousUnderflowTest)

#endif
