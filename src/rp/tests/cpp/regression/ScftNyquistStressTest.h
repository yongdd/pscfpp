#ifndef RPC_SCFT_NYQUIST_STRESS_TEST_H
#define RPC_SCFT_NYQUIST_STRESS_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>

#include <util/containers/FSArray.h>
#include <util/format/Dbl.h>
#include <util/math/Constants.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Stress of an orthogonal triclinic cell, on an even mesh, with a
* rough field (strong Nyquist components), for the bead model, whose
* stress is the exact derivative of -(phi/N) ln Q at fixed w.
*
* The minimum images of wavevectors with a Nyquist index switch side
* where a cell angle passes through 90 degrees, so that ln Q had a kink
* there and the stress at 90 degrees (taken from one side) disagreed
* with the finite difference. With kSq averaged over the signs of the
* Nyquist indices, ln Q is smooth.
*/
class ScftNyquistStressTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   void setupSystem(Rp::System<3,CPT>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

   /*
   * Stress n for cell parameters p, and the central finite difference
   * of -(phi/N) ln Q, at fixed w fields.
   */
   void stressAndDifference(Rp::System<3,CPT>& system,
                            FSArray<double, 6> const & p, int n,
                            double& stress, double& stressFD)
   {
      double h = 1.0E-6;
      system.setUnitCell(p);
      system.compute(true);
      stress = system.mixture().stress(n);

      FSArray<double, 6> pp = p;
      pp[n] += h;
      system.setUnitCell(pp);
      system.compute();
      double qp = system.mixture().polymer(0).q();

      FSArray<double, 6> pm = p;
      pm[n] -= h;
      system.setUnitCell(pm);
      system.compute();
      double qm = system.mixture().polymer(0).q();

      double phi = system.mixture().polymer(0).phi();
      double nBead = (double) system.mixture().polymer(0).nBead();
      stressFD = -(phi/nBead)*(std::log(qp) - std::log(qm))/(2.0*h);
   }

   /*
   * All six parameters of the orthogonal cell, and theta and gamma on
   * both sides of the orthogonal values: stress equals the finite
   * difference.
   */
   void testRightAngles()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftNyquistStress.log");

      Rp::System<3,CPT> system;
      setupSystem(system, "in/Scft/nyq_tric_param");
      system.w().readRGrid("in/Scft/nyq_tric_w.rf");

      // Triclinic parameters (a, b, c, phi, theta, gamma): the cell is
      // orthogonal for theta = 0 and gamma = 90 degrees
      double halfPi = 0.5*Constants::Pi;
      FSArray<double, 6> p;
      p.append(5.0);
      p.append(4.0);
      p.append(6.0);
      p.append(0.4);
      p.append(0.0);
      p.append(halfPi);

      double stress, stressFD, scale;
      double maxError = 0.0;
      for (int n = 0; n < 6; ++n) {
         stressAndDifference(system, p, n, stress, stressFD);
         scale = std::max(1.0E-2, std::abs(stressFD));
         Log::file() << "orthogonal  param " << n
                     << "  stress " << Dbl(stress, 20, 12)
                     << "  fin diff " << Dbl(stressFD, 20, 12)
                     << std::endl;
         maxError = std::max(maxError, std::abs(stress - stressFD)/scale);
      }

      // Tilt theta = +-1e-3 and gamma = 90 deg +- 1e-3
      double delta[2] = {-1.0E-3, 1.0E-3};
      for (int n = 4; n < 6; ++n) {
         for (int s = 0; s < 2; ++s) {
            FSArray<double, 6> q = p;
            q[n] += delta[s];
            stressAndDifference(system, q, n, stress, stressFD);
            scale = std::max(1.0E-2, std::abs(stressFD));
            Log::file() << "shift " << (s ? "+" : "-")
                        << " 1e-3  param " << n
                        << "  stress " << Dbl(stress, 20, 12)
                        << "  fin diff " << Dbl(stressFD, 20, 12)
                        << std::endl;
            maxError = std::max(maxError,
                                std::abs(stress - stressFD)/scale);
         }
      }
      Log::file() << "max relative error " << maxError << std::endl;
      TEST_ASSERT(maxError < 1.0E-6);
   }

};

TEST_BEGIN(ScftNyquistStressTest)
TEST_ADD(ScftNyquistStressTest, testRightAngles)
TEST_END(ScftNyquistStressTest)

#endif
