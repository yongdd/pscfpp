#ifndef RPC_SCFT_TIED_IMAGE_STRESS_TEST_H
#define RPC_SCFT_TIED_IMAGE_STRESS_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/scft/ScftThermo.h>
#include <rp/solvers/Mixture.h>
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
* Stress of a hexagonal state in a rhombic cell at exactly 120 degrees,
* 13x13 mesh (no Nyquist index).
*
* At 120 degrees, minimum images of some wavevectors lie on the edges of
* the hexagonal Voronoi cell of the lattice of image shifts and are tied
* with images whose |k|^2 has a different derivative. dKSq was taken from
* one of them, so that the angle stress of the converged hexagonal state
* (diblock, f = 0.3, chi = 20) was 6.3e-6 instead of the isotropic value
* (zero deviatoric stress), and a flexible cell started at 120 degrees
* moved by 1e-5 rad (with or without the group p_6_m_m). dKSq is now
* averaged over tied images: the stress is isotropic, and the cell stays
* at 120 degrees. (Off 120 degrees the minimum-image |k|^2 still has a
* kink there, which is not changed.)
*/
class ScftTiedImageStressTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   void setupSystem(Rp::System<2,CPT>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

   /*
   * Rigid cell: the stresses of a and gamma correspond to an isotropic
   * stress tensor (area a^2 sin(gamma)).
   */
   void testRigidStress()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftTiedImageRigid.log");

      Rp::System<2,CPT> system;
      setupSystem(system, "in/Scft/hexrh_rhombic13_param");
      system.w().readRGrid("in/Scft/hexrh_rhombic13_w.rf");
      system.compute(true);
      double a = system.domain().unitCell().parameter(0);
      double gamma = system.domain().unitCell().parameter(1);
      double sigmaA = system.mixture().stress(0);
      double sigmaGamma = system.mixture().stress(1);
      double deviatoric = sigmaGamma
                        - 0.5*a*std::cos(gamma)/std::sin(gamma)*sigmaA;
      Log::file() << "stress a " << Dbl(sigmaA, 20, 12)
                  << "  gamma " << Dbl(sigmaGamma, 20, 12)
                  << "  deviatoric " << Dbl(deviatoric, 20, 12)
                  << std::endl;
      TEST_ASSERT(std::abs(deviatoric) < 1.0E-9);
   }

   /*
   * Flexible cells started at 120 degrees, without a group and with
   * p_6_m_m: the cell stays at 120 degrees.
   */
   void testFlexible()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftTiedImageFlexible.log");
      double third = 2.0*Constants::Pi/3.0;
      {
         Rp::System<2,CPT> system;
         setupSystem(system, "in/Scft/hexrh_rhombic13_param");
         system.w().readRGrid("in/Scft/hexrh_rhombic13_w.rf");
         TEST_ASSERT(system.iterate() == 0);
         double gamma = system.domain().unitCell().parameter(1);
         Log::file() << "no group: gamma - 2pi/3 "
                     << Dbl(gamma - third, 20, 12) << std::endl;
         TEST_ASSERT(std::abs(gamma - third) < 1.0E-9);
      }
      {
         Rp::System<2,CPT> system;
         setupSystem(system, "in/Scft/hexrh_p6mm13_param");
         system.w().readBasis("in/Scft/hexrh_p6mm13_w.bf");
         TEST_ASSERT(system.iterate() == 0);
         double gamma = system.domain().unitCell().parameter(1);
         Log::file() << "p_6_m_m: gamma - 2pi/3 "
                     << Dbl(gamma - third, 20, 12) << std::endl;
         TEST_ASSERT(std::abs(gamma - third) < 1.0E-9);
      }
   }

};

TEST_BEGIN(ScftTiedImageStressTest)
TEST_ADD(ScftTiedImageStressTest, testRigidStress)
TEST_ADD(ScftTiedImageStressTest, testFlexible)
TEST_END(ScftTiedImageStressTest)

#endif
