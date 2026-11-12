#ifndef RPC_FTS_PERTURBATION_TEST_H
#define RPC_FTS_PERTURBATION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/perturbation/Perturbation.h>
#include <rp/field/WFields.h>

#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): EinsteinCrystalPerturbation under chi ramps.
*/
class FtsPerturbationTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   template <int D>
   void initSystem(Rp::System<D,CPT>& system, std::string filename)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile(filename, in);
      system.readParam(in);
      in.close();
   }

   template <int D>
   void initSimulator(Rp::BdSimulator<D,CPT>& simulator,
                      std::string filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      simulator.readParam(in);
      in.close();
   }

   /*
   * Simulate under a chi ramp with an EinsteinCrystalPerturbation of
   * lambda = 0 (the sampled H is the Einstein crystal Hamiltonian).
   *
   * Return the Hamiltonian of the final state computed with the
   * perturbation of the run (H), and after a new setup of the
   * perturbation for the final chi (Hsetup).
   */
   template <int D>
   void rampEinsteinCrystal(std::string systemFileName,
                            std::string simulatorFileName,
                            std::string wFileName,
                            int nStep, double& H, double& Hsetup)
   {
      Rp::System<D,CPT> system;
      initSystem(system, systemFileName);
      Rp::BdSimulator<D,CPT> simulator(system);
      initSimulator(simulator, simulatorFileName);
      system.w().readRGrid(wFileName);
      simulator.simulate(nStep);
      TEST_ASSERT(simulator.iStep() == nStep);

      // Hamiltonian of the final state for the final chi
      simulator.computeWc();
      simulator.computeHamiltonian();
      H = simulator.hamiltonian();

      // The same after a new setup of the perturbation
      simulator.perturbation().setup();
      simulator.computeHamiltonian();
      Hsetup = simulator.hamiltonian();
   }

   /*
   * Binary system, default epsilon (= chi) under a chi ramp (10 to 12).
   * The reference field differs from w, so that H_EC ~ 1/epsilon is
   * not zero.
   */
   void testEinsteinCrystalChiRampEpsilon()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsEinsteinCrystalChiRampEpsilon.log");

      double H, Hsetup;
      rampEinsteinCrystal<3>("in/Fts/param_system_small",
                             "in/Fts/param_BdSimulator_ecChiRamp",
                             "in/Fts/w_small.rf", 4, H, Hsetup);
      TEST_ASSERT(std::abs(Hsetup) > 1.0);
      TEST_ASSERT(std::abs(H - Hsetup) < 1.0E-8*std::abs(Hsetup));
   }

   /*
   * ABC triblock (3 monomers), given epsilon, under a ramp of chi_BC
   * (12 to 22) that rotates the chi eigenvectors. The reference field
   * is the initial w, and the fields stay there (mobility 1e-14), so
   * H_EC must stay close to zero.
   */
   void testEinsteinCrystalChiRampBasis()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsEinsteinCrystalChiRampBasis.log");

      double H, Hsetup;
      rampEinsteinCrystal<1>("in/Fts/param_system_abc_1D",
                             "in/Fts/param_BdSimulator_ecChiRamp_abc",
                             "in/Fts/w_abc_1D.rf", 5, H, Hsetup);
      TEST_ASSERT(std::abs(Hsetup) < 1.0E-6);
      TEST_ASSERT(std::abs(H) < 1.0E-6);
   }

};

TEST_BEGIN(FtsPerturbationTest)
TEST_ADD(FtsPerturbationTest, testEinsteinCrystalChiRampEpsilon)
TEST_ADD(FtsPerturbationTest, testEinsteinCrystalChiRampBasis)
TEST_END(FtsPerturbationTest)

#endif
