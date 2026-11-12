#ifndef RPG_FTS_BD_SIMULATOR_TEST_H
#define RPG_FTS_BD_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/field/cuda/RField.h>
#include <pscf/field/cuda/RFieldComparison.h>

#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): BdSimulator: ramps and final states.
*/
class FtsBdSimulatorTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   template <int D>
   void initSystem(Rp::System<D,CUT>& system, std::string filename)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile(filename, in);
      system.readParam(in);
      in.close();
   }

   template <int D>
   void initSimulator(Rp::BdSimulator<D,CUT>& simulator,
                      std::string filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      simulator.readParam(in);
      in.close();
   }

   /*
   * A LinearRamp of the phi of a closed polymer in a mixed ensemble
   * (closed diblock, open homopolymer). The changes of phi need not
   * add up to zero.
   */
   void testPhiRampMixedEnsemble()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsPhiRampMixedEnsemble.log");

      Rp::System<1,CUT> system;
      initSystem(system, "in/Fts/param_system_mixed_1D");
      Rp::BdSimulator<1,CUT> simulator(system);
      bool success = true;
      try {
         initSimulator(simulator, "in/Fts/param_BdSimulator_phiRamp");
      } catch (Exception& e) {
         success = false;
      }
      TEST_ASSERT(success);
      TEST_ASSERT(simulator.hasRamp());

      system.w().readRGrid("in/Fts/w_mixed_1D.rf");
      simulator.simulate(4);
      TEST_ASSERT(simulator.iStep() == 4);
      double phi = system.mixture().polymer(0).phi();
      TEST_ASSERT(std::abs(phi - 0.7) < 1.0E-10);
   }

   /*
   * A LinearRamp of vMonomer (0.01 to 0.002) with an LrCompressor. The
   * intramolecular correlation must follow the parameters, or the
   * compressor fails once vMonomer has changed enough.
   */
   void testVMonomerRampLrCompressor()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsVMonomerRampLrCompressor.log");

      Rp::System<3,CUT> system;
      initSystem(system, "in/Fts/param_system_small");
      Rp::BdSimulator<3,CUT> simulator(system);
      initSimulator(simulator, "in/Fts/param_BdSimulator_vMonomerRamp");
      system.w().readRGrid("in/Fts/w_small.rf");
      simulator.simulate(60);

      // Every step converged, so the ramp reached its end
      TEST_ASSERT(simulator.iStep() == 60);
      double vMonomer = system.mixture().vMonomer();
      TEST_ASSERT(std::abs(vMonomer - 0.002) < 1.0E-10);
   }

   /*
   * Allocate and copy an array of r-grid fields.
   */
   template <int D>
   void copyFields(DArray< RField<D,CUT> > const & in,
                   DArray< RField<D,CUT> >& out)
   {
      int n = in.capacity();
      out.allocate(n);
      for (int i = 0; i < n; ++i) {
         out[i].allocate(in[i].meshDimensions());
         out[i] = in[i];
      }
   }

   /*
   * Are the c fields of the system those of its current w fields?
   */
   template <int D>
   bool hasCurrentCFields(Rp::System<D,CUT>& system)
   {
      if (!system.c().hasData()) return false;
      DArray< RField<D,CUT> > c;
      copyFields(system.c().rgrid(), c);
      system.compute();
      RFieldComparison<D,CUT> comparison;
      comparison.compare(c, system.c().rgrid());
      return (comparison.maxDiff() < 1.0E-10);
   }

   /*
   * A chi ramp (10 to 12): after the last step, the stored forces dc
   * (used by the next step) and Hamiltonian must be those of the final
   * parameters.
   */
   void testChiRampState()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsBdChiRampState.log");

      Rp::System<3,CUT> system;
      initSystem(system, "in/Fts/param_system_small");
      Rp::BdSimulator<3,CUT> simulator(system);
      initSimulator(simulator, "in/Fts/param_BdSimulator_chiRamp");
      system.w().readRGrid("in/Fts/w_small.rf");
      simulator.simulate(4);
      TEST_ASSERT(simulator.iStep() == 4);
      TEST_ASSERT(std::abs(system.interaction().chi(0,1) - 12.0) < 1.0E-10);

      // Values stored by the simulator
      TEST_ASSERT(simulator.hasDc());
      DArray< RField<3,CUT> > dc;
      copyFields(simulator.dc(), dc);
      bool hasHamiltonian = simulator.hasHamiltonian();
      double hamiltonian = 0.0;
      if (hasHamiltonian) {
         hamiltonian = simulator.hamiltonian();
      }

      // Recompute them for the final state and parameters
      simulator.computeWc();
      simulator.computeCc();
      simulator.computeDc();
      simulator.computeHamiltonian();
      RFieldComparison<3,CUT> comparison;
      comparison.compare(dc, simulator.dc());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-8);
      if (hasHamiltonian) {
         double diff = std::abs(simulator.hamiltonian() - hamiltonian);
         TEST_ASSERT(diff < 1.0E-10*std::abs(hamiltonian));
      }
   }

   /*
   * A run whose last step fails (compressor with maxItr = 1): the c
   * fields must be those of the final (restored) w fields.
   */
   void testFailedLastStep()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsFailedLastStep.log");

      Rp::System<3,CUT> system;
      initSystem(system, "in/Fts/param_system_small");
      Rp::BdSimulator<3,CUT> simulator(system);
      initSimulator(simulator, "in/Fts/param_BdSimulator_failStep");
      system.w().readRGrid("in/Fts/w_small.rf");
      simulator.simulate(3);

      // Every step failed
      TEST_ASSERT(simulator.iTotalStep() == 3);
      TEST_ASSERT(simulator.iStep() == 0);

      TEST_ASSERT(hasCurrentCFields(system));
   }

   /*
   * The generic block label "Ramp" selects a LinearRamp.
   */
   void testRampLabel()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsRampLabel.log");

      Rp::System<3,CUT> system;
      initSystem(system, "in/Fts/param_system_small");
      Rp::BdSimulator<3,CUT> simulator(system);
      bool success = true;
      try {
         initSimulator(simulator, "in/Fts/param_BdSimulator_rampLabel");
      } catch (Exception& e) {
         success = false;
      }
      TEST_ASSERT(success);
      TEST_ASSERT(simulator.hasRamp());

      system.w().readRGrid("in/Fts/w_small.rf");
      simulator.simulate(2);
      TEST_ASSERT(simulator.iStep() == 2);
      TEST_ASSERT(std::abs(system.interaction().chi(0,1) - 12.0) < 1.0E-10);
   }

};

TEST_BEGIN(FtsBdSimulatorTest)
TEST_ADD(FtsBdSimulatorTest, testPhiRampMixedEnsemble)
TEST_ADD(FtsBdSimulatorTest, testVMonomerRampLrCompressor)
TEST_ADD(FtsBdSimulatorTest, testChiRampState)
TEST_ADD(FtsBdSimulatorTest, testFailedLastStep)
TEST_ADD(FtsBdSimulatorTest, testRampLabel)
TEST_END(FtsBdSimulatorTest)

#endif
