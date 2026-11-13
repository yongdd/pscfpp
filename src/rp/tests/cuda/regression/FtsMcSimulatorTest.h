#ifndef RPG_FTS_MC_SIMULATOR_TEST_H
#define RPG_FTS_MC_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/tests/CompressedInitialState.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/fts/montecarlo/McMove.h>
#include <rp/fts/montecarlo/McMoveManager.h>
#include <rp/fts/compressor/Compressor.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/field/cuda/RField.h>
#include <pscf/field/cuda/RFieldComparison.h>
#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/Reduce.h>

#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): McSimulator: ramps, final states and BdMove.
*/
class FtsMcSimulatorTest : public LogFileUnitTest
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
   void initSimulator(Rp::McSimulator<D,CUT>& simulator,
                      std::string filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      simulator.readParam(in);
      in.close();
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
   * Null moves (RealMove, sigma 1e-7) under a chi ramp (10 to 16). The
   * stored Hamiltonian must follow the ramp, so that every move is
   * accepted.
   */
   void testChiRampState()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsMcChiRampState.log");

      Rp::System<3,CUT> system;
      initSystem(system, "in/Fts/param_system_small");
      Rp::McSimulator<3,CUT> simulator(system);
      initSimulator(simulator, "in/Fts/param_McSimulator_chiRamp");
      system.w().readRGrid("in/Fts/w_small.rf");
      simulator.simulate(10);
      TEST_ASSERT(std::abs(system.interaction().chi(0,1) - 16.0) < 1.0E-10);

      Rp::McMove<3,CUT>& move = simulator.mcMoveManager()[0];
      TEST_ASSERT(move.nAttempt() == 10);
      TEST_ASSERT(move.nAccept() == 10);

      // The stored Hamiltonian is that of the final parameters
      TEST_ASSERT(simulator.hasHamiltonian());
      double hamiltonian = simulator.hamiltonian();
      simulator.computeWc();
      simulator.computeHamiltonian();
      double diff = std::abs(simulator.hamiltonian() - hamiltonian);
      TEST_ASSERT(diff < 1.0E-10*std::abs(hamiltonian));
   }

   /*
   * A run in which every move is rejected (RealMove, sigma 1.0): the c
   * fields must be those of the final (restored) w fields.
   */
   void testRejectedLastMove()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsRejectedLastMove.log");

      Rp::System<3,CUT> system;
      initSystem(system, "in/Fts/param_system_small");
      Rp::McSimulator<3,CUT> simulator(system);
      initSimulator(simulator, "in/Fts/param_McSimulator_reject");
      system.w().readRGrid("in/Fts/w_small.rf");
      simulator.simulate(3);

      Rp::McMove<3,CUT>& move = simulator.mcMoveManager()[0];
      TEST_ASSERT(move.nAttempt() == 3);
      TEST_ASSERT(move.nAccept() == 0);

      TEST_ASSERT(hasCurrentCFields(system));
   }

   /*
   * BdMove with a compressor that fails in every BD step (maxItr = 1),
   * so also in the last one. The McSimulator is that of the System, as
   * in pscf_rpc and pscf_rpg. The moves must complete, and leave the c
   * fields of the final w fields.
   */
   void testBdMoveFailedLastStep()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsBdMoveFailedLastStep.log");

      Rp::System<3,CUT> system;
      initSystem(system, "in/Fts/param_system_small_BdMove_fail");
      TEST_ASSERT(system.hasSimulator());
      system.w().readRGrid("in/Fts/w_small.rf");
      Rp::prepareCompressedInitialState(system, system.simulator());

      bool success = true;
      try {
         system.simulator().simulate(2);
      } catch (Exception& e) {
         success = false;
      }
      TEST_ASSERT(success);
      TEST_ASSERT(hasCurrentCFields(system));
   }

   /*
   * Noise of a BdMove with a single LM step (nStep = 1).
   *
   * The step adds eta - mobility*dc to the exchange component wc(0),
   * where eta, the sum of the "old" and "new" noise, must have the
   * Euler-Maruyama variance 2*mobility*meshSize/vSystem at each grid
   * point (half of it if the old noise has the usual variance).
   */
   void testBdMoveFirstStepNoise()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsBdMoveFirstStepNoise.log");

      Rp::System<3,CUT> system;
      initSystem(system, "in/Fts/param_system_disordered");
      Rp::McSimulator<3,CUT> simulator(system);
      initSimulator(simulator, "in/Fts/param_McSimulator_BdMove_noise");
      const double mobility = 1.0E-4; // value in the parameter file

      // Compressed initial state, its wc(0) and force dc(0)
      system.w().readRGrid("in/Fts/w_dis.rf");
      simulator.analyzeChi();
      TEST_ASSERT(simulator.compressor().compress() == 0);
      simulator.computeWc();
      simulator.computeCc();
      simulator.computeDc();

      // eta = wc(0) after - wc(0) before + mobility*dc(0) before
      RField<3,CUT> eta;
      eta.allocate(system.domain().mesh().dimensions());
      VecOp::eqV(eta, simulator.wc(0));
      VecOp::mulEqS(eta, -1.0);
      VecOp::addEqVc(eta, simulator.dc(0), mobility);
      simulator.simulate(1);
      simulator.computeWc();
      VecOp::addEqV(eta, simulator.wc(0));

      // Variance of eta over the 16^3 grid points (relative standard
      // error of the estimate about 0.022)
      const int meshSize = system.domain().mesh().size();
      const double vSystem = system.domain().unitCell().volume();
      double mean = Reduce::sum(eta)/double(meshSize);
      double variance = Reduce::sumSq(eta)/double(meshSize) - mean*mean;
      double expected = 2.0*mobility*double(meshSize)/vSystem;
      Log::file() << "Variance of eta / expected value = "
                  << variance/expected << std::endl;
      TEST_ASSERT(std::abs(variance/expected - 1.0) < 0.1);
   }

};

TEST_BEGIN(FtsMcSimulatorTest)
TEST_ADD(FtsMcSimulatorTest, testChiRampState)
TEST_ADD(FtsMcSimulatorTest, testRejectedLastMove)
TEST_ADD(FtsMcSimulatorTest, testBdMoveFailedLastStep)
TEST_ADD(FtsMcSimulatorTest, testBdMoveFirstStepNoise)
TEST_END(FtsMcSimulatorTest)

#endif
