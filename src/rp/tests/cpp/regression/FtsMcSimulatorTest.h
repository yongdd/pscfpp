#ifndef RPC_FTS_MC_SIMULATOR_TEST_H
#define RPC_FTS_MC_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/fts/montecarlo/McMove.h>
#include <rp/fts/montecarlo/McMoveManager.h>
#include <rp/fts/compressor/Compressor.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldComparison.h>

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
   void initSimulator(Rp::McSimulator<D,CPT>& simulator,
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
   void copyFields(DArray< RField<D,CPT> > const & in,
                   DArray< RField<D,CPT> >& out)
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
   bool hasCurrentCFields(Rp::System<D,CPT>& system)
   {
      if (!system.c().hasData()) return false;
      DArray< RField<D,CPT> > c;
      copyFields(system.c().rgrid(), c);
      system.compute();
      RFieldComparison<D,CPT> comparison;
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

      Rp::System<3,CPT> system;
      initSystem(system, "in/Fts/param_system_small");
      Rp::McSimulator<3,CPT> simulator(system);
      initSimulator(simulator, "in/Fts/param_McSimulator_chiRamp");
      system.w().readRGrid("in/Fts/w_small.rf");
      simulator.simulate(10);
      TEST_ASSERT(std::abs(system.interaction().chi(0,1) - 16.0) < 1.0E-10);

      Rp::McMove<3,CPT>& move = simulator.mcMoveManager()[0];
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

      Rp::System<3,CPT> system;
      initSystem(system, "in/Fts/param_system_small");
      Rp::McSimulator<3,CPT> simulator(system);
      initSimulator(simulator, "in/Fts/param_McSimulator_reject");
      system.w().readRGrid("in/Fts/w_small.rf");
      simulator.simulate(3);

      Rp::McMove<3,CPT>& move = simulator.mcMoveManager()[0];
      TEST_ASSERT(move.nAttempt() == 3);
      TEST_ASSERT(move.nAccept() == 0);

      TEST_ASSERT(hasCurrentCFields(system));
   }

};

TEST_BEGIN(FtsMcSimulatorTest)
TEST_ADD(FtsMcSimulatorTest, testChiRampState)
TEST_ADD(FtsMcSimulatorTest, testRejectedLastMove)
TEST_END(FtsMcSimulatorTest)

#endif
