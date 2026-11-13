#ifndef RPC_FTS_LINEAR_RAMP_TEST_H
#define RPC_FTS_LINEAR_RAMP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/tests/CompressedInitialState.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/field/WFields.h>
#include <pscf/interaction/Interaction.h>

#include <util/misc/Exception.h>
#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): LinearRamp parameters and end values.
*/
class FtsLinearRampTest : public LogFileUnitTest
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

   /*
   * Read a simulator parameter file, in which the string @PARAMETER@
   * is replaced by the given ramp parameter.
   */
   template <class SimulatorT>
   void initSimulator(SimulatorT& simulator, std::string filename,
                      std::string parameter = "")
   {
      std::ifstream file;
      openInputFile(filename, file);
      std::stringstream buffer;
      buffer << file.rdbuf();
      file.close();
      std::string text = buffer.str();
      std::string key = "@PARAMETER@";
      std::size_t pos = text.find(key);
      if (pos != std::string::npos) {
         text.replace(pos, key.size(), parameter);
      }

      std::istringstream in(text);
      simulator.readParam(in);
   }

   /*
   * Run a BD simulation of a diblock melt (nPolymer = 1, nBlock = 2,
   * nSolvent = 0, nMonomer = 2, cubic cell) with one ramp parameter.
   * Return the message of the Exception thrown, or "" if none.
   */
   std::string rampMessage(std::string const & parameter)
   {
      Rp::System<3,CPT> system;
      initSystem(system, "in/Fts/param_system_small");
      Rp::BdSimulator<3,CPT> simulator(system);
      initSimulator(simulator, "in/Fts/linRamp_ids", parameter);
      system.w().readRGrid("in/Fts/w_small.rf");
      std::string message;
      try {
         simulator.simulate(2);
      } catch (Exception& e) {
         message = e.message();
      }
      Log::file() << "Parameter " << parameter << ": message '"
                  << message << "'" << std::endl;
      return message;
   }

   /*
   * Does message contain the string s?
   */
   bool contains(std::string const & message, std::string const & s)
   {  return (message.find(s) != std::string::npos); }

   /*
   * Out-of-range indices of ramp parameters stop SIMULATE with a
   * message that names the index. They were not checked, and stopped
   * on unrelated assertions, or (block, last) read and wrote past the
   * array of blocks.
   */
   void testIdsOutOfRange()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsLinearRampIds.log");

      TEST_ASSERT(rampMessage("chi 0 1 0.5") == "");
      TEST_ASSERT(contains(rampMessage("phi_polymer 1 0.0"),
                           "polymer index 1 out of range"));
      TEST_ASSERT(contains(rampMessage("phi_solvent 0 0.0"),
                           "solvent index 0 out of range"));
      TEST_ASSERT(contains(rampMessage("kuhn 2 0.1"),
                           "monomer index 2 out of range"));
      TEST_ASSERT(contains(rampMessage("chi 0 2 0.5"),
                           "monomer index 2 out of range"));
      TEST_ASSERT(contains(rampMessage("cell_param 1 0.1"),
                           "lattice parameter index 1 out of range"));
      TEST_ASSERT(contains(rampMessage("block 0 2 0.1"),
                           "block index 2 out of range"));
   }

   /*
   * A closed blend (diblock + two homopolymers) with phi changes +0.05
   * and -0.04999999 that cancel only to 1e-8. The volume fractions
   * must be rescaled to add up to 1: with a sum of 1 + s*1e-8, the
   * compressor error stalled at s*1e-8, above epsilon = 1e-9, and 8
   * of the 10 steps failed.
   */
   void testClosedPhiSum()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsLinearRampSum.log");

      Rp::System<1,CPT> system;
      initSystem(system, "in/Fts/linRamp_mix3_system");
      Rp::BdSimulator<1,CPT> simulator(system);
      initSimulator(simulator, "in/Fts/linRamp_mix3_bd");
      system.w().readRGrid("in/Fts/linRamp_mix3_w.rf");
      simulator.simulate(10);

      TEST_ASSERT(simulator.iStep() == 10);
      double sum = 0.0;
      for (int i = 0; i < 3; ++i) {
         sum += system.mixture().polymer(i).phi();
      }
      if (verbose() > 0) {
         std::cout << "\nsum of phi - 1 = " << sum - 1.0;
      }
      TEST_ASSERT(std::abs(sum - 1.0) < 1.0E-14);
   }

   /*
   * BD steps that fail (compressor with maxItr = 1, so that every step
   * fails): the chi ramp (10 to 12) must still end at 12 after the
   * last step. It followed the accepted steps only, and stayed at 10.
   */
   void testFailedStepsBd()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsLinearRampFailBd.log");

      Rp::System<3,CPT> system;
      initSystem(system, "in/Fts/param_system_small");
      Rp::BdSimulator<3,CPT> simulator(system);
      initSimulator(simulator, "in/Fts/linRamp_failBd");
      system.w().readRGrid("in/Fts/w_small.rf");
      Rp::prepareCompressedInitialState(system, simulator);
      simulator.simulate(3);

      TEST_ASSERT(simulator.iTotalStep() == 3);
      TEST_ASSERT(simulator.iStep() == 0);
      double chi = system.interaction().chi(0,1);
      TEST_ASSERT(std::abs(chi - 12.0) < 1.0E-12);
   }

   /*
   * As testFailedStepsBd, for MC moves.
   */
   void testFailedStepsMc()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsLinearRampFailMc.log");

      Rp::System<3,CPT> system;
      initSystem(system, "in/Fts/param_system_small");
      Rp::McSimulator<3,CPT> simulator(system);
      initSimulator(simulator, "in/Fts/linRamp_failMc");
      system.w().readRGrid("in/Fts/w_small.rf");
      Rp::prepareCompressedInitialState(system, simulator);
      simulator.simulate(3);

      TEST_ASSERT(simulator.iTotalStep() == 3);
      TEST_ASSERT(simulator.iStep() == 0);
      double chi = system.interaction().chi(0,1);
      TEST_ASSERT(std::abs(chi - 12.0) < 1.0E-12);
   }

};

TEST_BEGIN(FtsLinearRampTest)
TEST_ADD(FtsLinearRampTest, testIdsOutOfRange)
TEST_ADD(FtsLinearRampTest, testClosedPhiSum)
TEST_ADD(FtsLinearRampTest, testFailedStepsBd)
TEST_ADD(FtsLinearRampTest, testFailedStepsMc)
TEST_END(FtsLinearRampTest)

#endif
