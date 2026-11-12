#ifndef RPC_FTS_BD_MOVE_FAIL_TEST_H
#define RPC_FTS_BD_MOVE_FAIL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/fts/montecarlo/McMove.h>
#include <rp/fts/montecarlo/McMoveManager.h>
#include <rp/fts/compressor/Compressor.h>
#include <rp/fts/compressor/LrCompressor.h>
#include <rp/fts/compressor/CompressorFactory.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/field/cpp/RField.h>
#include <pscf/backend/cpp/VecOp.h>
#include <pscf/backend/cpp/Reduce.h>

#include <util/param/Factory.h>
#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* LrCompressor that reports failure for chosen calls of compress().
*
* After setSchedule(fail), call i of compress() (counted from 0) runs
* the LrCompressor and then returns 1 if fail[i] is true. Later calls
* return the result of the LrCompressor.
*/
class FtsScheduledCompressor : public Rp::LrCompressor<3,CPT>
{

public:

   FtsScheduledCompressor(Rp::System<3,CPT>& system)
    : Rp::LrCompressor<3,CPT>(system),
      nCall_(0)
   {  ParamComposite::setClassName("FtsScheduledCompressor"); }

   void setSchedule(std::vector<bool> const & fail)
   {
      schedule_ = fail;
      nCall_ = 0;
   }

   int compress()
   {
      int error = Rp::LrCompressor<3,CPT>::compress();
      int i = nCall_;
      ++nCall_;
      if (i < (int)schedule_.size() && schedule_[i]) {
         return 1;
      }
      return error;
   }

private:

   std::vector<bool> schedule_;
   int nCall_;

};

/*
* Factory for FtsScheduledCompressor, a subfactory of the
* CompressorFactory of a simulator.
*/
class FtsScheduledCompressorFactory
 : public Factory< Rp::Compressor<3,CPT> >
{

public:

   FtsScheduledCompressorFactory(Rp::System<3,CPT>& system)
    : systemPtr_(&system)
   {}

   Rp::Compressor<3,CPT>* factory(std::string const & className) const
   {
      if (className == "FtsScheduledCompressor") {
         return new FtsScheduledCompressor(*systemPtr_);
      }
      return 0;
   }

private:

   Rp::System<3,CPT>* systemPtr_;

};

/*
* McSimulator with a public compressorFactory() accessor.
*/
class FtsScheduledMcSimulator : public Rp::McSimulator<3,CPT>
{

public:

   FtsScheduledMcSimulator(Rp::System<3,CPT>& system)
    : Rp::McSimulator<3,CPT>(system)
   {}

   using Rp::Simulator<3,CPT>::compressorFactory;

};

/*
* Regression tests (pscfpp fixes): BdMove after a failed BD step.
*/
class FtsBdMoveFailTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Read the system and a McSimulator with an FtsScheduledCompressor.
   */
   void initSimulator(Rp::System<3,CPT>& system,
                      FtsScheduledMcSimulator& simulator,
                      FtsScheduledCompressorFactory& factory)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile("in/Fts/param_system_disordered", in);
      system.readParam(in);
      in.close();

      simulator.compressorFactory().addSubfactory(factory);
      openInputFile("in/Fts/param_McSimulator_BdMove_sched", in);
      simulator.readParam(in);
      in.close();
   }

   /*
   * Read the contents of a file.
   */
   std::string readFile(std::string filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      std::stringstream buffer;
      buffer << in.rdbuf();
      return buffer.str();
   }

   /*
   * A BdMove of two BD steps whose first step fails.
   *
   * The second step continues from the restored state with the old
   * noise of the move (the first noise, of 3 times the usual variance)
   * and a new noise. It adds eta - mobility*dc to the exchange
   * component wc(0), where eta = etaOld + etaNew must have the full
   * Euler-Maruyama variance 2*mobility*meshSize/vSystem at each grid
   * point. If the new noise of the failed step became the old noise of
   * the second step, eta would have half of that variance. The move
   * is accepted, and counted as failed.
   */
   void testFailedStepNoise()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsBdMoveFailedStepNoise.log");

      Rp::System<3,CPT> system;
      FtsScheduledCompressorFactory factory(system);
      FtsScheduledMcSimulator simulator(system);
      initSimulator(system, simulator, factory);
      const double mobility = 1.0E-4; // value in the parameter file

      // Compressed initial state, its wc(0) and force dc(0)
      system.w().readRGrid("in/Fts/w_dis.rf");
      simulator.analyzeChi();
      TEST_ASSERT(simulator.compressor().compress() == 0);
      simulator.computeWc();
      simulator.computeCc();
      simulator.computeDc();

      // eta = wc(0) after - wc(0) before + mobility*dc(0) before
      RField<3,CPT> eta;
      eta.allocate(system.domain().mesh().dimensions());
      VecOp::eqV(eta, simulator.wc(0));
      VecOp::mulEqS(eta, -1.0);
      VecOp::addEqVc(eta, simulator.dc(0), mobility);

      // Calls of compress(): setup, step 1 (fails), step 2
      FtsScheduledCompressor& compressor
         = dynamic_cast<FtsScheduledCompressor&>(simulator.compressor());
      std::vector<bool> fail(3, false);
      fail[1] = true;
      compressor.setSchedule(fail);
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

      // Move statistics
      Rp::McMove<3,CPT>& move = simulator.mcMoveManager()[0];
      TEST_ASSERT(move.nAttempt() == 1);
      TEST_ASSERT(move.nAccept() == 1);
      TEST_ASSERT(move.nFail() == 1);
   }

   /*
   * Statistics of failed BD steps in three moves of two steps each,
   * in which steps 2, 5 and 6 fail. A move is counted as failed if
   * any of its steps failed, and the log reports the BD steps.
   */
   void testFailStatistics()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsBdMoveFailStatistics.log");

      Rp::System<3,CPT> system;
      FtsScheduledCompressorFactory factory(system);
      FtsScheduledMcSimulator simulator(system);
      initSimulator(system, simulator, factory);
      system.w().readRGrid("in/Fts/w_dis.rf");

      // Calls of compress(): setup, then two per move
      FtsScheduledCompressor& compressor
         = dynamic_cast<FtsScheduledCompressor&>(simulator.compressor());
      std::vector<bool> fail(7, false);
      fail[2] = true;
      fail[5] = true;
      fail[6] = true;
      compressor.setSchedule(fail);
      simulator.simulate(3);

      Rp::McMove<3,CPT>& move = simulator.mcMoveManager()[0];
      TEST_ASSERT(move.nAttempt() == 3);
      TEST_ASSERT(move.nAccept() == 3);
      TEST_ASSERT(move.nFail() == 2);

      // Line written by the move at the end of the simulation
      logFile_.flush();
      std::string log = readFile("out/FtsBdMoveFailStatistics.log");
      bool found = (log.find("BdMove: 6 BD steps, 3 failed")
                    != std::string::npos);
      TEST_ASSERT(found);
   }

};

TEST_BEGIN(FtsBdMoveFailTest)
TEST_ADD(FtsBdMoveFailTest, testFailedStepNoise)
TEST_ADD(FtsBdMoveFailTest, testFailStatistics)
TEST_END(FtsBdMoveFailTest)

#endif
