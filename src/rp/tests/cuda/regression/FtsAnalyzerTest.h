#ifndef RPG_FTS_ANALYZER_TEST_H
#define RPG_FTS_ANALYZER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/fts/analyzer/AnalyzerManager.h>
#include <rp/fts/trajectory/RGridTrajectoryReader.h>
#include <rp/field/WFields.h>
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
* Regression tests (pscfpp fixes): Analyzers, trajectory files and ANALYZE.
*/
class FtsAnalyzerTest : public LogFileUnitTest
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
   * Return the number of nonempty lines in a file.
   */
   int countLines(std::string filename)
   {
      std::ifstream file;
      openInputFile(filename, file);
      std::string line;
      int n = 0;
      while (std::getline(file, line)) {
         if (line.find_first_not_of(" \t\r") != std::string::npos) {
            ++n;
         }
      }
      file.close();
      return n;
   }

   /*
   * Two SIMULATE runs with averaging analyzers that write .dat files.
   *
   * The .dat files of an AverageListAnalyzer (HamiltonianAnalyzer) and
   * an AverageAnalyzer (MaxOrderParameter) must hold the samples of
   * both runs.
   */
   void testSimulateTwice()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsSimulateTwice.log");

      Rp::System<3,CUT> system;
      initSystem(system, "in/Fts/param_system_small");
      Rp::BdSimulator<3,CUT> simulator(system);
      initSimulator(simulator, "in/Fts/param_BdSimulator_twice");
      system.w().readRGrid("in/Fts/w_small.rf");

      bool success = true;
      try {
         simulator.simulate(3);
         simulator.simulate(3);
      } catch (Exception& e) {
         success = false;
      }
      TEST_ASSERT(success);

      // Steps 0 to 3 of each run, one line per sample
      TEST_ASSERT(countLines("out/FtsTwice_hamiltonian.dat") == 8);
      TEST_ASSERT(countLines("out/FtsTwice_maxOrder.dat") == 8);
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
   * Run ANALYZE with an RGridTrajectoryReader.
   *
   * Return false if it throws an Exception.
   */
   template <class SimulatorT>
   bool analyze(SimulatorT& simulator, int min, int max,
                std::string filename)
   {
      try {
         simulator.analyze(min, max, "RGridTrajectoryReader",
                           filePrefix() + filename);
      } catch (Exception& e) {
         return false;
      }
      return true;
   }

   /*
   * Copy a trajectory file, ending the copy after the first nLine
   * lines of frame iFrame (a damaged last frame).
   */
   void writeTruncatedTrajectory(std::string inFileName,
                                 std::string outFileName,
                                 int iFrame, int nLine)
   {
      std::ifstream in;
      openInputFile(inFileName, in);
      std::ofstream out;
      openOutputFile(outFileName, out);
      std::ostringstream label;
      label << "i = " << iFrame;
      std::string line;
      int n = -1;
      while (std::getline(in, line)) {
         if (line == label.str()) n = 0;
         if (n >= 0) {
            if (n == nLine) break;
            ++n;
         }
         out << line << "\n";
      }
      in.close();
      out.close();
   }

   /*
   * ANALYZE with baseInterval 2 (analyzer interval 2) of frames 0-3:
   * frames 0 and 2 are sampled.
   */
   template <class SimulatorT>
   void analyzeBaseInterval(std::string paramFileName,
                            std::string dataFileName)
   {
      Rp::System<3,CUT> system;
      initSystem(system, "in/Fts/param_system_small");
      SimulatorT simulator(system);
      initSimulator(simulator, paramFileName);
      TEST_ASSERT(analyze(simulator, 0, 3, "in/Fts/traj_small.rf"));
      TEST_ASSERT(countLines(dataFileName) == 2);
   }

   void testAnalyzeBaseInterval()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsAnalyzeBaseInterval.log");
      analyzeBaseInterval< Rp::BdSimulator<3,CUT> >(
                              "in/Fts/param_BdSimulator_analyzeBase",
                              "out/FtsAnalyzeBaseBd_hamiltonian.dat");
      analyzeBaseInterval< Rp::McSimulator<3,CUT> >(
                              "in/Fts/param_McSimulator_analyzeBase",
                              "out/FtsAnalyzeBaseMc_hamiltonian.dat");
   }

   /*
   * ANALYZE 0 2 reads frames 0 to 2 only. The system is left with
   * frame 2, and a damaged frame 3 does not matter.
   */
   template <class SimulatorT>
   void analyzeLastFrame(std::string paramFileName)
   {
      // Frame 2, read directly
      DArray< RField<3,CUT> > w2;
      {
         Rp::System<3,CUT> system;
         initSystem(system, "in/Fts/param_system_small");
         Rp::RGridTrajectoryReader<3,CUT> reader(system);
         reader.open(filePrefix() + "in/Fts/traj_small.rf");
         reader.readHeader();
         for (int i = 0; i <= 2; ++i) {
            TEST_ASSERT(reader.readFrame());
         }
         reader.close();
         copyFields(system.w().rgrid(), w2);
      }

      // Complete trajectory, then the one with a damaged frame 3
      std::string fileNames[2] = {"in/Fts/traj_small.rf",
                                  "out/FtsTrajDamaged.rf"};
      for (int j = 0; j < 2; ++j) {
         Rp::System<3,CUT> system;
         initSystem(system, "in/Fts/param_system_small");
         SimulatorT simulator(system);
         initSimulator(simulator, paramFileName);
         TEST_ASSERT(analyze(simulator, 0, 2, fileNames[j]));
         RFieldComparison<3,CUT> comparison;
         comparison.compare(w2, system.w().rgrid());
         TEST_ASSERT(comparison.maxDiff() < 1.0E-12);
      }
   }

   void testAnalyzeLastFrame()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsAnalyzeLastFrame.log");
      writeTruncatedTrajectory("in/Fts/traj_small.rf",
                               "out/FtsTrajDamaged.rf", 3, 100);
      analyzeLastFrame< Rp::BdSimulator<3,CUT> >(
                              "in/Fts/param_BdSimulator_analyze");
      analyzeLastFrame< Rp::McSimulator<3,CUT> >(
                              "in/Fts/param_McSimulator_analyze");
   }

   /*
   * ANALYZE with no frame in [min, max] (the trajectory has frames
   * 0-3) must fail.
   */
   template <class SimulatorT>
   void analyzeNoFrame(std::string paramFileName)
   {
      Rp::System<3,CUT> system;
      initSystem(system, "in/Fts/param_system_small");
      SimulatorT simulator(system);
      initSimulator(simulator, paramFileName);
      TEST_ASSERT(!analyze(simulator, 10, 20, "in/Fts/traj_small.rf"));
   }

   void testAnalyzeNoFrame()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsAnalyzeNoFrame.log");
      analyzeNoFrame< Rp::BdSimulator<3,CUT> >(
                              "in/Fts/param_BdSimulator_analyze");
      analyzeNoFrame< Rp::McSimulator<3,CUT> >(
                              "in/Fts/param_McSimulator_analyze");
   }

   /*
   * ANALYZE of frames 1-3 with an analyzer of interval 4, which
   * samples none of them, must fail.
   */
   template <class SimulatorT>
   void analyzeNoSample(std::string paramFileName)
   {
      Rp::System<3,CUT> system;
      initSystem(system, "in/Fts/param_system_small");
      SimulatorT simulator(system);
      initSimulator(simulator, paramFileName);
      TEST_ASSERT(!analyze(simulator, 1, 3, "in/Fts/traj_small.rf"));
   }

   void testAnalyzeNoSample()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsAnalyzeNoSample.log");
      analyzeNoSample< Rp::BdSimulator<3,CUT> >(
                              "in/Fts/param_BdSimulator_analyzeInterval");
      analyzeNoSample< Rp::McSimulator<3,CUT> >(
                              "in/Fts/param_McSimulator_analyzeInterval");
   }

   /*
   * Timers of an analysis-only simulator, which has no Compressor.
   */
   template <class SimulatorT>
   void timersWithoutCompressor(std::string paramFileName)
   {
      Rp::System<3,CUT> system;
      initSystem(system, "in/Fts/param_system_small");
      SimulatorT simulator(system);
      initSimulator(simulator, paramFileName);
      TEST_ASSERT(!simulator.hasCompressor());
      TEST_ASSERT(analyze(simulator, 0, 3, "in/Fts/traj_small.rf"));

      bool success = true;
      try {
         std::ostringstream out;
         simulator.outputTimers(out);
         simulator.clearTimers();
      } catch (Exception& e) {
         success = false;
      }
      TEST_ASSERT(success);
   }

   void testTimersWithoutCompressor()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsTimersWithoutCompressor.log");
      timersWithoutCompressor< Rp::BdSimulator<3,CUT> >(
                              "in/Fts/param_BdSimulator_analyze");
      timersWithoutCompressor< Rp::McSimulator<3,CUT> >(
                              "in/Fts/param_McSimulator_analyze");
   }

};

TEST_BEGIN(FtsAnalyzerTest)
TEST_ADD(FtsAnalyzerTest, testSimulateTwice)
TEST_ADD(FtsAnalyzerTest, testAnalyzeBaseInterval)
TEST_ADD(FtsAnalyzerTest, testAnalyzeLastFrame)
TEST_ADD(FtsAnalyzerTest, testAnalyzeNoFrame)
TEST_ADD(FtsAnalyzerTest, testAnalyzeNoSample)
TEST_ADD(FtsAnalyzerTest, testTimersWithoutCompressor)
TEST_END(FtsAnalyzerTest)

#endif
