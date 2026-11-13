#ifndef RPG_FTS_STRUCTURE_FACTOR_TEST_H
#define RPG_FTS_STRUCTURE_FACTOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/analyzer/BinaryStructureFactor.h>
#include <rp/field/WFields.h>

#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): BinaryStructureFactor wave bunches.
*/
class FtsStructureFactorTest : public LogFileUnitTest
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
   * Rectangular cell 3.2 x 3.2000000064 on a 4 x 4 mesh. Waves that are
   * equivalent in a square cell have k^2 values 1.5e-8 apart (separate
   * bunches of the WaveList) but |k| values only 4e-9 apart. The
   * BinaryStructureFactor setup (SIMULATE) must accept this cell.
   */
   void testWaveBunchesRectangular()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsWaveBunchesRectangular.log");

      Rp::System<2,CUT> system;
      initSystem(system, "in/Fts/param_system_rect_2D");
      Rp::BdSimulator<2,CUT> simulator(system);
      initSimulator(simulator, "in/Fts/param_BdSimulator_structureFactor");
      system.w().readRGrid("in/Fts/w_rect_2D.rf");

      bool success = true;
      try {
         simulator.simulate(2);
      } catch (Exception& e) {
         success = false;
      }
      TEST_ASSERT(success);

      // S(q) output: one line per bunch, in order of wavenumber
      // (wavenumbers written with 8 significant digits)
      std::ifstream file;
      openInputFile("out/FtsStructureFactor", file);
      std::vector<double> q;
      double x, s;
      while (file >> x >> s) {
         q.push_back(x);
      }
      file.close();
      TEST_ASSERT(q.size() == 9);
      TEST_ASSERT(std::abs(q[0]) < 1.0E-10);
      TEST_ASSERT(q[1] > 1.0);
      for (size_t i = 1; i < q.size(); ++i) {
         TEST_ASSERT(q[i] >= q[i-1]);
      }
   }

   void testSampleCell()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsSampleCell.log");
      Rp::System<3,CUT> system;
      initSystem(system, "in/Fts/param_system_small");
      Rp::BdSimulator<3,CUT> simulator(system);
      initSimulator(simulator, "in/Fts/param_BdSimulator_structureFactor");
      system.w().readRGrid("in/Fts/w_small.rf");
      simulator.analyzeChi();
      Rp::BinaryStructureFactor<3,CUT> analyzer(simulator, system);
      std::istringstream param("BinaryStructureFactor{ interval 2 "
                              "outputFileName out/FtsSampleCell }");
      analyzer.readParam(param);
      analyzer.setup();
      FSArray<double, 6> cell = system.domain().unitCell().parameters();
      cell[0] = 4.0;
      system.setUnitCell(cell);
      analyzer.sample(1); // Not sampled: changing the cell is still permitted.
      analyzer.sample(2);
      analyzer.output();
      std::ifstream file;
      openInputFile("out/FtsSampleCell", file);
      double q, value;
      file >> q >> value;
      TEST_ASSERT(std::abs(q) < 1.0E-10);
      file >> q >> value;
      TEST_ASSERT(std::abs(q - 0.5*std::acos(-1.0)) < 1.0E-7);
      file.close();
      analyzer.sample(4); // Same cell: accumulate normally.
      cell[0] = 2.0;
      system.setUnitCell(cell);
      bool rejected = false;
      try { analyzer.sample(6); }
      catch (Exception& e) {
         rejected = e.message().find("same unit cell") != std::string::npos;
      }
      TEST_ASSERT(rejected);
      // A separate accumulation can use a different cell.
      analyzer.setup();
      analyzer.sample(8);
      analyzer.output();
      openInputFile("out/FtsSampleCell", file);
      file >> q >> value >> q >> value;
      TEST_ASSERT(std::abs(q - std::acos(-1.0)) < 1.0E-7);
   }

};

TEST_BEGIN(FtsStructureFactorTest)
TEST_ADD(FtsStructureFactorTest, testWaveBunchesRectangular)
TEST_ADD(FtsStructureFactorTest, testSampleCell)
TEST_END(FtsStructureFactorTest)

#endif
