#ifndef RPG_FTS_STRUCTURE_FACTOR_TEST_H
#define RPG_FTS_STRUCTURE_FACTOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
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

};

TEST_BEGIN(FtsStructureFactorTest)
TEST_ADD(FtsStructureFactorTest, testWaveBunchesRectangular)
TEST_END(FtsStructureFactorTest)

#endif
