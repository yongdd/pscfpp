#ifndef RPC_FTS_ANALYZE_CHI_TEST_H
#define RPC_FTS_ANALYZE_CHI_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>

#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <string>

#if defined(__GLIBC__)
#include <malloc.h>
#endif

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): Simulator::analyzeChi frees its GSL
* workspace.
*/
class FtsAnalyzeChiTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   #if defined(__GLIBC__)
   /*
   * Number of bytes allocated with malloc and not yet freed.
   */
   long heapInUse()
   {
      #if (__GLIBC__ > 2) || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 33)
      struct mallinfo2 info = mallinfo2();
      #else
      struct mallinfo info = mallinfo();
      #endif
      return (long) info.uordblks;
   }
   #endif

   /*
   * analyzeChi() runs at every step of a ramp that changes chi. Each
   * call must free the memory it allocates (with the GSL leak, 1000
   * calls for 2 monomers left about 370 kB allocated).
   */
   void testNoLeak()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsAnalyzeChiNoLeak.log");

      Rp::System<3,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile("in/Fts/param_system_small", in);
      system.readParam(in);
      in.close();

      Rp::BdSimulator<3,CPT> simulator(system);
      openInputFile("in/Fts/param_BdSimulator_analyzeChi", in);
      simulator.readParam(in);
      in.close();

      #if defined(__GLIBC__)
      simulator.analyzeChi();
      long before = heapInUse();
      for (int i = 0; i < 1000; ++i) {
         simulator.analyzeChi();
      }
      long growth = heapInUse() - before;
      Log::file() << "Heap growth in 1000 calls = " << growth 
                  << " bytes" << std::endl;
      TEST_ASSERT(growth < 10000);
      #else
      Log::file() << "Heap use not measured (no glibc)" << std::endl;
      #endif
   }

};

TEST_BEGIN(FtsAnalyzeChiTest)
TEST_ADD(FtsAnalyzeChiTest, testNoLeak)
TEST_END(FtsAnalyzeChiTest)

#endif
