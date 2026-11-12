#ifndef R1_COMMAND_TEST_H
#define R1_COMMAND_TEST_H

#include "R1TestBase.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/field/Domain.h>

#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

class R1CommandTest : public R1TestBase
{

public:

   /*
   * Commands with arguments that are missing or cannot be read.
   */
   void testBadArgument()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1CommandTestBadArgument.log");

      System sys;
      readSystem(sys, "in/R1/planar_nr2.prm");
      runCommands(sys, "READ_W in/R1/planar.w\n"
                       "WRITE_W out/badArgumentRef.w\nFINISH\n");
      std::string ref = fileContents("out/badArgumentRef.w");

      // Each command file writes out/badArgument.w, then gives a command
      // with a bad argument. The command loop must throw an Exception
      // without executing that command, which would otherwise use the
      // file name of the previous command.
      const int nCase = 5;
      std::string commands[nCase] = {
         // Last arguments missing at the end of the command file
         "REMESH_W 201\n",
         "WRITE_C\n",
         "COMPARE_HOMOGENEOUS\n",
         // Arguments that cannot be read
         "REMESH_W out/badArgumentRemesh.w 201\nFINISH\n",
         "COMPARE_HOMOGENEOUS out/badArgument.homo\nFINISH\n"
      };
      for (int i = 0; i < nCase; ++i) {
         Log::file() << "\nCommand: " << commands[i] << std::endl;
         bool thrown = false;
         try {
            runCommands(sys, "WRITE_W out/badArgument.w\n" + commands[i]);
         } catch (Exception& e) {
            thrown = true;
         }
         TEST_ASSERT(thrown);
         TEST_ASSERT(fileContents("out/badArgument.w") == ref);
      }
   }

   /*
   * Command file whose last command is not followed by a newline.
   */
   void testNoFinalNewline()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1CommandTestNoFinalNewline.log");

      System sys;
      readSystem(sys, "in/R1/planar_nr2.prm");
      int nx = sys.domain().nx();
      for (int i = 0; i < nx; ++i) {
         sys.cField(0)[i] = 0.0;
         sys.cField(1)[i] = 0.0;
      }

      // COMPUTE ends the command file, without FINISH or a newline
      runCommands(sys, "READ_W in/R1/planar.w\nCOMPUTE");

      // The concentrations were computed: they add up to 1 on average
      double sum = sys.domain().spatialAverage(sys.cField(0))
                 + sys.domain().spatialAverage(sys.cField(1));
      Log::file() << "average of c(0) + c(1) = " << sum << std::endl;
      TEST_ASSERT(std::abs(sum - 1.0) < 1.0E-6);
   }

   /*
   * Free energy and pressure after COMPUTE.
   */
   void testComputeThermo()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1CommandTestComputeThermo.log");

      System sys;
      readSystem(sys, "in/R1/planar_nr2.prm");

      // Converged solution
      runCommands(sys, "READ_W in/R1/planar.w\nITERATE\nFINISH\n");
      double fIterate = sys.fHelmholtz();

      // MDE solution for the initial fields, without iteration
      runCommands(sys, "READ_W in/R1/planar.w\nCOMPUTE\n"
                       "WRITE_THERMO out/computeThermo\nFINISH\n");
      double f = sys.fHelmholtz();
      double p = sys.pressure();

      // Values computed explicitly for the current fields
      sys.computeFreeEnergy();
      Log::file() << "fHelmholtz: ITERATE " << fIterate
                  << ", COMPUTE " << f
                  << ", current fields " << sys.fHelmholtz() << std::endl;
      TEST_ASSERT(std::abs(sys.fHelmholtz() - fIterate) > 1.0E-3);
      TEST_ASSERT(eq(f, sys.fHelmholtz()));
      TEST_ASSERT(eq(p, sys.pressure()));
   }

};

TEST_BEGIN(R1CommandTest)
TEST_ADD(R1CommandTest, testBadArgument)
TEST_ADD(R1CommandTest, testNoFinalNewline)
TEST_ADD(R1CommandTest, testComputeThermo)
TEST_END(R1CommandTest)

#endif
