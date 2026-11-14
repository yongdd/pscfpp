#ifndef R1_ITERATOR_TEST_H
#define R1_ITERATOR_TEST_H

#include "R1TestBase.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/field/Domain.h>
#include <r1/iterator/Iterator.h>
#include <pscf/interaction/Interaction.h>

#include <sstream>
#include <cmath>
#include <limits>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

class R1IteratorTest : public R1TestBase
{

public:

   /*
   * NrIterator on fields that give a NaN residual.
   */
   void testNrNan()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1IteratorTestNrNan.log");

      System sys;
      readSystem(sys, "in/R1/planar_nr2.prm");
      runCommands(sys, "READ_W in/R1/planar.w\nFINISH\n");

      // A NaN at one grid point makes the propagators, concentrations
      // and residuals NaN
      int nx = sys.domain().nx();
      sys.wField(0)[nx/2] = std::numeric_limits<double>::quiet_NaN();

      // The iterator must report a failure, not convergence
      int error = sys.iterator().solve();
      TEST_ASSERT(error != 0);
   }

   /*
   * Solve, and return the error code of the iterator, or -1 if the
   * solve threw an Exception.
   */
   int trySolve(System& sys)
   {
      try {
         return sys.iterator().solve();
      } catch (Exception const &) {
         // The exception has written its message to the log file
         return -1;
      }
   }

   /*
   * Check that the iterator of a parameter file reports a failed solve
   * (error code 1, not an exception) for a field that is not finite at
   * one grid point, and converges again from valid fields afterwards.
   */
   void checkNonfiniteFields(char const * paramFile, char const * logFile)
   {
      openLogFile(logFile);
      System sys;
      readSystem(sys, paramFile);
      const int nx = sys.domain().nx();
      const double infinity = std::numeric_limits<double>::infinity();
      const double invalid[] = {infinity, -infinity,
                               std::numeric_limits<double>::quiet_NaN()};
      // An interior point, and the canonical gauge reference (last
      // element), which the iterators subtract from all fields
      const int points[] = {nx/2, nx-1};
      for (double value : invalid) {
         for (int monomer = 0; monomer < 2; ++monomer) {
            for (int point : points) {
               runCommands(sys, "READ_W in/R1/planar.w\nFINISH\n");
               sys.wField(monomer)[point] = value;
               TEST_ASSERT(trySolve(sys) == 1);
            }
         }
      }
      runCommands(sys, "READ_W in/R1/planar.w\nFINISH\n");
      TEST_ASSERT(trySolve(sys) == 0);
   }

   void testNrNonfiniteFieldsAndRecovery()
   {
      printMethod(TEST_FUNC);
      checkNonfiniteFields("in/R1/planar_nr2.prm",
                           "out/R1IteratorTestNonfinite.log");
   }

   void testAmNonfiniteFieldsAndRecovery()
   {
      printMethod(TEST_FUNC);
      checkNonfiniteFields("in/R1/planar_am.prm",
                           "out/R1IteratorTestAmNonfinite.log");
   }

   void testBinaryRelaxNonfiniteFieldsAndRecovery()
   {
      printMethod(TEST_FUNC);
      checkNonfiniteFields("in/R1/planar_br.prm",
                           "out/R1IteratorTestBinaryRelaxNonfinite.log");
   }

   /*
   * BinaryRelaxIterator with a nonzero diagonal chi element.
   */
   void testBinaryRelaxDiagonalChi()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1IteratorTestBinaryRelaxDiagonalChi.log");

      // Diblock with chi(0,0) = 6 and chi(0,1) = 20
      System sys;
      readSystem(sys, "in/R1/planar_br_diag.prm");
      runCommands(sys, "READ_W in/R1/planar.w\nFINISH\n");

      int error = sys.iterator().solve();
      TEST_ASSERT(error == 0);

      // SCF equations: w(i) = sum_j chi(i,j) c(j) + xi, c(0) + c(1) = 1
      Interaction const & interaction = sys.interaction();
      double chi00 = interaction.chi(0,0);
      double chi01 = interaction.chi(0,1);
      double chi11 = interaction.chi(1,1);
      double c0, c1, r;
      double maxExchange = 0.0;
      double maxIncompress = 0.0;
      int nx = sys.domain().nx();
      for (int i = 0; i < nx; ++i) {
         c0 = sys.cField(0)[i];
         c1 = sys.cField(1)[i];
         r = sys.wField(1)[i] - sys.wField(0)[i]
             - (chi01 - chi00)*c0 - (chi11 - chi01)*c1;
         if (std::abs(r) > maxExchange) maxExchange = std::abs(r);
         r = c0 + c1 - 1.0;
         if (std::abs(r) > maxIncompress) maxIncompress = std::abs(r);
      }
      Log::file() << "max exchange residual = " << maxExchange << "\n"
                  << "max incompressibility residual = "
                  << maxIncompress << std::endl;
      TEST_ASSERT(maxExchange < 1.0E-6);
      TEST_ASSERT(maxIncompress < 1.0E-6);
   }

   /*
   * BinaryRelaxIterator after a solve that reduced its step size.
   */
   void testBinaryRelaxLambda()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1IteratorTestBinaryRelaxLambda.log");

      // lambdaPlus = lambdaMinus = 6 is too large at the start of this
      // solve, so the iterator reduces its step
      System sys;
      readSystem(sys, "in/R1/planar_br.prm");
      runCommands(sys, "READ_W in/R1/planar.w\nFINISH\n");
      int error = sys.iterator().solve();
      TEST_ASSERT(error == 0);

      // The parameters written afterwards are those that were read
      std::ostringstream out;
      sys.iterator().writeParam(out);
      Log::file() << out.str();
      TEST_ASSERT(eq(paramValue(out.str(), "lambdaPlus"), 6.0));
      TEST_ASSERT(eq(paramValue(out.str(), "lambdaMinus"), 6.0));
   }

};

TEST_BEGIN(R1IteratorTest)
TEST_ADD(R1IteratorTest, testNrNan)
TEST_ADD(R1IteratorTest, testNrNonfiniteFieldsAndRecovery)
TEST_ADD(R1IteratorTest, testAmNonfiniteFieldsAndRecovery)
TEST_ADD(R1IteratorTest, testBinaryRelaxNonfiniteFieldsAndRecovery)
TEST_ADD(R1IteratorTest, testBinaryRelaxDiagonalChi)
TEST_ADD(R1IteratorTest, testBinaryRelaxLambda)
TEST_END(R1IteratorTest)

#endif
