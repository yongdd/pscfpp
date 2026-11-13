#ifndef R1_LINEAR_SWEEP_TEST_H
#define R1_LINEAR_SWEEP_TEST_H

#include "R1TestBase.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <util/misc/Exception.h>
#include <util/misc/FileMaster.h>

#include <cmath>
#include <sstream>
#include <string>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

/*
* Regression tests (pscfpp fixes): LinearSweep parameters and states.
*/
class R1LinearSweepTest : public R1TestBase
{

public:

   /*
   * Sweep a planar blend (nPolymer = 2, with 2 and 1 blocks, nSolvent
   * = 0, nMonomer = 2) with one sweep parameter. Return the message of
   * the Exception thrown, or "" if none.
   */
   std::string sweepMessage(std::string const & parameter)
   {
      std::string text = fileContents("in/R1/linSweep_ids.prm");
      std::string key = "@PARAMETER@";
      std::size_t pos = text.find(key);
      TEST_ASSERT(pos != std::string::npos);
      text.replace(pos, key.size(), parameter);

      System sys;
      std::istringstream in(text);
      sys.readParam(in);
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      std::string message;
      try {
         runCommands(sys, "READ_W in/R1/linSweep_blend.w\nSWEEP\nFINISH\n");
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
   * Out-of-range indices of sweep parameters stop the sweep with a
   * message that names the index. They were not checked, and stopped
   * on unrelated assertions, or (negative and block indices, last)
   * read and wrote past the arrays. MixtureTmpl::polymer(-1) must
   * throw.
   */
   void testIdsOutOfRange()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1LinearSweepIds.log");

      TEST_ASSERT(sweepMessage("chi 0 1 0.5") == "");
      TEST_ASSERT(contains(sweepMessage("phi_polymer 2 0.0"),
                           "polymer index 2 out of range"));
      TEST_ASSERT(contains(sweepMessage("solvent 0 0.1"),
                           "solvent index 0 out of range"));
      TEST_ASSERT(contains(sweepMessage("kuhn 2 0.1"),
                           "monomer index 2 out of range"));
      TEST_ASSERT(contains(sweepMessage("chi 0 2 1.0"),
                           "monomer index 2 out of range"));
      TEST_ASSERT(contains(sweepMessage("phi_polymer -1 0.0"),
                           "polymer index -1 out of range"));
      TEST_ASSERT(contains(sweepMessage("block 1 1 0.1"),
                           "block index 1 out of range"));
      TEST_ASSERT(contains(sweepMessage("block 0 2 0.1"),
                           "block index 2 out of range"));

      // Negative species index
      System sys;
      readSystem(sys, "in/R1/linSweep_end.prm");
      bool thrown = false;
      try {
         sys.mixture().polymer(-1);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(thrown);
   }

   /*
   * A sweep that removes a homopolymer (phi 0.1 -> 0) with ns = 9.
   * Repeated addition of ds = 1/9 gave a last s = 1 + 2.2e-16, so
   * phi = -2.8e-17 at the last state, and the sweep stopped there.
   */
   void testEndsAtOne()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1LinearSweepEnd.log");

      System sys;
      readSystem(sys, "in/R1/linSweep_end.prm");
      bool thrown = false;
      try {
         runCommands(sys, "READ_W in/R1/linSweep_blend.w\nSWEEP\nFINISH\n");
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
      TEST_ASSERT(sys.mixture().polymer(0).phi() == 1.0);
      TEST_ASSERT(sys.mixture().polymer(1).phi() == 0.0);
   }

   /*
   * In a closed blend, phi changes +0.05 and -0.0499995 cancel only to
   * 5e-7. The volume fractions must be rescaled to add up to 1: they
   * added up to 1 + 5e-7, and the iterator failed to converge at s > 0
   * until the sweep stopped ("Sweep decreased ds too many times").
   */
   void testClosedPhiSum()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1LinearSweepSum.log");

      System sys;
      readSystem(sys, "in/R1/linSweep_sum.prm");
      bool thrown = false;
      try {
         runCommands(sys, "READ_W in/R1/linSweep_blend.w\nSWEEP\nFINISH\n");
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);

      double phi0 = sys.mixture().polymer(0).phi();
      double phi1 = sys.mixture().polymer(1).phi();
      Log::file() << "phi0 + phi1 - 1 = " << phi0 + phi1 - 1.0
                  << std::endl;
      TEST_ASSERT(std::abs(phi0 + phi1 - 1.0) < 1.0E-14);
      TEST_ASSERT(std::abs(phi0 - 0.95/1.0000005) < 1.0E-13);
   }

   void testSummaryRows()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1LinearSweepSummary.log");
      TEST_ASSERT(sweepMessage("chi 0 1 0.5") == "");
      std::istringstream in(fileContents("out/R1LinearSweepIds_log"));
      int index, count = 0;
      double coordinate, freeEnergy, pressure;
      while (in >> index >> coordinate >> freeEnergy >> pressure) {
         TEST_ASSERT(index == count);
         TEST_ASSERT(std::abs(coordinate - 0.5*count) < 1.0E-12);
         TEST_ASSERT(std::isfinite(freeEnergy) && std::isfinite(pressure));
         ++count;
      }
      TEST_ASSERT(count == 3);
   }

};

TEST_BEGIN(R1LinearSweepTest)
TEST_ADD(R1LinearSweepTest, testIdsOutOfRange)
TEST_ADD(R1LinearSweepTest, testEndsAtOne)
TEST_ADD(R1LinearSweepTest, testClosedPhiSum)
TEST_ADD(R1LinearSweepTest, testSummaryRows)
TEST_END(R1LinearSweepTest)

#endif
