#ifndef RPG_SCFT_LINEAR_SWEEP_TEST_H
#define RPG_SCFT_LINEAR_SWEEP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <pscf/crystal/UnitCell.h>

#include <util/misc/Exception.h>
#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): LinearSweep parameters and states.
*/
class ScftLinearSweepTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Read a System from a parameter file in which the strings
   * @FLEXIBLE@ and @PARAMETER@ are replaced by the given values.
   */
   void readSystem(Rp::System<1,CUT>& system, std::string filename,
                   std::string flexible = "0",
                   std::string parameter = "")
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream file;
      openInputFile(filename, file);
      std::stringstream buffer;
      buffer << file.rdbuf();
      file.close();
      std::string text = buffer.str();
      replace(text, "@FLEXIBLE@", flexible);
      replace(text, "@PARAMETER@", parameter);

      std::istringstream in(text);
      system.readParam(in);
   }

   /*
   * Replace the first occurrence of key in text by value.
   */
   void replace(std::string& text, std::string const & key,
                std::string const & value)
   {
      std::size_t pos = text.find(key);
      if (pos != std::string::npos) {
         text.replace(pos, key.size(), value);
      }
   }

   /*
   * Sweep a lamellar diblock (nPolymer = 1, nBlock = 2, nSolvent = 0,
   * nMonomer = 2, one lattice parameter) with one sweep parameter.
   * Return the message of the Exception thrown, or "" if none.
   */
   std::string sweepMessage(std::string const & parameter,
                            std::string const & flexible = "0")
   {
      Rp::System<1,CUT> system;
      readSystem(system, "in/Scft/linSweep_lam", flexible, parameter);
      system.w().readBasis("in/Scft/linSweep_lam_w.bf");
      std::string message;
      try {
         system.sweep();
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
   * message that names the index. They were not checked: block and
   * negative polymer indices read and wrote past the arrays, a lattice
   * parameter index gave a sweep that changed nothing, and the others
   * stopped on unrelated assertions. Negative and block indices come
   * last, since they may crash without the check.
   */
   void testIdsOutOfRange()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftLinearSweepIds.log");

      TEST_ASSERT(sweepMessage("block 0 1 +0.02") == "");
      TEST_ASSERT(contains(sweepMessage("mu_polymer 1 0.1"),
                           "polymer index 1 out of range"));
      TEST_ASSERT(contains(sweepMessage("block 1 0 +0.02"),
                           "polymer index 1 out of range"));
      TEST_ASSERT(contains(sweepMessage("phi_solvent 0 0.0"),
                           "solvent index 0 out of range"));
      TEST_ASSERT(contains(sweepMessage("kuhn 2 0.1"),
                           "monomer index 2 out of range"));
      TEST_ASSERT(contains(sweepMessage("chi 0 2 1.0"),
                           "monomer index 2 out of range"));
      TEST_ASSERT(contains(sweepMessage("cell_param 1 +0.1"),
                           "lattice parameter index 1 out of range"));
      TEST_ASSERT(contains(sweepMessage("phi_polymer -1 0.0"),
                           "polymer index -1 out of range"));
      TEST_ASSERT(contains(sweepMessage("cell_param -1 +0.1"),
                           "lattice parameter index -1 out of range"));
      TEST_ASSERT(contains(sweepMessage("block 0 2 +0.02"),
                           "block index 2 out of range"));
      TEST_ASSERT(contains(sweepMessage("block -1 0 +0.02"),
                           "polymer index -1 out of range"));
   }

   /*
   * A sweep that removes a homopolymer (phi 0.1 -> 0) with ns = 9.
   * Repeated addition of ds = 1/9 gave a last s = 1 + 2.2e-16, so
   * phi = -2.8e-17 at the last state, and the sweep stopped there.
   * The last state must be at s = 1, with the exact final values.
   */
   void testEndsAtOne()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftLinearSweepEnd.log");

      std::string lastFile = filePrefix()
                           + "out/ScftLinearSweepEnd_9.stt";
      std::remove(lastFile.c_str());

      Rp::System<1,CUT> system;
      readSystem(system, "in/Scft/linSweep_blend_end");
      system.w().readBasis("in/Scft/linSweep_blend_w.bf");
      bool thrown = false;
      try {
         system.sweep();
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
      TEST_ASSERT(system.mixture().polymer(0).phi() == 1.0);
      TEST_ASSERT(system.mixture().polymer(1).phi() == 0.0);

      // The last state (s = 1) was written
      std::ifstream f;
      f.open(lastFile.c_str());
      TEST_ASSERT(f.good());
   }

   /*
   * In a closed blend, phi changes +0.05 and -0.0499993 cancel only to
   * 7e-7 (accepted by readParameters). The volume fractions must be
   * rescaled to add up to 1: they added up to 1 + 7e-7, and the SCFT
   * iterator converged to a mean total concentration of 1 + 7e-7.
   */
   void testClosedPhiSum()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftLinearSweepSum.log");

      Rp::System<1,CUT> system;
      readSystem(system, "in/Scft/linSweep_blend_sum");
      system.w().readBasis("in/Scft/linSweep_blend_w.bf");
      system.sweep();

      double phi0 = system.mixture().polymer(0).phi();
      double phi1 = system.mixture().polymer(1).phi();
      if (verbose() > 0) {
         std::cout << "\nphi0 + phi1 - 1 = " << phi0 + phi1 - 1.0;
      }
      TEST_ASSERT(std::abs(phi0 + phi1 - 1.0) < 1.0E-14);
      TEST_ASSERT(std::abs(phi0 - 0.95/1.0000007) < 1.0E-13);

      // Mean total concentration (coefficient of basis function 0)
      double cTot = system.c().basis(0)[0] + system.c().basis(1)[0];
      if (verbose() > 0) {
         std::cout << "\nmean c - 1 = " << cTot - 1.0;
      }
      TEST_ASSERT(std::abs(cTot - 1.0) < 1.0E-10);
   }

   /*
   * A cell_param sweep of a flexible lattice parameter stops with a
   * message: the iterator optimized the cell and the extrapolation
   * replaced the swept value, so every state was the same cell. A
   * sweep of the same parameter with a rigid cell changes the cell.
   */
   void testCellParamFlexible()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftLinearSweepCell.log");

      // Rigid cell
      {
         Rp::System<1,CUT> system;
         readSystem(system, "in/Scft/linSweep_lam", "0",
                    "cell_param 0 +0.1");
         system.w().readBasis("in/Scft/linSweep_lam_w.bf");
         double a0 = system.domain().unitCell().parameter(0);
         system.sweep();
         double a1 = system.domain().unitCell().parameter(0);
         TEST_ASSERT(std::abs(a1 - a0 - 0.1) < 1.0E-12);
      }

      // Flexible cell
      std::string message = sweepMessage("cell_param 0 +0.1", "1");
      TEST_ASSERT(contains(message, "lattice parameter is flexible"));
   }

};

TEST_BEGIN(ScftLinearSweepTest)
TEST_ADD(ScftLinearSweepTest, testIdsOutOfRange)
TEST_ADD(ScftLinearSweepTest, testEndsAtOne)
TEST_ADD(ScftLinearSweepTest, testClosedPhiSum)
TEST_ADD(ScftLinearSweepTest, testCellParamFlexible)
TEST_END(ScftLinearSweepTest)

#endif
