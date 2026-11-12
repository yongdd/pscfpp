#ifndef RPC_COMMAND_TEST_H
#define RPC_COMMAND_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <rp/field/Mask.h>
#include <rp/solvers/Mixture.h>
#include <rp/environment/Environment.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/compressor/Compressor.h>

#include <pscf/chem/PolymerModel.h>
#include <pscf/field/cpp/RField.h>

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
* Regression tests of the command script loop, System::readCommands.
*/
class CommandTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * ANALYZE in a system without a Simulator throws an Exception.
   */
   void testAnalyzeNoSimulator()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testAnalyzeNoSimulator.log");

      Rp::System<1,CPT> system;
      initSystem(system, "in/Command/lam_param", "in/Command/lam_w.bf");
      TEST_ASSERT(!system.hasSimulator());

      TEST_ASSERT(readCommandsThrows(system,
                  "ANALYZE 0 10 HamiltonianAnalyzer out/none\nFINISH\n"));
   }

   /*
   * COMPUTE in a film system without an Iterator.
   */
   void testFilmNoIterator()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testFilmNoIterator.log");

      Rp::System<1,CPT> system;
      initSystem(system, "in/Command/film_param", "in/Command/film_w.bf");
      TEST_ASSERT(system.hasEnvironment());
      TEST_ASSERT(!system.hasIterator());

      // The FilmEnvironment generates the mask and the wall fields
      TEST_ASSERT(!readCommandsThrows(system, "COMPUTE\nFINISH\n"));
      TEST_ASSERT(system.mask().hasData());
      TEST_ASSERT(system.h().hasData());
      TEST_ASSERT(system.c().hasData());
   }

   /*
   * WRITE_STRESS in a film system without an Iterator.
   */
   void testFilmNoIteratorStress()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testFilmNoIteratorStress.log");

      Rp::System<1,CPT> system;
      initSystem(system, "in/Command/film_param", "in/Command/film_w.bf");
      TEST_ASSERT(!system.hasIterator());

      std::string stressFile = "out/testFilmNoIteratorStress_stress";
      removeFile(stressFile);
      TEST_ASSERT(!readCommandsThrows(system,
                  "COMPUTE\nWRITE_STRESS " + stressFile + "\nFINISH\n"));
      TEST_ASSERT(system.mixture().hasStress());
      TEST_ASSERT(system.environment().hasStress());
      TEST_ASSERT(fileExists(stressFile));
   }

   /*
   * A last command not followed by a newline is executed.
   */
   void testLastCommandNoNewline()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLastCommandNoNewline.log");

      Rp::System<1,CPT> system;
      initSystem(system, "in/Command/lam_param", "in/Command/lam_w.bf");
      TEST_ASSERT(!system.c().hasData());

      // Script without FINISH that ends right after the command name
      TEST_ASSERT(!readCommandsThrows(system, "COMPUTE"));
      TEST_ASSERT(system.c().hasData());
   }

   /*
   * A command with an unreadable argument throws an Exception before
   * it is executed.
   */
   void testBadArgument()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testBadArgument.log");

      Rp::System<1,CPT> system;
      initSystem(system, "in/Command/lam_param", "in/Command/lam_w.bf");
      double length = system.domain().unitCell().parameter(0);

      // Unreadable cell parameter: the unit cell is not changed
      TEST_ASSERT(readCommandsThrows(system,
                  "SET_UNIT_CELL lamellar x\nFINISH\n"));
      TEST_ASSERT(eq(system.domain().unitCell().parameter(0), length));

      // Unreadable direction id: no propagator file is written
      std::string qFile = "out/testBadArgument_q.rf";
      removeFile(qFile);
      TEST_ASSERT(!readCommandsThrows(system, "COMPUTE\nFINISH\n"));
      TEST_ASSERT(readCommandsThrows(system,
                  "WRITE_Q_TAIL " + qFile + " 0 0 x\nFINISH\n"));
      TEST_ASSERT(!fileExists(qFile));
   }

   /*
   * WRITE_STRESS needs a solution for the current w fields and cell.
   */
   void testWriteStressNoSolution()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testWriteStressNoSolution.log");

      Rp::System<1,CPT> system;
      initSystem(system, "in/Command/lam_param", "in/Command/lam_w.bf");
      std::string root = "out/testWriteStressNoSolution_";
      for (int i = 1; i <= 4; ++i) {
         removeFile(root + std::to_string(i));
      }

      // Stress of the solution for the current fields and cell
      TEST_ASSERT(!readCommandsThrows(system,
                  "COMPUTE\nWRITE_STRESS " + root + "1\nFINISH\n"));
      TEST_ASSERT(fileExists(root + "1"));

      // After SET_UNIT_CELL, the propagators belong to the old cell
      TEST_ASSERT(!readCommandsThrows(system,
                  "SET_UNIT_CELL lamellar 1.7\nFINISH\n"));
      TEST_ASSERT(readCommandsThrows(system,
                  "WRITE_STRESS " + root + "2\nFINISH\n"));
      TEST_ASSERT(!fileExists(root + "2"));

      // A new solution can be written
      TEST_ASSERT(!readCommandsThrows(system,
                  "COMPUTE\nWRITE_STRESS " + root + "3\nFINISH\n"));
      TEST_ASSERT(fileExists(root + "3"));

      // After READ_W_BASIS, the solution belongs to the old fields
      TEST_ASSERT(!readCommandsThrows(system,
                  "READ_W_BASIS in/Command/lam_w.bf\nFINISH\n"));
      TEST_ASSERT(readCommandsThrows(system,
                  "WRITE_STRESS " + root + "4\nFINISH\n"));
      TEST_ASSERT(!fileExists(root + "4"));
   }

   /*
   * WRITE_BLOCK_C_RGRID and WRITE_Q* need a solution for the current
   * w fields and unit cell.
   */
   void testWriteBlockCAndQNoSolution()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testWriteBlockCAndQNoSolution.log");

      Rp::System<1,CPT> system;
      initSystem(system, "in/Command/lam_param", "in/Command/lam_w.bf");

      // Commands, and the (first) file that each one writes
      std::string root = "out/testWriteBlockCAndQNoSolution_";
      std::vector<std::string> commands, files;
      commands.push_back("WRITE_BLOCK_C_RGRID " + root + "blockC.rf");
      files.push_back(root + "blockC.rf");
      commands.push_back("WRITE_Q_SLICE " + root + "qSlice.rf 0 0 0 10");
      files.push_back(root + "qSlice.rf");
      commands.push_back("WRITE_Q_TAIL " + root + "qTail.rf 0 0 1");
      files.push_back(root + "qTail.rf");
      commands.push_back("WRITE_Q " + root + "q.rq 0 1 0");
      files.push_back(root + "q.rq");
      commands.push_back("WRITE_Q_ALL " + root + "qAll");
      files.push_back(root + "qAll_0_0_0.rq");
      int n = commands.size();

      // w fields have been read, but not yet solved for
      for (int i = 0; i < n; ++i) {
         removeFile(files[i]);
         TEST_ASSERT(readCommandsThrows(system,
                                        commands[i] + "\nFINISH\n"));
         TEST_ASSERT(!fileExists(files[i]));
      }

      // Solution for the current fields and cell
      TEST_ASSERT(!readCommandsThrows(system, "COMPUTE\nFINISH\n"));
      for (int i = 0; i < n; ++i) {
         TEST_ASSERT(!readCommandsThrows(system,
                                         commands[i] + "\nFINISH\n"));
         TEST_ASSERT(fileExists(files[i]));
      }

      // After SET_UNIT_CELL, the solution belongs to the old cell
      TEST_ASSERT(!readCommandsThrows(system,
                  "SET_UNIT_CELL lamellar 1.7\nFINISH\n"));
      for (int i = 0; i < n; ++i) {
         removeFile(files[i]);
         TEST_ASSERT(readCommandsThrows(system,
                                        commands[i] + "\nFINISH\n"));
         TEST_ASSERT(!fileExists(files[i]));
      }
   }

   /*
   * WRITE_PARAM of a bead model system writes a file that is read back
   * as the bead model.
   */
   void testWriteParamBead()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testWriteParamBead.log");

      std::string paramFile = "out/testWriteParamBead_param";
      {
         Rp::System<1,CPT> system;
         initSystem(system, "in/Command/bead_param");
         TEST_ASSERT(PolymerModel::isBead());
         TEST_ASSERT(!readCommandsThrows(system,
                     "WRITE_PARAM " + paramFile + "\nFINISH\n"));
      }

      // Read the written parameter file (as a thread model file,
      // it would lack the required ds)
      Rp::System<1,CPT> system;
      bool isRead = true;
      try {
         initSystem(system, paramFile);
      } catch (Exception& e) {
         isRead = false;
      }
      TEST_ASSERT(isRead);
      TEST_ASSERT(PolymerModel::isBead());
   }

   /*
   * COMPRESS with an LrCompressor after a change of the unit cell.
   */
   void testCompressCellChange()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testCompressCellChange.log");

      Rp::System<1,CPT> system;
      initSystem(system, "in/Command/fts_param");
      TEST_ASSERT(!readCommandsThrows(system,
                  "READ_W_RGRID in/Command/fts_w.rf\nCOMPRESS\n"
                  "SET_UNIT_CELL lamellar 10.0\nCOMPRESS\nFINISH\n"));
      TEST_ASSERT(eq(system.domain().unitCell().parameter(0), 10.0));

      // The second COMPRESS gives incompressible fields in the new cell
      TEST_ASSERT(system.c().hasData());
      double error = maxIncompressibilityError(system);
      if (verbose() > 0) {
         std::cout << "\nmax |sum of c - 1| = " << error;
      }
      TEST_ASSERT(error < 1.0E-6);
   }

   /*
   * Commands after a failed COMPRESS are not executed.
   */
   void testCompressFailure()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testCompressFailure.log");

      // The compressor of this system stops after 2 iterations
      Rp::System<1,CPT> system;
      initSystem(system, "in/Command/fts_param_fail");
      std::string wFile = "out/testCompressFailure_w.rf";
      removeFile(wFile);
      TEST_ASSERT(!readCommandsThrows(system,
                  "READ_W_RGRID in/Command/fts_w.rf\nCOMPRESS\n"
                  "WRITE_W_RGRID " + wFile + "\nFINISH\n"));
      TEST_ASSERT(system.simulator().compressor().mdeCounter() > 0);
      TEST_ASSERT(!fileExists(wFile));
   }

   // Utility functions

   /*
   * Read a parameter file and, if a file name is given, w fields in
   * basis format.
   */
   template <int D>
   void initSystem(Rp::System<D,CPT>& system,
                   std::string const & paramFileName,
                   std::string const & wFileName = "")
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile(paramFileName, in);
      system.readParam(in);
      in.close();

      if (!wFileName.empty()) {
         system.w().readBasis(wFileName);
      }
   }

   /*
   * Execute a command script given as a string, as from a command file.
   *
   * Returns true if a Util::Exception was thrown.
   */
   template <int D>
   bool readCommandsThrows(Rp::System<D,CPT>& system,
                           std::string const & commands)
   {
      std::istringstream in(commands);
      try {
         system.readCommands(in);
      } catch (Exception& e) {
         return true;
      }
      return false;
   }

   /*
   * Does a file exist? (name relative to the test directory)
   */
   bool fileExists(std::string const & filename)
   {
      std::ifstream file((filePrefix() + filename).c_str());
      return file.is_open();
   }

   /*
   * Remove an output file left by an earlier run.
   */
   void removeFile(std::string const & filename)
   {  std::remove((filePrefix() + filename).c_str()); }

   /*
   * Maximum deviation of the total concentration from 1 (NaN if any
   * total concentration is not a number).
   */
   template <int D>
   double maxIncompressibilityError(Rp::System<D,CPT> const & system)
   {
      int nMonomer = system.mixture().nMonomer();
      int meshSize = system.domain().mesh().size();
      double maxError = 0.0;
      for (int j = 0; j < meshSize; ++j) {
         double sum = 0.0;
         for (int i = 0; i < nMonomer; ++i) {
            sum += system.c().rgrid(i)[j];
         }
         if (std::isnan(sum)) return sum;
         double error = std::abs(sum - 1.0);
         if (error > maxError) maxError = error;
      }
      return maxError;
   }

};

TEST_BEGIN(CommandTest)
TEST_ADD(CommandTest, testAnalyzeNoSimulator)
TEST_ADD(CommandTest, testFilmNoIterator)
TEST_ADD(CommandTest, testFilmNoIteratorStress)
TEST_ADD(CommandTest, testLastCommandNoNewline)
TEST_ADD(CommandTest, testBadArgument)
TEST_ADD(CommandTest, testWriteStressNoSolution)
TEST_ADD(CommandTest, testWriteBlockCAndQNoSolution)
TEST_ADD(CommandTest, testWriteParamBead)
TEST_ADD(CommandTest, testCompressCellChange)
TEST_ADD(CommandTest, testCompressFailure)
TEST_END(CommandTest)

#endif
