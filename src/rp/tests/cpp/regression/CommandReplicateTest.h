#ifndef RPC_COMMAND_REPLICATE_TEST_H
#define RPC_COMMAND_REPLICATE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>

#include <util/misc/Exception.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>

using namespace Util;
using namespace Pscf;

/*
* Regression tests of the REPLICATE_UNIT_CELL command: zero or negative
* numbers of replicas must be rejected before any file is opened. They
* were rejected only after the output file had been opened, which was
* left behind empty.
*/
class CommandReplicateTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   void testReplicate1D()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testReplicate1D.log");
      std::string const prefix = "out/testReplicate1D_";

      Rp::System<1,CPT> system;
      initSystem(system, "in/Command/lam_param");
      std::string const inFile = prefix + "w.rf";
      TEST_ASSERT(!readCommandsThrows(system,
                  "BASIS_TO_RGRID in/Command/lam_w.bf " + inFile
                  + "\nFINISH\n"));

      TEST_ASSERT(replicateThrows(system, inFile, prefix + "r1", "0"));
      TEST_ASSERT(replicateThrows(system, inFile, prefix + "r2", "-2"));

      // Valid input
      std::string const outFile = prefix + "ok.rf";
      removeFile(outFile);
      TEST_ASSERT(!readCommandsThrows(system,
                  "REPLICATE_UNIT_CELL " + inFile + " " + outFile
                  + " 3\nFINISH\n"));
      TEST_ASSERT(fileExists(outFile));
      Rp::System<1,CPT> system3;
      initSystem(system3, param("96", "lamellar"));
      TEST_ASSERT(!readCommandsThrows(system3,
                  "READ_W_RGRID " + outFile + "\nFINISH\n"));
      TEST_ASSERT(system3.w().hasData());
   }

   void testReplicate2D()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testReplicate2D.log");
      std::string const prefix = "out/testReplicate2D_";
      std::string const inFile = "../field/in/w_hex.rf";

      Rp::System<2,CPT> system;
      initSystem(system, param("32 32", "hexagonal", "p_6_m_m"));
      TEST_ASSERT(replicateThrows(system, inFile, prefix + "r1", "0 0"));
      TEST_ASSERT(replicateThrows(system, inFile, prefix + "r2", "2 0"));
      TEST_ASSERT(replicateThrows(system, inFile, prefix + "r3", "0 2"));
      TEST_ASSERT(replicateThrows(system, inFile, prefix + "r4",
                                  "-1 -1"));

      // Valid input
      std::string const outFile = prefix + "ok.rf";
      removeFile(outFile);
      TEST_ASSERT(!readCommandsThrows(system,
                  "REPLICATE_UNIT_CELL " + inFile + " " + outFile
                  + " 2 2\nFINISH\n"));
      TEST_ASSERT(fileExists(outFile));
   }

   // Utility functions

   /*
   * Does REPLICATE_UNIT_CELL with these numbers of replicas throw an
   * Exception without writing outFile?
   */
   template <int D>
   bool replicateThrows(Rp::System<D,CPT>& system,
                        std::string const & inFile,
                        std::string const & outFile,
                        std::string const & replicas)
   {
      removeFile(outFile);
      bool thrown = readCommandsThrows(system,
                       "REPLICATE_UNIT_CELL " + inFile + " " + outFile
                       + " " + replicas + "\nFINISH\n");
      if (verbose() > 0) {
         std::cout << "\nreplicas " << replicas << ": thrown = "
                   << thrown << ", file = " << fileExists(outFile);
      }
      return thrown && !fileExists(outFile);
   }

   /*
   * Parameter file of a diblock melt, with an optional space group.
   */
   std::string param(std::string const & mesh,
                     std::string const & lattice,
                     std::string const & group = "")
   {
      std::string groupLine = "";
      if (!group.empty()) {
         groupLine = "     groupName   " + group + "\n";
      }
      return std::string("System{\n")
             + "  Mixture{\n"
               "    nMonomer  2\n"
               "    monomers[\n"
               "              1.0\n"
               "              1.0\n"
               "    ]\n"
               "    nPolymer  1\n"
               "    Polymer{\n"
               "       type    linear\n"
               "       nBlock  2\n"
               "       blocks[\n"
               "               0  0.3\n"
               "               1  0.7\n"
               "       ]\n"
               "       phi     1.0\n"
               "    }\n"
               "    ds   0.01\n"
               "  }\n"
               "  Interaction{\n"
               "    chi(\n"
               "         1   0   20.0\n"
               "    )\n"
               "  }\n"
               "  Domain{\n"
               "     mesh        " + mesh + "\n"
               "     lattice     " + lattice + "\n"
             + groupLine
             + "  }\n"
               "}\n";
   }

   /*
   * Read a parameter file, given as a file name or as a string that
   * starts with "System{".
   */
   template <int D>
   void initSystem(Rp::System<D,CPT>& system, std::string const & param)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      if (param.compare(0, 7, "System{") == 0) {
         std::istringstream in(param);
         system.readParam(in);
      } else {
         std::ifstream in;
         openInputFile(param, in);
         system.readParam(in);
         in.close();
      }
   }

   /*
   * Execute a command script given as a string. Returns true if a
   * Util::Exception was thrown.
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

};

TEST_BEGIN(CommandReplicateTest)
TEST_ADD(CommandReplicateTest, testReplicate1D)
TEST_ADD(CommandReplicateTest, testReplicate2D)
TEST_END(CommandReplicateTest)

#endif
