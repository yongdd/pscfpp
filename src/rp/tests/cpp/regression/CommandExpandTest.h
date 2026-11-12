#ifndef RPC_COMMAND_EXPAND_TEST_H
#define RPC_COMMAND_EXPAND_TEST_H

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
* Regression tests of the EXPAND_RGRID_DIMENSION command.
*/
class CommandExpandTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Expanding a 1D field: zero or negative numbers of grid points in
   * the added dimensions, and invalid dimensions, are rejected before
   * any file is written (a broken file was written with exit code 0).
   */
   void testExpand1D()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testExpand1D.log");
      std::string const prefix = "out/testExpand1D_";

      Rp::System<1,CPT> system;
      initSystem(system, "in/Command/lam_param");
      std::string const inFile = prefix + "w.rf";
      TEST_ASSERT(!readCommandsThrows(system,
                  "BASIS_TO_RGRID in/Command/lam_w.bf " + inFile
                  + "\nFINISH\n"));

      TEST_ASSERT(expandThrows(system, inFile, prefix + "g1", "3 0 5"));
      TEST_ASSERT(expandThrows(system, inFile, prefix + "g2", "3 -3 5"));
      TEST_ASSERT(expandThrows(system, inFile, prefix + "g3", "3 4 0"));
      TEST_ASSERT(expandThrows(system, inFile, prefix + "g4", "2 0"));
      TEST_ASSERT(expandThrows(system, inFile, prefix + "g5", "2 -1"));
      TEST_ASSERT(expandThrows(system, inFile, prefix + "g6", "1 4"));

      // Valid input
      std::string const outFile = prefix + "ok.rf";
      removeFile(outFile);
      TEST_ASSERT(!readCommandsThrows(system,
                  "EXPAND_RGRID_DIMENSION " + inFile + " " + outFile
                  + " 3 4 5\nFINISH\n"));
      TEST_ASSERT(fileExists(outFile));
      Rp::System<3,CPT> system3;
      initSystem(system3, param3D("32 4 5", "orthorhombic"));
      TEST_ASSERT(!readCommandsThrows(system3,
                  "READ_W_RGRID " + outFile + "\nFINISH\n"));
      TEST_ASSERT(system3.w().hasData());
   }

   /*
   * Expanding a 2D field: zero or negative numbers of grid points are
   * rejected before any file is written.
   */
   void testExpand2D()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testExpand2D.log");
      std::string const prefix = "out/testExpand2D_";
      std::string const inFile = "../field/in/w_hex.rf";

      Rp::System<2,CPT> system;
      initSystem(system, hexParam());
      TEST_ASSERT(expandThrows(system, inFile, prefix + "g1", "3 0"));
      TEST_ASSERT(expandThrows(system, inFile, prefix + "g2", "3 -4"));
      TEST_ASSERT(expandThrows(system, inFile, prefix + "g3", "4 4"));

      // Valid input
      std::string const outFile = prefix + "ok.rf";
      removeFile(outFile);
      TEST_ASSERT(!readCommandsThrows(system,
                  "EXPAND_RGRID_DIMENSION " + inFile + " " + outFile
                  + " 3 4\nFINISH\n"));
      TEST_ASSERT(fileExists(outFile));
   }

   // Utility functions

   /*
   * Does EXPAND_RGRID_DIMENSION with these arguments (d and the numbers
   * of grid points) throw an Exception without writing outFile?
   */
   template <int D>
   bool expandThrows(Rp::System<D,CPT>& system,
                     std::string const & inFile,
                     std::string const & outFile,
                     std::string const & args)
   {
      removeFile(outFile);
      bool thrown = readCommandsThrows(system,
                       "EXPAND_RGRID_DIMENSION " + inFile + " " + outFile
                       + " " + args + "\nFINISH\n");
      if (verbose() > 0) {
         std::cout << "\nargs " << args << ": thrown = " << thrown
                   << ", file = " << fileExists(outFile);
      }
      return thrown && !fileExists(outFile);
   }

   /*
   * Parameter file of a hexagonal diblock melt.
   */
   std::string hexParam()
   {
      return std::string("System{\n") + mixture()
             + "  Domain{\n"
               "     mesh        32    32\n"
               "     lattice     hexagonal\n"
               "     groupName   p_6_m_m\n"
               "  }\n"
               "}\n";
   }

   /*
   * Parameter file of a 3D system without space group.
   */
   std::string param3D(std::string const & mesh,
                       std::string const & lattice)
   {
      return std::string("System{\n") + mixture()
             + "  Domain{\n"
               "     mesh        " + mesh + "\n"
               "     lattice     " + lattice + "\n"
               "  }\n"
               "}\n";
   }

   /*
   * Mixture and Interaction blocks of a diblock melt.
   */
   std::string mixture()
   {
      return "  Mixture{\n"
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
             "  }\n";
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

TEST_BEGIN(CommandExpandTest)
TEST_ADD(CommandExpandTest, testExpand1D)
TEST_ADD(CommandExpandTest, testExpand2D)
TEST_END(CommandExpandTest)

#endif
