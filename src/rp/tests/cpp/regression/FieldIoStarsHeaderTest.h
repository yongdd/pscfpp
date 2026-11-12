#ifndef RPC_FIELD_IO_STARS_HEADER_TEST_H
#define RPC_FIELD_IO_STARS_HEADER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/field/Domain.h>

#include <pscf/crystal/UnitCell.h>

#include <util/misc/Exception.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests of the headers of WRITE_STARS and WRITE_WAVES.
*/
class FieldIoStarsHeaderTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * A run that only converts field files (as in the hex 'convert'
   * example) writes the unit cell from which the basis was built, with
   * the same listing as a run that reads w first. The header had
   * cell_param 0, since the system unit cell was never set. Once the
   * system unit cell is set, the header gives it, as before.
   */
   void testConvertOnly()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testStarsHeaderConvertOnly.log");
      std::string const prefix = "out/testStarsHeaderConvertOnly_";
      std::string const wFile = "../field/in/w_hex.bf";

      // Conversions only
      Rp::System<2,CPT> system;
      initSystem(system, hexParam());
      TEST_ASSERT(!readCommandsThrows(system,
                  "BASIS_TO_RGRID " + wFile + " " + prefix + "w.rf\n"
                  "RGRID_TO_BASIS " + prefix + "w.rf " + prefix + "w.bf\n"
                  "WRITE_STARS " + prefix + "stars\n"
                  "WRITE_WAVES " + prefix + "waves\nFINISH\n"));
      TEST_ASSERT(!system.domain().unitCell().isInitialized());
      double a = 1.6651597308;
      TEST_ASSERT(std::abs(cellParam(prefix + "stars") - a) < 1.0E-10);
      TEST_ASSERT(std::abs(cellParam(prefix + "waves") - a) < 1.0E-10);

      // Reference: READ_W_BASIS first
      Rp::System<2,CPT> ref;
      initSystem(ref, hexParam());
      TEST_ASSERT(!readCommandsThrows(ref,
                  "READ_W_BASIS " + wFile + "\n"
                  "WRITE_STARS " + prefix + "stars_ref\n"
                  "WRITE_WAVES " + prefix + "waves_ref\nFINISH\n"));
      TEST_ASSERT(sameFiles(prefix + "stars", prefix + "stars_ref"));
      TEST_ASSERT(sameFiles(prefix + "waves", prefix + "waves_ref"));

      // Once the system unit cell is set, the header gives it
      TEST_ASSERT(!readCommandsThrows(ref,
                  "SET_UNIT_CELL hexagonal 1.7\n"
                  "WRITE_STARS " + prefix + "stars_set\n"
                  "WRITE_WAVES " + prefix + "waves_set\nFINISH\n"));
      TEST_ASSERT(std::abs(cellParam(prefix + "stars_set") - 1.7) < 1.0E-10);
      TEST_ASSERT(std::abs(cellParam(prefix + "waves_set") - 1.7) < 1.0E-10);
   }

   // Utility functions

   /*
   * Parameter file of a hexagonal diblock melt, without unit cell.
   */
   std::string hexParam()
   {
      return "System{\n"
             "  Mixture{\n"
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
             "     mesh        32    32\n"
             "     lattice     hexagonal\n"
             "     groupName   p_6_m_m\n"
             "  }\n"
             "}\n";
   }

   /*
   * Read a parameter file given as a string.
   */
   template <int D>
   void initSystem(Rp::System<D,CPT>& system, std::string const & text)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream in(text);
      system.readParam(in);
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
   * Read the first cell parameter from the header of a file.
   */
   double cellParam(std::string const & filename)
   {
      std::ifstream in((filePrefix() + filename).c_str());
      UTIL_CHECK(in.is_open());
      std::string word;
      while (in >> word) {
         if (word == "cell_param") {
            double value;
            in >> value;
            UTIL_CHECK(!in.fail());
            return value;
         }
      }
      UTIL_THROW("No cell_param in file header");
      return 0.0;
   }

   /*
   * Do two files have the same content?
   */
   bool sameFiles(std::string const & name1, std::string const & name2)
   {
      std::ifstream in1((filePrefix() + name1).c_str());
      std::ifstream in2((filePrefix() + name2).c_str());
      UTIL_CHECK(in1.is_open());
      UTIL_CHECK(in2.is_open());
      std::stringstream s1, s2;
      s1 << in1.rdbuf();
      s2 << in2.rdbuf();
      return (s1.str() == s2.str()) && (s1.str().size() > 0);
   }

};

TEST_BEGIN(FieldIoStarsHeaderTest)
TEST_ADD(FieldIoStarsHeaderTest, testConvertOnly)
TEST_END(FieldIoStarsHeaderTest)

#endif
