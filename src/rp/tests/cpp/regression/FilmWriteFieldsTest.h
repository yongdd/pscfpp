#ifndef RPC_FILM_WRITE_FIELDS_TEST_H
#define RPC_FILM_WRITE_FIELDS_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>
#include <rp/field/WFields.h>
#include <rp/field/Mask.h>
#include <rp/environment/Environment.h>

#include <pscf/crystal/UnitCell.h>

#include <util/containers/DArray.h>
#include <util/misc/Exception.h>
#include <util/tests/LogFileUnitTest.h>

#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests for WRITE_MASK_* and WRITE_H_* after a change of the
* unit cell of a thin film.
*
* The film is 1D (group P_-1, walls with chi = 3 for monomer 0), with
* the w fields of in/Film/sweep1D_w.bf (L = 1.5). The written mask and
* wall fields must be those of the cell in the file header.
*/
class FilmWriteFieldsTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Read the parameter file and the w fields (L = 1.5).
   */
   void initSystem(Rp::System<1,CPT>& system)
   {
      std::ostringstream param;
      param << "System{\n"
            << "  Mixture{\n"
            << "    nMonomer  2\n"
            << "    monomers[  1.0  1.0  ]\n"
            << "    nPolymer  1\n"
            << "    Polymer{\n"
            << "      type    linear\n"
            << "      nBlock  2\n"
            << "      blocks[\n"
            << "              0  0.40\n"
            << "              1  0.60\n"
            << "      ]\n"
            << "      phi     1.0\n"
            << "    }\n"
            << "    ds   0.01\n"
            << "  }\n"
            << "  Interaction{\n"
            << "    chi(\n"
            << "         1   0   8.0\n"
            << "    )\n"
            << "  }\n"
            << "  Domain{\n"
            << "    mesh           32\n"
            << "    lattice        lamellar\n"
            << "    groupName      P_-1\n"
            << "  }\n"
            << "  FilmEnvironment{\n"
            << "    normalVecId        0\n"
            << "    interfaceThickness 0.4\n"
            << "    excludedThickness  0.5\n"
            << "    chiBottom[  3.0  0.0  ]\n"
            << "    chiTop[  3.0  0.0  ]\n"
            << "  }\n"
            << "  AmIteratorBasis{\n"
            << "    epsilon      1.0e-9\n"
            << "    maxItr       100\n"
            << "    maxHist      25\n"
            << "    isFlexible   0\n"
            << "  }\n"
            << "}\n";

      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream in(param.str());
      system.readParam(in);
      TEST_ASSERT(system.hasEnvironment());
      system.w().readBasis("in/Film/sweep1D_w.bf");
   }

   /*
   * Execute a command script given as a string.
   *
   * Returns true if a Util::Exception was thrown.
   */
   bool readCommandsThrows(Rp::System<1,CPT>& system,
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

   // Maximum difference of two arrays of the same size
   double maxDiff(DArray<double> const & a, DArray<double> const & b)
   {
      TEST_ASSERT(a.capacity() == b.capacity());
      double diff = 0.0;
      for (int i = 0; i < a.capacity(); ++i) {
         double d = std::abs(a[i] - b[i]);
         if (d > diff) diff = d;
      }
      return diff;
   }

   /*
   * After SET_UNIT_CELL (L = 2.0 or 1.5), each of WRITE_MASK_BASIS,
   * WRITE_MASK_RGRID, WRITE_H_BASIS and WRITE_H_RGRID writes the mask
   * and wall fields generated for the new cell.
   */
   void testWriteAfterSetUnitCell()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmWriteFieldsSetUnitCell.log");
      std::string root = "out/FilmWriteFields_";

      // Reference mask and h for L = 1.5 (k = 0) and L = 2.0 (k = 1)
      double length[2] = {1.5, 2.0};
      DArray<double> mask[2], h0[2];
      for (int k = 0; k < 2; ++k) {
         Rp::System<1,CPT> ref;
         initSystem(ref);
         std::ostringstream cmd;
         cmd << "SET_UNIT_CELL lamellar " << length[k] << "\n"
             << "COMPUTE\nFINISH\n";
         TEST_ASSERT(!readCommandsThrows(ref, cmd.str()));
         mask[k] = ref.mask().basis();
         h0[k] = ref.h().basis(0);
      }
      TEST_ASSERT(std::abs(mask[1][0] - mask[0][0]) > 1.0E-2);

      Rp::System<1,CPT> system;
      initSystem(system);
      TEST_ASSERT(!readCommandsThrows(system, "COMPUTE\nFINISH\n"));
      TEST_ASSERT(maxDiff(system.mask().basis(), mask[0]) < 1.0E-12);

      std::string commands[4] = {"WRITE_MASK_BASIS", "WRITE_MASK_RGRID",
                                 "WRITE_H_BASIS", "WRITE_H_RGRID"};
      std::string files[4] = {"mask.bf", "mask.rf", "h.bf", "h.rf"};
      for (int i = 0; i < 4; ++i) {

         // Change of the cell, then write
         int k = (i % 2 == 0) ? 1 : 0;
         std::ostringstream cmd;
         cmd << "SET_UNIT_CELL lamellar " << length[k] << "\n"
             << commands[i] << " " << root << files[i] << "\n"
             << "FINISH\n";
         TEST_ASSERT(!readCommandsThrows(system, cmd.str()));

         // Fields stored in the system (and written)
         TEST_ASSERT(maxDiff(system.mask().basis(), mask[k]) < 1.0E-12);
         TEST_ASSERT(maxDiff(system.h().basis(0), h0[k]) < 1.0E-12);
      }

      // Read back the files in basis format: header cell and data
      {
         DArray<double> field;
         UnitCell<1> cell;
         system.domain().fieldIo().readFieldBasis(root + files[0],
                                                  field, cell);
         TEST_ASSERT(eq(cell.parameter(0), length[1]));
         TEST_ASSERT(maxDiff(field, mask[1]) < 1.0E-9);
      }
      {
         DArray< DArray<double> > fields;
         UnitCell<1> cell;
         system.domain().fieldIo().readFieldsBasis(root + files[2],
                                                   fields, cell);
         TEST_ASSERT(eq(cell.parameter(0), length[1]));
         TEST_ASSERT(maxDiff(fields[0], h0[1]) < 1.0E-9);
      }
   }

};

TEST_BEGIN(FilmWriteFieldsTest)
TEST_ADD(FilmWriteFieldsTest, testWriteAfterSetUnitCell)
TEST_END(FilmWriteFieldsTest)

#endif
