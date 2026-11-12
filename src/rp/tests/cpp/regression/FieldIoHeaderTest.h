#ifndef RPC_FIELD_IO_HEADER_TEST_H
#define RPC_FIELD_IO_HEADER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Propagator.h>
#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>
#include <rp/field/WFields.h>
#include <rp/field/Mask.h>

#include <pscf/field/cpp/RField.h>
#include <pscf/crystal/UnitCell.h>

#include <util/containers/FSArray.h>
#include <util/misc/Exception.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests for the headers of r-grid files of single fields:
* propagator slices and tails (WRITE_Q_SLICE, WRITE_Q_TAIL) and masks
* (WRITE_MASK_RGRID), for fields with and without space group symmetry.
*/
class FieldIoHeaderTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Read a parameter file to initialize a System.
   */
   template <int D>
   void setupSystem(Rp::System<D,CPT>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

   /*
   * Return the first word of a text file.
   */
   std::string firstWord(std::string filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      std::string word;
      in >> word;
      in.close();
      return word;
   }

   /*
   * Does a text file contain a string?
   */
   bool fileContains(std::string filename, std::string const & text)
   {
      std::ifstream in;
      openInputFile(filename, in);
      std::string line;
      bool found = false;
      while (std::getline(in, line)) {
         if (line.find(text) != std::string::npos) {
            found = true;
         }
      }
      in.close();
      return found;
   }

   /*
   * Maximum absolute difference between two r-grid fields.
   */
   template <int D>
   double maxDiff(RField<D,CPT> const & a, RField<D,CPT> const & b)
   {
      UTIL_CHECK(a.capacity() == b.capacity());
      double diff;
      double max = 0.0;
      for (int i = 0; i < a.capacity(); ++i) {
         diff = std::abs(a[i] - b[i]);
         if (diff > max) max = diff;
      }
      return max;
   }

   /*
   * Check an r-grid file of a single field written by the system: it
   * must start with a full header, name the space group iff hasGroup,
   * and read back to the expected field.
   */
   template <int D>
   void checkRGridFile(Rp::System<D,CPT>& system,
                       std::string filename,
                       RField<D,CPT> const & expected,
                       bool hasGroup)
   {
      TEST_ASSERT(firstWord(filename) == "format");
      TEST_ASSERT(fileContains(filename, "mesh"));
      TEST_ASSERT(fileContains(filename, "group_name") == hasGroup);

      RField<D,CPT> field;
      field.allocate(system.domain().mesh().dimensions());
      UnitCell<D> unitCell;
      bool success = true;
      try {
         system.domain().fieldIo().readFieldRGrid(filename, field,
                                                  unitCell);
      } catch (Exception& e) {
         success = false;
      }
      TEST_ASSERT(success);
      TEST_ASSERT(maxDiff(field, expected) < 1.0E-10);
      TEST_ASSERT(eq(unitCell.parameter(0),
                     system.domain().unitCell().parameter(0)));
   }

   /*
   * Write a slice and the tail of propagator (block 1, direction 0) of
   * polymer 0 and check both files.
   */
   template <int D>
   void checkQSliceTail(Rp::System<D,CPT>& system, std::string root,
                        bool hasGroup)
   {
      auto const & propagator = system.mixture().polymer(0).propagator(1, 0);
      system.mixture().writeQSlice(root + "_qslice.rf", 0, 1, 0, 3);
      checkRGridFile(system, root + "_qslice.rf", propagator.q(3),
                     hasGroup);
      system.mixture().writeQTail(root + "_qtail.rf", 0, 1, 0);
      checkRGridFile(system, root + "_qtail.rf", propagator.tail(),
                     hasGroup);
   }

   /*
   * Propagator slice and tail of a branched (star) polymer in a system
   * without a space group.
   */
   void testQSliceTailNoGroup_star()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FieldIoHeaderQStar.log");

      Rp::System<1,CPT> system;
      setupSystem(system, "in/FieldIo/hdr_star_param");
      system.w().readRGrid("in/FieldIo/hdr_w.rf");
      system.compute();
      checkQSliceTail(system, "out/FieldIoHeader_star", false);
   }

   /*
   * Propagator slice and tail of a linear polymer in a system without
   * a space group.
   */
   void testQSliceTailNoGroup_linear()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FieldIoHeaderQLinear.log");

      Rp::System<1,CPT> system;
      setupSystem(system, "in/FieldIo/hdr_linear_param");
      system.w().readRGrid("in/FieldIo/hdr_w.rf");
      system.compute();
      checkQSliceTail(system, "out/FieldIoHeader_linear", false);
   }

   /*
   * Propagator slice and tail in a system with a space group (P_-1):
   * symmetric fields (w read in basis format) name the group, while
   * fields computed from w read in r-grid format do not.
   */
   void testQSliceTailGroup()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FieldIoHeaderQGroup.log");

      Rp::System<1,CPT> system;
      setupSystem(system, "in/FieldIo/hdr_sym_param");
      system.w().readBasis("in/FieldIo/hdr_w.bf");
      system.compute();
      checkQSliceTail(system, "out/FieldIoHeader_sym", true);

      system.w().writeRGrid("out/FieldIoHeader_sym_w.rf");
      system.w().readRGrid("out/FieldIoHeader_sym_w.rf");
      system.compute();
      checkQSliceTail(system, "out/FieldIoHeader_asym", false);
   }

   /*
   * Masks in r-grid format: a mask without symmetry (no space group)
   * and a symmetric mask (group P_-1, read in basis format).
   */
   void testMaskRGrid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FieldIoHeaderMask.log");
      FSArray<double, 6> parameters;
      parameters.append(2.0);

      Rp::System<1,CPT> system;
      setupSystem(system, "in/FieldIo/hdr_linear_param");
      system.setUnitCell(parameters);
      system.mask().readRGrid("in/FieldIo/hdr_mask.rf");
      system.mask().writeRGrid("out/FieldIoHeader_mask.rf");
      checkRGridFile(system, "out/FieldIoHeader_mask.rf",
                     system.mask().rgrid(), false);

      Rp::System<1,CPT> systemSym;
      setupSystem(systemSym, "in/FieldIo/hdr_sym_param");
      systemSym.setUnitCell(parameters);
      systemSym.mask().readBasis("in/FieldIo/hdr_mask.bf");
      systemSym.mask().writeRGrid("out/FieldIoHeader_mask_sym.rf");
      checkRGridFile(systemSym, "out/FieldIoHeader_mask_sym.rf",
                     systemSym.mask().rgrid(), true);
   }

};

TEST_BEGIN(FieldIoHeaderTest)
TEST_ADD(FieldIoHeaderTest, testQSliceTailNoGroup_star)
TEST_ADD(FieldIoHeaderTest, testQSliceTailNoGroup_linear)
TEST_ADD(FieldIoHeaderTest, testQSliceTailGroup)
TEST_ADD(FieldIoHeaderTest, testMaskRGrid)
TEST_END(FieldIoHeaderTest)

#endif
