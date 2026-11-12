#ifndef RPG_FILM_SPACE_GROUP_TEST_H
#define RPG_FILM_SPACE_GROUP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/Mask.h>
#include <rp/scft/ScftThermo.h>
#include <rp/environment/Environment.h>

#include <util/containers/DArray.h>
#include <util/misc/Exception.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests for the space group of thin-film systems.
*
* A film needs no space group, and the generated mask and wall fields
* are symmetric whenever the system has a space group, with or without
* an Iterator. All systems are 1D films (mesh 32) with the w fields of
* in/Film/sweep1D_w.bf.
*/
class FilmSpaceGroupTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Parameter file of a 1D film.
   *
   * \param group  group name ("" for none)
   * \param chiBottom  chi of the bottom wall with monomers 0 and 1
   * \param chiTop  chi of the top wall with monomers 0 and 1
   * \param iterator  "AmIteratorBasis", "AmIteratorGrid" or "" (none)
   */
   std::string filmParam(std::string const & group,
                         std::string const & chiBottom,
                         std::string const & chiTop,
                         std::string const & iterator)
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
            << "    lattice        lamellar\n";
      if (!group.empty()) {
         param << "    groupName      " << group << "\n";
      }
      param << "  }\n"
            << "  FilmEnvironment{\n"
            << "    normalVecId        0\n"
            << "    interfaceThickness 0.4\n"
            << "    excludedThickness  0.5\n"
            << "    chiBottom[  " << chiBottom << "  ]\n"
            << "    chiTop[  " << chiTop << "  ]\n"
            << "  }\n";
      if (!iterator.empty()) {
         param << "  " << iterator << "{\n"
               << "    epsilon      1.0e-9\n"
               << "    maxItr       1000\n"
               << "    maxHist      25\n"
               << "    isFlexible   0\n"
               << "  }\n";
      }
      param << "}\n";
      return param.str();
   }

   /*
   * Read a parameter file given as a string.
   */
   void initSystem(Rp::System<1,CUT>& system, std::string const & param)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream in(param);
      system.readParam(in);
      TEST_ASSERT(system.hasEnvironment());
   }

   /*
   * Execute a command script given as a string.
   *
   * Returns true if a Util::Exception was thrown.
   */
   bool readCommandsThrows(Rp::System<1,CUT>& system,
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

   // Helmholtz free energy of the current solution
   double fHelmholtz(Rp::System<1,CUT>& system)
   {
      system.scft().compute();
      return system.scft().fHelmholtz();
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

   // Does a file (name relative to the test directory) contain a string?
   bool fileContains(std::string const & filename, std::string const & s)
   {
      std::ifstream file((filePrefix() + filename).c_str());
      std::stringstream buffer;
      buffer << file.rdbuf();
      return (buffer.str().find(s) != std::string::npos);
   }

   /*
   * A film without a space group, with walls that are chemically
   * identical (checked by the mask generator) or not (also checked by
   * the wall field generator), gives the same solution as with a group.
   */
   void testNoSpaceGroup()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmSpaceGroupNone.log");
      std::string wFile = "out/FilmSpaceGroupNone_w.rf";
      std::string wFileNG = "out/FilmSpaceGroupNone_wNG.rf";

      // Reference with group P_-1 and identical walls
      Rp::System<1,CUT> ref;
      initSystem(ref, filmParam("P_-1", "3.0  0.0", "3.0  0.0", ""));
      TEST_ASSERT(!readCommandsThrows(ref,
                  "READ_W_BASIS in/Film/sweep1D_w.bf\n"
                  "WRITE_W_RGRID " + wFile + "\nCOMPUTE\nFINISH\n"));
      double fRef = fHelmholtz(ref);

      // No group, identical walls
      {
         Rp::System<1,CUT> system;
         initSystem(system, filmParam("", "3.0  0.0", "3.0  0.0",
                                      "AmIteratorGrid"));
         TEST_ASSERT(!system.domain().hasGroup());
         TEST_ASSERT(!readCommandsThrows(system,
                     "READ_W_RGRID " + wFile + "\nCOMPUTE\n"
                     "WRITE_W_RGRID " + wFileNG + "\nFINISH\n"));
         TEST_ASSERT(system.mask().hasData());
         TEST_ASSERT(!system.mask().isSymmetric());
         TEST_ASSERT(system.h().hasData());
         TEST_ASSERT(std::abs(fHelmholtz(system) - fRef) < 1.0E-10);

         // ITERATE converges
         TEST_ASSERT(system.iterate() == 0);
      }

      // Reference with group P_1 and different walls
      Rp::System<1,CUT> ref1;
      initSystem(ref1, filmParam("P_1", "5.0  0.0", "2.0  10.0", ""));
      TEST_ASSERT(!readCommandsThrows(ref1,
                  "READ_W_RGRID " + wFileNG + "\nCOMPUTE\nFINISH\n"));
      double fRef1 = fHelmholtz(ref1);

      // No group, different walls
      {
         Rp::System<1,CUT> system;
         initSystem(system, filmParam("", "5.0  0.0", "2.0  10.0",
                                      "AmIteratorGrid"));
         TEST_ASSERT(!readCommandsThrows(system,
                     "READ_W_RGRID " + wFileNG + "\nCOMPUTE\nFINISH\n"));
         TEST_ASSERT(system.h().hasData());
         TEST_ASSERT(!system.h().isSymmetric());
         TEST_ASSERT(std::abs(fHelmholtz(system) - fRef1) < 1.0E-10);
      }
   }

   /*
   * In a film system with a space group, the mask and wall fields are
   * symmetric without an Iterator and with AmIteratorGrid, as with
   * AmIteratorBasis, so WRITE_MASK_BASIS and WRITE_H_BASIS work.
   */
   void testSymmetricFields()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmSpaceGroupSymmetric.log");
      std::string root = "out/FilmSpaceGroupSymmetric_";
      std::string commands = "READ_W_BASIS in/Film/sweep1D_w.bf\nCOMPUTE\n";

      // Reference with AmIteratorBasis
      Rp::System<1,CUT> ref;
      initSystem(ref, filmParam("P_-1", "3.0  0.0", "3.0  0.0",
                                "AmIteratorBasis"));
      TEST_ASSERT(!readCommandsThrows(ref, commands + "FINISH\n"));
      TEST_ASSERT(ref.mask().isSymmetric());
      TEST_ASSERT(ref.h().isSymmetric());

      for (int k = 0; k < 2; ++k) {
         std::string iterator = (k == 0) ? "" : "AmIteratorGrid";
         std::string prefix = root + std::to_string(k);
         Rp::System<1,CUT> system;
         initSystem(system, filmParam("P_-1", "3.0  0.0", "3.0  0.0",
                                      iterator));
         TEST_ASSERT(system.hasIterator() == (k == 1));
         TEST_ASSERT(!readCommandsThrows(system, commands
                     + "WRITE_MASK_BASIS " + prefix + "_mask.bf\n"
                     + "WRITE_H_BASIS " + prefix + "_h.bf\n"
                     + "WRITE_MASK_RGRID " + prefix + "_mask.rf\n"
                     + "FINISH\n"));
         TEST_ASSERT(system.mask().isSymmetric());
         TEST_ASSERT(system.h().isSymmetric());
         TEST_ASSERT(maxDiff(system.mask().basis(), ref.mask().basis())
                     < 1.0E-12);
         for (int i = 0; i < 2; ++i) {
            TEST_ASSERT(maxDiff(system.h().basis(i), ref.h().basis(i))
                        < 1.0E-12);
         }
         TEST_ASSERT(fileContains(prefix + "_mask.rf", "P_-1"));
      }
   }

};

TEST_BEGIN(FilmSpaceGroupTest)
TEST_ADD(FilmSpaceGroupTest, testNoSpaceGroup)
TEST_ADD(FilmSpaceGroupTest, testSymmetricFields)
TEST_END(FilmSpaceGroupTest)

#endif
