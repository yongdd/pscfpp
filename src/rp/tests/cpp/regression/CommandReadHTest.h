#ifndef RPC_COMMAND_READ_H_TEST_H
#define RPC_COMMAND_READ_H_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <pscf/field/cpp/RField.h>

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
* Regression tests of the commands that read external fields, in a
* system without an Environment (the h basis arrays are allocated on
* first use).
*/
class CommandReadHTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * READ_H_BASIS after READ_H_RGRID gives the fields of READ_H_BASIS
   * alone (it aborted in WFields::allocateBasis: '!hasData_').
   */
   void testReadHBasisAfterRGrid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testReadHBasisAfterRGrid.log");

      // Reference: READ_H_BASIS only
      Rp::System<1,CPT> ref;
      initSystem(ref);
      TEST_ASSERT(!readCommandsThrows(ref,
                  "READ_W_BASIS in/Command/lam_w.bf\n"
                  "READ_H_BASIS in/Command/lam_h.bf\nFINISH\n"));
      TEST_ASSERT(ref.h().hasData());
      TEST_ASSERT(ref.h().isSymmetric());

      // READ_H_RGRID, then READ_H_BASIS
      Rp::System<1,CPT> system;
      initSystem(system);
      TEST_ASSERT(!readCommandsThrows(system,
                  "READ_W_BASIS in/Command/lam_w.bf\n"
                  "READ_H_RGRID in/Command/lam_h.rf\nFINISH\n"));
      TEST_ASSERT(system.h().hasData());
      TEST_ASSERT(!system.h().isSymmetric());
      TEST_ASSERT(!readCommandsThrows(system,
                  "READ_H_BASIS in/Command/lam_h.bf\nFINISH\n"));
      checkSameH(system, ref);

      // READ_H_RGRID as the first command, then READ_H_BASIS
      Rp::System<1,CPT> first;
      initSystem(first);
      TEST_ASSERT(!readCommandsThrows(first,
                  "READ_H_RGRID in/Command/lam_h.rf\n"
                  "READ_W_BASIS in/Command/lam_w.bf\n"
                  "READ_H_BASIS in/Command/lam_h.bf\nFINISH\n"));
      checkSameH(first, ref);

      // The external fields enter the solution
      TEST_ASSERT(!readCommandsThrows(ref, "COMPUTE\nFINISH\n"));
      TEST_ASSERT(!readCommandsThrows(first, "COMPUTE\nFINISH\n"));
      TEST_ASSERT(maxDiff(first.c().rgrid(0), ref.c().rgrid(0)) < 1.0E-10);
   }

   /*
   * WFields::setBasis after an asymmetric setRGrid (the same lazy
   * allocation of the basis arrays).
   */
   void testSetBasisAfterRGrid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testSetBasisAfterRGrid.log");

      Rp::System<1,CPT> ref;
      DArray< DArray<double> > basis;
      DArray< RField<1,CPT> > rgrid;
      initReference(ref, basis, rgrid);

      Rp::System<1,CPT> system;
      initAsymmetricH(system, rgrid);
      bool thrown = false;
      try {
         system.h().setBasis(basis);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
      checkSameH(system, ref);
   }

   /*
   * A symmetric WFields::setRGrid after an asymmetric one.
   */
   void testSetSymmetricRGridAfterRGrid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testSetSymmetricRGridAfterRGrid.log");

      Rp::System<1,CPT> ref;
      DArray< DArray<double> > basis;
      DArray< RField<1,CPT> > rgrid;
      initReference(ref, basis, rgrid);

      Rp::System<1,CPT> system;
      initAsymmetricH(system, rgrid);
      bool thrown = false;
      try {
         system.h().setRGrid(rgrid, true);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
      checkSameH(system, ref);
   }

   /*
   * WFields::symmetrize after an asymmetric setRGrid (it did not
   * allocate the basis arrays).
   */
   void testSymmetrizeAfterRGrid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testSymmetrizeAfterRGrid.log");

      Rp::System<1,CPT> ref;
      DArray< DArray<double> > basis;
      DArray< RField<1,CPT> > rgrid;
      initReference(ref, basis, rgrid);

      Rp::System<1,CPT> system;
      initAsymmetricH(system, rgrid);
      bool thrown = false;
      try {
         system.h().symmetrize();
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
      checkSameH(system, ref);
   }

   // Utility functions

   /*
   * Read the parameter file of a lamellar system (group P_-1).
   */
   void initSystem(Rp::System<1,CPT>& system)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile("in/Command/lam_param", in);
      system.readParam(in);
      in.close();
   }

   /*
   * Reference system with w fields and the external fields of
   * lam_h.bf, and copies of those external fields in basis and r-grid
   * form.
   */
   void initReference(Rp::System<1,CPT>& ref,
                      DArray< DArray<double> >& basis,
                      DArray< RField<1,CPT> >& rgrid)
   {
      initSystem(ref);
      TEST_ASSERT(!readCommandsThrows(ref,
                  "READ_W_BASIS in/Command/lam_w.bf\n"
                  "READ_H_BASIS in/Command/lam_h.bf\nFINISH\n"));
      int nMonomer = ref.mixture().nMonomer();
      basis.allocate(nMonomer);
      rgrid.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         basis[i] = ref.h().basis(i);
         rgrid[i].allocate(ref.domain().mesh().dimensions());
         rgrid[i] = ref.h().rgrid(i);
      }
   }

   /*
   * System with w fields, and with external fields set in r-grid form
   * as not symmetric.
   */
   void initAsymmetricH(Rp::System<1,CPT>& system,
                        DArray< RField<1,CPT> > const & rgrid)
   {
      initSystem(system);
      TEST_ASSERT(!readCommandsThrows(system,
                  "READ_W_BASIS in/Command/lam_w.bf\nFINISH\n"));
      system.h().setRGrid(rgrid, false);
      TEST_ASSERT(system.h().hasData());
      TEST_ASSERT(!system.h().isSymmetric());
   }

   /*
   * Execute a command script given as a string. Returns true if a
   * Util::Exception was thrown.
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

   /*
   * Check that system has the symmetric external fields of ref, in
   * basis and r-grid form.
   */
   void checkSameH(Rp::System<1,CPT> const & system,
                   Rp::System<1,CPT> const & ref)
   {
      TEST_ASSERT(system.h().hasData());
      TEST_ASSERT(system.h().isSymmetric());
      int nMonomer = ref.mixture().nMonomer();
      int nBasis = ref.domain().basis().nBasis();
      for (int i = 0; i < nMonomer; ++i) {
         DArray<double> const & a = system.h().basis(i);
         DArray<double> const & b = ref.h().basis(i);
         TEST_ASSERT(a.capacity() == nBasis);
         for (int j = 0; j < nBasis; ++j) {
            TEST_ASSERT(std::abs(a[j] - b[j]) < 1.0E-12);
         }
         TEST_ASSERT(maxDiff(system.h().rgrid(i), ref.h().rgrid(i))
                     < 1.0E-12);
      }
   }

   /*
   * Maximum absolute difference between two r-grid fields.
   */
   double maxDiff(RField<1,CPT> const & a, RField<1,CPT> const & b)
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

};

TEST_BEGIN(CommandReadHTest)
TEST_ADD(CommandReadHTest, testReadHBasisAfterRGrid)
TEST_ADD(CommandReadHTest, testSetBasisAfterRGrid)
TEST_ADD(CommandReadHTest, testSetSymmetricRGridAfterRGrid)
TEST_ADD(CommandReadHTest, testSymmetrizeAfterRGrid)
TEST_END(CommandReadHTest)

#endif
