#ifndef RPG_FIELD_READ_INVALIDATION_TEST_H
#define RPG_FIELD_READ_INVALIDATION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <rp/system/System.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <rp/field/Domain.h>
#include <rp/scft/ScftThermo.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/chem/PolymerModel.h>
#include <util/tests/LogFileUnitTest.h>
#include <util/param/BracketPolicy.h>
#include <util/misc/Exception.h>
#include <sstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

class FieldReadInvalidationTest : public LogFileUnitTest
{
public:
   void setUp()
   {
      BracketPolicy::set(BracketPolicy::Optional);
      PolymerModel::setModel(PolymerModel::Thread);
   }

   void testMalformedBasisReplacement()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FieldReadInvalidationBasis.log");
      for (int failure = 0; failure < 3; ++failure) check(true, failure);
   }

   void testMalformedRGridReplacement()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FieldReadInvalidationRGrid.log");
      for (int failure = 0; failure < 3; ++failure) check(false, failure);
   }

private:
   typedef CUT T;

   std::string header(double length)
   {
      std::ostringstream out;
      out << "format 1 0\ndim 1\ncrystal_system lamellar\nN_cell_param 1\n"
          << "cell_param " << length << "\ngroup_name P_-1\nN_monomer 2\n";
      return out.str();
   }

   std::string field(bool basis, double value)
   {
      std::ostringstream out;
      if (basis) out << "N_basis 1\n" << value << " " << value << " 0 1\n";
      else {
         out << "mesh 16\n";
         for (int i = 0; i < 16; ++i) out << value << " " << value << "\n";
      }
      return out.str();
   }

   void read(Rp::System<1,T>& system, bool basis, std::string const& text)
   {
      std::istringstream input(text);
      if (basis) system.w().readBasis(input);
      else system.w().readRGrid(input, true);
   }

   void check(bool basis, int failure)
   {
      Rp::System<1,T> system;
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1 1 ] nPolymer 1 "
         "Polymer{ type linear nBlock 2 blocks[ 0 .5 1 .5 ] phi 1 } ds .05 } "
         "Interaction{ chi( 0 1 15 ) } "
         "Domain{ mesh 16 lattice lamellar groupName P_-1 } "
         "AmIteratorBasis{ epsilon 1e-8 maxItr 10 maxHist 4 isFlexible 0 } }");
      system.readParam(param);
      read(system, true, header(4.0) + field(true, 5.0));
      system.compute(true);
      system.scft().compute();
      TEST_ASSERT(std::abs(system.scft().fHelmholtz() - 2.75) < 1e-11);

      // Fail before a header is accepted, after partial data, and after
      // both a changed cell and partial data. No iteration is required.
      const double length = failure == 2 ? 6.0 : 4.0;
      std::string bad = "BROKEN\n";
      if (failure != 0) {
         bad = header(length) + (basis ? "N_basis 2\n7 7 0 1\nBROKEN\n"
                                       : "mesh 16\n7 7\nBROKEN\n");
      }
      bool threw = false;
      try { read(system, basis, bad); }
      catch (Exception&) { threw = true; }
      TEST_ASSERT(threw);
      TEST_ASSERT(!system.w().hasData());
      TEST_ASSERT(!system.w().isSymmetric());
      TEST_ASSERT(!system.c().hasData());
      TEST_ASSERT(!system.scft().hasData());
      if (failure == 2) {
         TEST_ASSERT(std::abs(system.domain().unitCell().parameter(0) - 6.0) < 1e-12);
      }
      threw = false;
      try { system.compute(); }
      catch (Exception&) { threw = true; }
      TEST_ASSERT(threw);
      std::ostringstream output;
      threw = false;
      try { system.w().writeBasis(output); }
      catch (Exception&) { threw = true; }
      TEST_ASSERT(threw);

      // Successful replacement recovers coherent fields and the analytic
      // homogeneous free energy, independent of a common field shift.
      read(system, basis, header(length) + field(basis, 7.0));
      TEST_ASSERT(system.w().hasData());
      TEST_ASSERT(system.w().isSymmetric());
      system.compute(true);
      system.scft().compute();
      TEST_ASSERT(std::abs(system.scft().fHelmholtz() - 2.75) < 1e-11);
   }
};

TEST_BEGIN(FieldReadInvalidationTest)
TEST_ADD(FieldReadInvalidationTest, testMalformedBasisReplacement)
TEST_ADD(FieldReadInvalidationTest, testMalformedRGridReplacement)
TEST_END(FieldReadInvalidationTest)

#endif
