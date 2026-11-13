#ifndef RPC_MASK_READ_TRANSACTION_TEST_H
#define RPC_MASK_READ_TRANSACTION_TEST_H

#include <test/UnitTestRunner.h>
#include <rp/system/System.h>
#include <rp/field/WFields.h>
#include <rp/field/Mask.h>
#include <rp/field/CFields.h>
#include <rp/field/Domain.h>
#include <rp/scft/ScftThermo.h>
#include <pscf/backend/HostRecvArray.h>
#include <pscf/backend/Reduce.h>
#include <pscf/chem/PolymerModel.h>
#include <util/tests/LogFileUnitTest.h>
#include <util/param/BracketPolicy.h>
#include <util/signal/Signal.h>
#include <util/misc/Exception.h>
#include <sstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

class MaskReadTransactionTest : public LogFileUnitTest
{
public:
   void setUp()
   {
      BracketPolicy::set(BracketPolicy::Optional);
      PolymerModel::setModel(PolymerModel::Thread);
   }

   void testMalformedBasisPreservesMask()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/MaskReadTransactionBasis.log");
      for (int failure = 0; failure < 4; ++failure) check(true, failure);
   }

   void testMalformedRGridPreservesMask()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/MaskReadTransactionRGrid.log");
      for (int failure = 0; failure < 4; ++failure) check(false, failure);
   }

   void testFirstReadFailureThenRecovery()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/MaskReadTransactionFirst.log");
      Rp::System<1,T> system; initialize(system);
      std::istringstream bad(header(1) + "mesh 8\n.25\nBROKEN\n");
      bool threw = false;
      try { system.mask().readRGrid(bad, true); }
      catch (Exception&) { threw = true; }
      TEST_ASSERT(threw);
      TEST_ASSERT(!system.mask().hasData());
      std::istringstream good(header(1) + grid(.25));
      system.mask().readRGrid(good, true);
      assertMask(system, .25);
      system.compute();
      TEST_ASSERT(std::abs(concentrationMean(system) - .25) < tolerance());
   }

private:
   typedef CPT T;
   struct Counter { int n = 0; void tick() { ++n; } };

   double tolerance() const
   {
#ifdef SINGLE_PRECISION
      return 1e-6;
#else
      return 1e-12;
#endif
   }

   std::string header(int nMonomer, double length = 4.0)
   {
      std::ostringstream out;
      out << "format 1 0\ndim 1\ncrystal_system lamellar\n"
          << "N_cell_param 1\ncell_param " << length
          << "\ngroup_name P_-1\nN_monomer "
          << nMonomer << "\n";
      return out.str();
   }

   std::string grid(double value)
   {
      std::ostringstream out; out << "mesh 8\n";
      for (int i = 0; i < 8; ++i) out << value << "\n";
      return out.str();
   }

   void initialize(Rp::System<1,T>& system)
   {
      std::istringstream parameter(
         "System{ Mixture{ nMonomer 2 monomers[ 1 1 ] nPolymer 1 nSolvent 2 "
         "Polymer{ type linear nBlock 1 blocks[ 0 1 ] phi 0 } "
         "Solvent{ monomerId 0 size 1 phi .5 } "
         "Solvent{ monomerId 1 size 1 phi .5 } ds .05 } "
         "Interaction{ chi( 0 1 10 ) } "
         "Domain{ mesh 8 lattice lamellar groupName P_-1 } }");
      system.readParam(parameter);
      std::istringstream w(header(2) + "N_basis 1\n0 0 0 1\n");
      system.w().readBasis(w);
   }

   void assertMask(Rp::System<1,T>& system, double expected)
   {
      TEST_ASSERT(system.mask().hasData());
      TEST_ASSERT(system.mask().isSymmetric());
      TEST_ASSERT(std::abs(system.mask().basis()[0] - expected) < tolerance());
      TEST_ASSERT(std::abs(system.mask().phiTot() - expected) < tolerance());
      HostRecvArray<typename T::Real,T> host(system.mask().rgrid());
      for (int i = 0; i < host.capacity(); ++i) {
         TEST_ASSERT(std::abs(host[i] - expected) < tolerance());
      }
   }

   double concentrationMean(Rp::System<1,T>& system)
   {
      return (Reduce::sum(system.c().rgrid(0))
            + Reduce::sum(system.c().rgrid(1))) / 8.0;
   }

   void check(bool basis, int failure)
   {
      Rp::System<1,T> system; initialize(system);
      std::istringstream initial(header(1) + "N_basis 1\n.5 0 1\n");
      system.mask().readBasis(initial);
      system.compute(); system.scft().compute();
      assertMask(system, .5);
      Counter count; system.mask().signal().addObserver(count, &Counter::tick);
      std::string data = "BROKEN\n";
      if (failure == 1) {
         data = header(1) + (basis ? "N_basis 2\n.25 0 1\nBROKEN\n"
                                  : "mesh 8\n.25\nBROKEN\n");
      } else if (failure == 2) {
         data = header(1) + (basis ? "N_basis 2\n.25 0 1\n"
                                  : "mesh 8\n.25\n.25\n.25\n.25\n");
      } else if (failure == 3) {
         data = header(1, 6.0) + (basis ? "N_basis 2\n.25 0 1\nBROKEN\n"
                                       : "mesh 8\n.25\nBROKEN\n");
      }
      std::istringstream bad(data);
      bool threw = false;
      try {
         if (basis) system.mask().readBasis(bad);
         else system.mask().readRGrid(bad, true);
      } catch (Exception&) { threw = true; }
      TEST_ASSERT(threw);
      assertMask(system, .5);
      TEST_ASSERT(system.domain().unitCell().parameter(0) == 4.0);
      TEST_ASSERT(count.n == 0);
      TEST_ASSERT(system.c().hasData());
      TEST_ASSERT(system.scft().hasData());
      TEST_ASSERT(std::abs(concentrationMean(system) - .5) < tolerance());
      system.compute();
      TEST_ASSERT(std::abs(concentrationMean(system) - .5) < tolerance());

      std::istringstream good(header(1) + (basis ? "N_basis 1\n.25 0 1\n"
                                               : grid(.25)));
      if (basis) system.mask().readBasis(good);
      else system.mask().readRGrid(good, true);
      assertMask(system, .25);
      TEST_ASSERT(count.n == 1);
      TEST_ASSERT(!system.c().hasData());
      TEST_ASSERT(!system.scft().hasData());
      system.compute(); system.scft().compute();
      TEST_ASSERT(std::abs(concentrationMean(system) - .25) < tolerance());
   }
};

TEST_BEGIN(MaskReadTransactionTest)
TEST_ADD(MaskReadTransactionTest, testMalformedBasisPreservesMask)
TEST_ADD(MaskReadTransactionTest, testMalformedRGridPreservesMask)
TEST_ADD(MaskReadTransactionTest, testFirstReadFailureThenRecovery)
TEST_END(MaskReadTransactionTest)

#endif
