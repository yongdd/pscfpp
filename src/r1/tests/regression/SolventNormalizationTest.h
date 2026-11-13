#ifndef R1_SOLVENT_NORMALIZATION_TEST_H
#define R1_SOLVENT_NORMALIZATION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <r1/solvers/Solvent.h>
#include <r1/field/Domain.h>
#include <util/containers/DArray.h>
#include <util/tests/LogFileUnitTest.h>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace Pscf;
using namespace Pscf::R1;

// Nine planar nodes give equal measures to alternating fields, including endpoints.
class SolventNormalizationTest : public LogFileUnitTest
{
protected:
   void check(double size, double contrast, double gauge, bool open,
              double amount)
   {
      Domain domain;
      domain.setPlanarParameters(0.0, 1.0, 9);
      Solvent solvent;
      std::ostringstream text;
      text << std::setprecision(17)
           << "Solvent{ monomerId 0 size " << size
           << (open ? " mu " : " phi ")
           << (open ? size*gauge + amount : amount) << " }";
      std::istringstream input(text.str());
      solvent.readParam(input);
      solvent.setDiscretization(domain);
      Util::DArray<double> w;
      w.allocate(9);
      for (int i = 0; i < 9; ++i) {
         w[i] = gauge + (i % 2 ? contrast : 0.0);
      }
      solvent.compute(w);
      Solvent::FieldT const & host = solvent.cField();

      // Only nonpositive exponentials occur in this independent oracle.
      double high = std::exp(-size*contrast);
      double qScaled = (1.0 + high)/2.0;
      double expectedLnQ = std::log(qScaled) - size*gauge;
      double expectedPhi = open ? std::exp(amount)*qScaled : amount;
      TEST_ASSERT(std::isfinite(solvent.lnQ()));
      TEST_ASSERT(std::abs(solvent.lnQ() - expectedLnQ) < 1e-10);
      TEST_ASSERT(std::isfinite(solvent.phi()));
      TEST_ASSERT(std::abs(solvent.phi() - expectedPhi) < 1e-10);
      if (open || amount > 0.0) {
         double expectedMu = open ? size*gauge + amount
                                  : std::log(amount) - expectedLnQ;
         TEST_ASSERT(std::isfinite(solvent.mu()));
         TEST_ASSERT(std::abs(solvent.mu() - expectedMu) < 1e-10);
      } else {
         TEST_ASSERT(std::isinf(solvent.mu()) && solvent.mu() < 0.0);
      }
      for (int i = 0; i < 9; ++i) {
         double weight = i % 2 ? high : 1.0;
         double expected = expectedPhi*weight/qScaled;
         TEST_ASSERT(std::isfinite(host[i]));
         TEST_ASSERT(std::abs(host[i] - expected) < 1e-10);
      }
   }

public:
   void testClosed()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SolventNormalizationClosed.log");
      check(1.0, 2.0, 0.0, false, 0.5);
      check(700.0, 2.0, 0.0, false, 0.5);
      check(1000.0, 2.0, 0.0, false, 0.5);
      check(1000.0, 0.0, 1.0, false, 0.5);
      check(1000.0, 2.0, 20.0, false, 0.5);
      check(1000.0, 2.0, -20.0, false, 0.5);
      check(1000.0, 2.0, 0.0, false, 0.0);
   }

   void testOpen()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SolventNormalizationOpen.log");
      check(1.0, 2.0, 0.0, true, 0.0);
      check(1000.0, 2.0, 0.0, true, 0.0);
      check(1000.0, 2.0, 20.0, true, 0.0);
      check(1000.0, 2.0, -20.0, true, 0.0);
      check(1000.0, 0.0, 1.0, true, std::log(0.5));
      check(1000.0, 2.0, 0.0, true, -800.0);
   }

   void testFieldSize()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SolventNormalizationFieldSize.log");
      Domain domain;
      domain.setPlanarParameters(0.0, 1.0, 9);
      Solvent solvent;
      std::istringstream input("Solvent{ monomerId 0 size 1000 phi .5 }");
      solvent.readParam(input);
      solvent.setDiscretization(domain);
      for (int size : {8, 10}) {
         Util::DArray<double> w;
         w.allocate(size);
         for (int i = 0; i < size; ++i) w[i] = 0.0;
         bool caught = false;
         try {
            solvent.compute(w);
         } catch (Util::Exception const &) {
            caught = true;
         }
         TEST_ASSERT(caught);
      }
      // A rejected field must not prevent a later valid computation.
      Util::DArray<double> w;
      w.allocate(9);
      for (int i = 0; i < 9; ++i) w[i] = 0.0;
      solvent.compute(w);
      TEST_ASSERT(std::abs(solvent.lnQ()) < 1e-12);
      for (int i = 0; i < 9; ++i) {
         TEST_ASSERT(std::abs(solvent.cField()[i] - .5) < 1e-12);
      }
   }
};

TEST_BEGIN(SolventNormalizationTest)
TEST_ADD(SolventNormalizationTest, testClosed)
TEST_ADD(SolventNormalizationTest, testOpen)
TEST_ADD(SolventNormalizationTest, testFieldSize)
TEST_END(SolventNormalizationTest)

#endif
