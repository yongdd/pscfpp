#ifndef RP_SOLVENT_NORMALIZATION_CHECKS_H
#define RP_SOLVENT_NORMALIZATION_CHECKS_H

#include <test/UnitTest.h>
#include <rp/solvers/Solvent.h>
#include <pscf/field/RField.h>
#include <pscf/backend/VecOp.h>
#include <pscf/mesh/Mesh.h>
#include <util/containers/DArray.h>
#include <util/tests/LogFileUnitTest.h>
#include <sstream>
#include <iomanip>
#include <cmath>

namespace Pscf {
namespace Rp {

// Exact two-level Boltzmann weights, including an occupied-volume correction.
template <class T>
class SolventNormalizationChecks : public LogFileUnitTest
{
protected:
   void check(double size, double contrast, double gauge, bool open,
              double amount, double phiTot = 1.0)
   {
      IntVec<1> dimensions(8);
      Mesh<1> mesh(dimensions);
      Solvent<1,T> solvent;
      std::ostringstream text;
      text << std::setprecision(17)
           << "Solvent{ monomerId 0 size " << size
           << (open ? " mu " : " phi ")
           << (open ? size*gauge + amount : amount) << " }";
      std::istringstream input(text.str());
      solvent.readParam(input);
      solvent.associate(mesh);
      solvent.allocate();
      RField<1,T> w;
      w.allocate(dimensions);
      Util::DArray<typename T::Real> host;
      host.allocate(8);
      for (int i = 0; i < 8; ++i) {
         host[i] = gauge + (i % 2 ? contrast : 0.0);
      }
      VecOp::eqV(w, host);
      solvent.compute(w, phiTot);
      VecOp::eqV(host, solvent.cField());

      // Only nonpositive exponentials occur in this independent oracle.
      double high = std::exp(-size*contrast);
      double qScaled = (1.0 + high)/(2.0*phiTot);
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
      for (int i = 0; i < 8; ++i) {
         double weight = i % 2 ? high : 1.0;
         double expected = expectedPhi*weight/qScaled;
         TEST_ASSERT(std::isfinite(host[i]));
         TEST_ASSERT(std::abs(host[i] - expected) < 1e-10);
      }
   }

   void checkClosed()
   {
      openLogFile("out/SolventNormalizationClosed.log");
      check(1.0, 2.0, 0.0, false, 0.5);
      check(700.0, 2.0, 0.0, false, 0.5);
      check(1000.0, 2.0, 0.0, false, 0.5);
      check(1000.0, 0.0, 1.0, false, 0.5);
      check(1000.0, 2.0, 20.0, false, 0.5);
      check(1000.0, 2.0, -20.0, false, 0.5);
      check(1000.0, 2.0, 0.0, false, 0.0);
      check(1000.0, 2.0, 0.0, false, 0.5, 0.25);
      closeLogFile();
   }

   void checkOpen()
   {
      openLogFile("out/SolventNormalizationOpen.log");
      check(1.0, 2.0, 0.0, true, 0.0);
      check(1000.0, 2.0, 0.0, true, 0.0);
      check(1000.0, 2.0, 20.0, true, 0.0);
      check(1000.0, 2.0, -20.0, true, 0.0);
      check(1000.0, 0.0, 1.0, true, std::log(0.5));
      check(1000.0, 2.0, 0.0, true, -800.0);
      check(1000.0, 2.0, 0.0, true, 0.0, 0.25);
      closeLogFile();
   }
};

}
}
#endif
