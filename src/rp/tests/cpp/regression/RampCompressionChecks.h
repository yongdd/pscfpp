#ifndef RP_RAMP_COMPRESSION_CHECKS_H
#define RP_RAMP_COMPRESSION_CHECKS_H

#include <test/UnitTest.h>
#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <pscf/field/RField.h>
#include <pscf/backend/VecOp.h>
#include <util/containers/FSArray.h>
#include <util/tests/LogFileUnitTest.h>
#include <sstream>
#include <iomanip>
#include <cmath>

namespace Pscf {
namespace Rp {

// Homogeneous open homopolymers make the partial saddle analytic.
template <class T>
class RampCompressionChecks : public LogFileUnitTest
{
protected:
   void check(bool bd, int nStep, bool failed, bool zeroChange)
   {
      double change = zeroChange ? 0.0 : std::log(3.0);
      System<1,T> system;
      std::ostringstream text;
      text << std::setprecision(17)
           << "System{ Mixture{ nMonomer 2 monomers[ 1 1 ] nPolymer 2 "
           << "Polymer{ type linear nBlock 1 blocks[ 0 1 ] mu "
           << std::log(0.5) << " } "
           << "Polymer{ type linear nBlock 1 blocks[ 1 1 ] mu "
           << std::log(0.5) << " } vMonomer .01 ds .05 } "
           << "Interaction{ chi( 0 1 10 ) } "
           << "Domain{ mesh 8 lattice lamellar } "
           << (bd ? "BdSimulator{ seed 1234 ExplicitBdStep{ mobility 0 } "
                  : "McSimulator{ seed 1234 McMoveManager{ "
                    "RealMove{ probability 1 sigma 0 } } ")
           << "LrCompressor{ epsilon 1e-10 maxItr "
           << (failed || zeroChange ? 1 : 100)
           << " verbose 0 errorType max } "
           << "LinearRamp{ nParameter 1 parameters[ mu_polymer 0 "
           << change << " ] } } }";
      std::istringstream input(text.str());
      system.readParam(input);
      Util::FSArray<double,6> cell;
      cell.append(1.0);
      system.setUnitCell(cell);
      Util::DArray<RField<1,T> > w;
      w.allocate(2);
      IntVec<1> dimensions(8);
      for (int i = 0; i < 2; ++i) {
         w[i].allocate(dimensions);
         VecOp::eqS(w[i], 0.0);
      }
      system.w().setRGrid(w);
      bool caught = false;
      try {
         system.simulate(nStep);
      } catch (Util::Exception const &) {
         caught = true;
      }
      TEST_ASSERT(caught == failed);
      if (failed) {
         TEST_ASSERT(system.simulator().iStep() == 0);
         TEST_ASSERT(!system.c().hasData());
         TEST_ASSERT(!system.simulator().hasHamiltonian());
         TEST_ASSERT(!system.simulator().hasWc());
         TEST_ASSERT(!system.simulator().hasCc());
         TEST_ASSERT(!system.simulator().hasDc());
         // The ramp still uses attempted-step number, not accepted count.
         TEST_ASSERT(std::abs(system.mixture().polymer(0).mu()
                             - std::log(0.5) - change/nStep) < 1e-12);
      } else {
         TEST_ASSERT(system.simulator().iStep() == nStep);
         TEST_ASSERT(system.simulator().iTotalStep() == nStep);
         TEST_ASSERT(system.c().hasData());
         double pressure = std::log((std::exp(change) + 1.0)/2.0);
         double expectedH = 100.0*(-1.0 - pressure + 2.5);
         TEST_ASSERT(system.simulator().hasHamiltonian());
         TEST_ASSERT(std::abs(system.simulator().hamiltonian() - expectedH)
                     < 1e-8);
         Util::DArray<typename T::Real> host;
         host.allocate(8);
         for (int j = 0; j < 2; ++j) {
            VecOp::eqV(host, system.c().rgrid(j));
            double expected = std::exp(std::log(0.5)
                                      + (j == 0 ? change : 0.0) - pressure);
            for (int i = 0; i < 8; ++i) {
               TEST_ASSERT(std::abs(host[i] - expected) < 1e-8);
            }
         }
      }
   }

   void checkSimulator(bool bd)
   {
      // Expected compressor failures write to the global Log stream. Own a
      // live stream instead of relying on a previous test's log lifetime.
      openLogFile(bd ? "out/RampCompressionFailureBd.log"
                     : "out/RampCompressionFailureMc.log");
      for (int nStep = 1; nStep <= 2; ++nStep) {
         check(bd, nStep, true, false);
         check(bd, nStep, false, false);
         check(bd, nStep, false, true);
      }
      closeLogFile();
   }
};

}
}
#endif
