#ifndef RP_FTS_RANGE_FAILURE_CHECKS_H
#define RP_FTS_RANGE_FAILURE_CHECKS_H

#include <util/tests/LogFileUnitTest.h>
#include <util/param/BracketPolicy.h>
#include <util/misc/Exception.h>
#include <util/containers/DArray.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <rp/fts/compressor/AmCompressor.h>
#include <rp/fts/compressor/LrAmCompressor.h>
#include <rp/fts/compressor/LrCompressor.h>
#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/fts/montecarlo/McMove.h>
#include <rp/fts/montecarlo/McMoveManager.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <pscf/field/RField.h>
#include <pscf/backend/VecOp.h>
#include <pscf/chem/PolymerModel.h>
#include <sstream>
#include <string>
#include <cmath>

namespace Pscf {
namespace Rp {

/*
* A solver result that is not representable during compression must
* end compress() with failure, as a NaN error does, and not with an
* exception: a MC move or BD step then fails, and the simulation goes
* on from the restored state. A bead homopolymer of 1000 beads in a
* field of magnitude 1e304 or more has a ln Q beyond the double range,
* whatever field shift the solver uses.
*/
template <class T>
class FtsRangeFailureChecks : public LogFileUnitTest
{

   Util::BracketPolicy::Type bracketPolicy_;

   /*
   * Bead homopolymer A (1000 beads, phi = phiA) and solvent B (size 1),
   * chi_AB = 2, mesh 8, lamellar cell of length 4. A homogeneous w is
   * compressed for these closed species.
   */
   void read(System<1,T>& system, double phiA)
   {
      std::ostringstream text;
      text << "System{ polymerModel bead Mixture{ nMonomer 2 "
           << "monomers[ 1 1 ] nPolymer 1 nSolvent 1 "
           << "Polymer{ type linear nBlock 1 blocks[ 0 1000 ] phi "
           << phiA << " } Solvent{ monomerId 1 size 1 phi "
           << 1.0 - phiA << " } } Interaction{ chi( 0 1 2 ) } "
           << "Domain{ mesh 8 lattice lamellar } }";
      std::istringstream in(text.str());
      system.readParam(in);
      Util::FSArray<double,6> cell;
      cell.append(4.0);
      system.setUnitCell(cell);
   }

   /*
   * Set homogeneous w fields.
   */
   void setFields(System<1,T>& system, double wA, double wB)
   {
      Util::DArray<RField<1,T> > fields;
      fields.allocate(2);
      Util::DArray<typename T::Real> host;
      host.allocate(8);
      for (int m = 0; m < 2; ++m) {
         fields[m].allocate(IntVec<1>(8));
         for (int i = 0; i < 8; ++i) host[i] = m ? wB : wA;
         VecOp::eqV(fields[m], host);
      }
      system.w().setRGrid(fields);
   }

   /*
   * Are the c fields those of the current w fields?
   */
   bool hasCurrentCFields(System<1,T>& system)
   {
      if (!system.c().hasData()) return false;
      Util::DArray<typename T::Real> before, after;
      before.allocate(16);
      after.allocate(16);
      for (int m = 0; m < 2; ++m) {
         Util::DArray<typename T::Real> host;
         host.allocate(8);
         VecOp::eqV(host, system.c().rgrid(m));
         for (int i = 0; i < 8; ++i) before[8*m + i] = host[i];
      }
      system.compute();
      for (int m = 0; m < 2; ++m) {
         Util::DArray<typename T::Real> host;
         host.allocate(8);
         VecOp::eqV(host, system.c().rgrid(m));
         for (int i = 0; i < 8; ++i) after[8*m + i] = host[i];
      }
      for (int i = 0; i < 16; ++i) {
         if (!(std::abs(before[i] - after[i]) < 1.0E-10)) return false;
      }
      return true;
   }

public:

   void setUp()
   {
      bracketPolicy_ = Util::BracketPolicy::get();
      Util::BracketPolicy::set(Util::BracketPolicy::Optional);
   }

   void tearDown()
   {
      LogFileUnitTest::tearDown();
      PolymerModel::setModel(PolymerModel::Thread);
      Util::BracketPolicy::set(bracketPolicy_);
   }

   /*
   * compress() returns 1 for a range failure of the MDE solution,
   * also when repeated, and converges afterwards from valid fields.
   */
   template <class CompressorT>
   void checkCompressor(std::string const & block)
   {
      System<1,T> system;
      read(system, 0.5);
      CompressorT compressor(system);
      std::istringstream in(block);
      compressor.readParam(in);
      for (int i = 0; i < 2; ++i) {
         setFields(system, 1.0E306, 0.0);
         int error = 0;
         try {
            error = compressor.compress();
         } catch (Util::Exception& e) {
            error = -1;
         }
         TEST_ASSERT(error == 1);
      }
      setFields(system, 0.0, 0.0);
      TEST_ASSERT(compressor.compress() == 0);
      TEST_ASSERT(hasCurrentCFields(system));
   }

   void checkCompressors()
   {
      checkCompressor< AmCompressor<1,T> >(
                         "AmCompressor{ epsilon 1e-8 maxItr 20 }");
      checkCompressor< LrAmCompressor<1,T> >(
                         "LrAmCompressor{ epsilon 1e-8 maxItr 20 }");
      checkCompressor< LrCompressor<1,T> >(
                         "LrCompressor{ epsilon 1e-8 maxItr 20 }");
   }

   /*
   * MC moves of magnitude 1e306 fail with range failures in the
   * compressor. Every move fails, and the simulation completes with
   * the c fields of the restored w fields.
   */
   void checkMcMove(std::string const & move)
   {
      System<1,T> system;
      read(system, 0.3);
      McSimulator<1,T> simulator(system);
      std::istringstream in("McSimulator{ seed 1784106389 "
                            "McMoveManager{ " + move + " } "
                            "AmCompressor{ epsilon 1e-8 maxItr 20 } }");
      simulator.readParam(in);
      setFields(system, 0.0, 0.0);
      bool thrown = false;
      try {
         simulator.simulate(3);
      } catch (Util::Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
      McMove<1,T>& mcMove = simulator.mcMoveManager()[0];
      TEST_ASSERT(mcMove.nAttempt() == 3);
      TEST_ASSERT(mcMove.nFail() == 3);
      TEST_ASSERT(hasCurrentCFields(system));
   }

   void checkMcMoves()
   {
      checkMcMove("RealMove{ probability 1.0 sigma 1e306 }");
      checkMcMove("ForceBiasMove{ probability 1.0 mobility 1e306 }");
   }

   /*
   * As checkMcMove, for BD steps.
   */
   void checkBdStep()
   {
      System<1,T> system;
      read(system, 0.3);
      BdSimulator<1,T> simulator(system);
      std::istringstream in("BdSimulator{ seed 5 "
                            "LMBdStep{ mobility 1e306 } "
                            "AmCompressor{ epsilon 1e-8 maxItr 20 } }");
      simulator.readParam(in);
      setFields(system, 0.0, 0.0);
      bool thrown = false;
      try {
         simulator.simulate(2);
      } catch (Util::Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
      TEST_ASSERT(simulator.iStep() == 0);
      TEST_ASSERT(hasCurrentCFields(system));
   }

};

}
}
#endif
