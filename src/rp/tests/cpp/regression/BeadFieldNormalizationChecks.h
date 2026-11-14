#ifndef RP_BEAD_FIELD_NORMALIZATION_CHECKS_H
#define RP_BEAD_FIELD_NORMALIZATION_CHECKS_H

#include <util/tests/LogFileUnitTest.h>
#include <util/param/BracketPolicy.h>
#include <util/misc/Exception.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Block.h>
#include <rp/solvers/Propagator.h>
#include <rp/scft/ScftThermo.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <pscf/backend/VecOp.h>
#include <pscf/chem/PolymerModel.h>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>

namespace Pscf {
namespace Rp {

// A direct, real-space eight-site Gaussian transfer matrix is an independent
// oracle for the FFT solver. Long double keeps the two-level Boltzmann weights
// representable even for contrast 2000. The oracle uses unshifted fields.
template <class T>
class BeadFieldNormalizationChecks : public LogFileUnitTest
{
   typedef std::vector<long double> Row;
   struct Exact { long double q; Row c; };
   Util::BracketPolicy::Type bracketPolicy_;

   Row bond(Row const& input, double length, double fraction = 1.0)
   {
      Row result(8, 0.0L);
      const long double pi = std::acos(-1.0L);
      for (int i = 0; i < 8; ++i) {
         for (int j = 0; j < 8; ++j) {
            long double kernel = 0.0L;
            for (int m = 0; m < 8; ++m) {
               const int k = m <= 4 ? m : m - 8;
               const long double wave = 2*pi*k/length;
               kernel += std::exp(-fraction*wave*wave/6)
                       * std::cos(2*pi*k*(i-j)/8)/8;
            }
            result[i] += kernel*input[j];
         }
      }
      return result;
   }

   Exact exact(int count, double contrast, double length, bool branched)
   {
      Row weight(8), density(8, 0.0L);
      for (int i = 0; i < 8; ++i) {
         weight[i] = std::exp(-(i % 2 ? (long double)contrast : 0.0L));
      }
      long double q = 0.0L;
      if (branched) {
         // Three arms attached to one junction: (1,1,1) or (2,1,1).
         // Integrate over the junction independently, with one half bond
         // on each arm and an extra full bond in the two-bead arm.
         Row arm = bond(weight, length, .5);
         Row longEnd = bond(weight, length);
         for (int i = 0; i < 8; ++i) longEnd[i] *= weight[i];
         Row longArm = count == 4 ? bond(longEnd, length, .5) : arm;
         Row other(8);
         for (int i = 0; i < 8; ++i) {
            q += longArm[i]*arm[i]*arm[i]/8;
            other[i] = arm[i]*arm[i];
         }
         other = bond(other, length, .5);
         if (count == 3) {
            for (int i = 0; i < 8; ++i) density[i] = weight[i]*other[i]/q;
         } else {
            Row terminal(8), shortArm(8);
            for (int i = 0; i < 8; ++i) {
               density[i] = longEnd[i]*other[i];
               terminal[i] = weight[i]*other[i];
               shortArm[i] = longArm[i]*arm[i];
            }
            terminal = bond(terminal, length);
            shortArm = bond(shortArm, length, .5);
            for (int i = 0; i < 8; ++i) {
               density[i] = (density[i] + weight[i]*(terminal[i]+2*shortArm[i]))/(4*q);
            }
         }
      } else {
         std::vector<Row> forward(count, Row(8));
         forward[0] = weight;
         for (int j = 1; j < count; ++j) {
            forward[j] = bond(forward[j-1], length);
            for (int i = 0; i < 8; ++i) forward[j][i] *= weight[i];
         }
         for (int i = 0; i < 8; ++i) q += forward[count-1][i]/8;
         for (int j = 0; j < count; ++j) {
            for (int i = 0; i < 8; ++i) {
               density[i] += forward[j][i]*forward[count-1-j][i]
                            /weight[i]/(count*q);
            }
         }
      }
      return Exact{q, density};
   }

   void check(int count, double contrast, double gauge, bool open,
              double phi, bool branched = false, bool batched = false)
   {
      const double length = 4.0;
      Exact expected = exact(count, contrast, length, branched);
      const double lnQ = (double)std::log(expected.q) - count*gauge;
      std::ostringstream text;
      text << std::setprecision(17)
           << "System{ polymerModel bead Mixture{ nMonomer 2 monomers[ 1 1 ] "
           << "nPolymer 1 nSolvent 1 Polymer{ type "
           << (branched ? "branched nBlock 3 blocks[ 0 " + std::to_string(count-2)
                             + " 0 1 0 1 0 2 0 1 0 3 ] "
                        : "linear nBlock 1 blocks[ 0 " + std::to_string(count) + " ] ")
           << (open ? "mu " : "phi ") << (open ? std::log(phi)-lnQ : phi)
           << " } Solvent{ monomerId 1 size 1 phi " << 1.0-phi
           << " } useBatchedFFT " << int(batched)
           << " } Interaction{ chi( 0 1 2 ) } "
           << "Domain{ mesh 8 lattice lamellar } }";
      System<1,T> system;
      std::istringstream in(text.str()); system.readParam(in);
      Util::FSArray<double,6> cell; cell.append(length); system.setUnitCell(cell);
      Util::DArray<RField<1,T> > fields; fields.allocate(2);
      Util::DArray<typename T::Real> host; host.allocate(8);
      for (int m = 0; m < 2; ++m) {
         fields[m].allocate(IntVec<1>(8));
         for (int i = 0; i < 8; ++i) {
            host[i] = m ? 0.0 : gauge + (i % 2 ? contrast : 0.0);
         }
         VecOp::eqV(fields[m], host);
      }
      system.w().setRGrid(fields);
      system.compute();
      auto const& polymer = system.mixture().polymer(0);
      TEST_ASSERT(std::isfinite(polymer.lnQ()));
      TEST_ASSERT(std::abs(polymer.lnQ() - lnQ) < 2e-9);
      TEST_ASSERT(std::abs(polymer.phi() - phi) < 2e-10);
      VecOp::eqV(host, system.c().rgrid(0));
      for (int i = 0; i < 8; ++i) {
         TEST_ASSERT(std::isfinite(host[i]));
         TEST_ASSERT(std::abs(host[i] - phi*(double)expected.c[i]) < 2e-10);
      }
      VecOp::eqV(host, system.c().rgrid(1));
      for (int i = 0; i < 8; ++i) TEST_ASSERT(std::abs(host[i]-(1-phi)) < 1e-12);

      // Every junction/direction must give the same actual shifted Q.
      double lnShift = 0.0;
      for (int b = 0; b < polymer.nBlock(); ++b) {
         lnShift += polymer.block(b).wShift()*polymer.block(b).nBead();
      }
      for (int b = 0; b < polymer.nBlock(); ++b) {
         for (int d = 0; d < 2; ++d) {
            double q;
            polymer.block(b).propagator(d).computeQ(q);
            TEST_ASSERT(std::isfinite(q) && q > 0.0);
            TEST_ASSERT(std::abs(std::log(q)-lnQ-lnShift) < 2e-9);
         }
      }
      if (count == 1) {
         VecOp::eqV(host, polymer.block(0).propagator(0).q(1));
         for (int i = 0; i < 8; ++i) {
            const double value = std::exp(-(gauge + (i % 2 ? contrast : 0.0)
                                           -polymer.block(0).wShift()));
            TEST_ASSERT(std::abs(host[i]-value) <= 1e-12*std::max(1.0,value));
         }
      }
      system.computeStress();
      const double h = 1e-4;
      const double derivative = (double)(std::log(exact(count, contrast, length+h, branched).q)
                                       -std::log(exact(count, contrast, length-h, branched).q))/(2*h);
      const double stress = -phi*derivative/count;
      TEST_ASSERT(std::isfinite(polymer.stress(0)));
      TEST_ASSERT(std::abs(polymer.stress(0)-stress) < 2e-8);
      system.scft().compute();
      TEST_ASSERT(std::isfinite(system.scft().fHelmholtz()));
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
   void checkLinear()
   {
      this->openLogFile("out/BeadFieldNormalizationLinear.log");
      // Ordinary, either side of the product guard, and overflow fields.
      for (double contrast : {2.0, 450.0, 455.0, 600.0, 800.0, 2000.0}) {
         check(1, contrast, 0.0, false, .5);
      }
      check(1, 2000.0, 0.0, false, 0.0);
      check(1, 2000.0, 1000.0, true, .5);
      check(1, 2000.0, -1000.0, false, .5);
      check(3, 2.0, 0.0, false, .5);
      check(3, 600.0, 0.0, false, .5);
      check(20, 100.0, 0.0, false, .5);
      check(3, 2000.0, 1000.0, true, .5);
      // The shifted Q underflows; the shifts are raised until it is one.
      check(2000, 100.0, 0.0, false, .5);
   }
   void checkBranched()
   {
      this->openLogFile("out/BeadFieldNormalizationBranched.log");
      for (bool batched : {false, true}) {
         check(3, 2.0, 0.0, false, .5, true, batched);
         check(3, 800.0, 0.0, false, .5, true, batched);
         check(3, 2000.0, 1000.0, true, .5, true, batched);
         check(4, 2000.0, 0.0, false, .5, true, batched);
      }
   }
   void checkRangeAndRecovery()
   {
      this->openLogFile("out/BeadFieldNormalizationRange.log");
      for (int mode = 1; mode < 2; ++mode) {
         // Finite fields but an unrepresentable log Q. (A shifted Q
         // below the supported range, 2000 beads in 0|100, is now
         // computed by raising the shifts; see checkLinear.)
         System<1,T> system;
         std::ostringstream text;
         text << "System{ polymerModel bead Mixture{ nMonomer 2 monomers[ 1 1 ] "
              << "nPolymer 1 nSolvent 1 Polymer{ type linear nBlock 1 blocks[ 0 "
              << (mode ? 3 : 2000) << " ] phi .5 } "
              << "Solvent{ monomerId 1 size 1 phi .5 } } "
              << "Interaction{ chi( 0 1 2 ) } Domain{ mesh 8 lattice lamellar } }";
         std::istringstream in(text.str()); system.readParam(in);
         Util::FSArray<double,6> cell; cell.append(4.0); system.setUnitCell(cell);
         Util::DArray<RField<1,T> > fields; fields.allocate(2);
         Util::DArray<typename T::Real> host; host.allocate(8);
         for (int m = 0; m < 2; ++m) {
            fields[m].allocate(IntVec<1>(8));
            for (int i = 0; i < 8; ++i) {
               host[i] = m ? 0.0 : (mode ? 1e308 : (i % 2 ? 100.0 : 0.0));
            }
            VecOp::eqV(fields[m], host);
         }
         system.w().setRGrid(fields);
         bool threw = false;
         try { system.compute(); }
         catch (Util::Exception& e) {
            threw = e.message().find("normalization exceeds floating-point range")
                  != std::string::npos;
         }
         TEST_ASSERT(threw);
         TEST_ASSERT(!system.c().hasData());
         for (int m = 0; m < 2; ++m) VecOp::eqS(fields[m], 0.0);
         system.w().setRGrid(fields);
         system.compute();
         system.computeStress();
         for (int m = 0; m < 2; ++m) {
            VecOp::eqV(host, system.c().rgrid(m));
            for (int i = 0; i < 8; ++i) TEST_ASSERT(std::abs(host[i]-.5) < 1e-10);
         }
      }
   }
};

}
}
#endif
