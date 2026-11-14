#ifndef RP_SCFT_FIELD_SHIFT_CHECKS_H
#define RP_SCFT_FIELD_SHIFT_CHECKS_H

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
#include <string>
#include <vector>
#include <cmath>

namespace Pscf {
namespace Rp {

/*
* Field shifts of the polymer MDE solver (Polymer::compute).
*
* Each block is solved in w - s, s = min(w) - ln<exp(-n(w - min(w)))>/n,
* with n the total length or number of beads; if the shifted partition
* function is too small, all shifts are raised so that it becomes one.
* Thread references for mu, c and f come from an independent numpy
* pseudo-spectral (RQM4) solver in the unshifted fields; the stresses come
* from the unfixed devel code in the unshifted fields, which agrees with
* that solver to all 12 printed digits of mu, c and f. The other thread
* references come from the solver with the average shift (alternating
* field) or from devel (field 0|2|4). Bead references come from a direct
* long double transfer matrix.
*/
template <class T>
class ScftFieldShiftChecks : public LogFileUnitTest
{

   typedef std::vector<long double> Row;

   Util::BracketPolicy::Type bracketPolicy_;

   /*
   * Read a 1D lamellar system with one homopolymer (monomer 0) and
   * either a second homopolymer or a solvent (monomer 1).
   */
   void makeSystem(System<1,T>& system, bool bead, double length,
                   double ds, int mesh, double cell, bool solvent,
                   double length1 = -1.0)
   {
      if (length1 < 0.0) length1 = length;
      std::ostringstream text;
      text << std::setprecision(17) << "System{ ";
      if (bead) text << "polymerModel bead ";
      text << "Mixture{ nMonomer 2 monomers[ 1 1 ] "
           << "nPolymer " << (solvent ? 1 : 2) << " "
           << (solvent ? "nSolvent 1 " : "")
           << "Polymer{ type linear nBlock 1 blocks[ 0 " << length
           << " ] phi 0.5 } ";
      if (solvent) {
         text << "Solvent{ monomerId 1 size 1 phi 0.5 } ";
      } else {
         text << "Polymer{ type linear nBlock 1 blocks[ 1 " << length1
              << " ] phi 0.5 } ";
      }
      if (!bead) text << "ds " << ds << " ";
      text << "} Interaction{ chi( 0 1 0 ) } "
           << "Domain{ mesh " << mesh << " lattice lamellar } }";
      std::istringstream in(text.str());
      system.readParam(in);
      Util::FSArray<double,6> cell1;
      cell1.append(cell);
      system.setUnitCell(cell1);
   }

   /*
   * Set the w fields of a system from host values.
   */
   void setFields(System<1,T>& system, std::vector<double> const & w0,
                  std::vector<double> const & w1)
   {
      const int n = (int) w0.size();
      Util::DArray<RField<1,T> > fields;
      fields.allocate(2);
      Util::DArray<typename T::Real> host;
      host.allocate(n);
      for (int m = 0; m < 2; ++m) {
         fields[m].allocate(IntVec<1>(n));
         for (int i = 0; i < n; ++i) host[i] = m ? w1[i] : w0[i];
         VecOp::eqV(fields[m], host);
      }
      system.w().setRGrid(fields);
   }

   /*
   * The documented shift, min(w) - ln<exp(-n(w - min(w)))>/n.
   */
   long double rule(std::vector<double> const & w, double n)
   {
      long double wMin = w[0], sum = 0.0L;
      for (size_t i = 1; i < w.size(); ++i) {
         if (w[i] < wMin) wMin = w[i];
      }
      for (size_t i = 0; i < w.size(); ++i) {
         sum += std::exp(-(long double)n*((long double)w[i] - wMin));
      }
      return wMin - std::log(sum/(long double)w.size())/(long double)n;
   }

   /*
   * Gaussian bond on a periodic 1D mesh (Fourier symbol exp(-k^2/6)).
   */
   Row bond(Row const & input, double length)
   {
      const int n = (int) input.size();
      const long double pi = std::acos(-1.0L);
      Row result(n, 0.0L);
      for (int i = 0; i < n; ++i) {
         for (int j = 0; j < n; ++j) {
            long double kernel = 0.0L;
            for (int m = 0; m < n; ++m) {
               const int k = m <= n/2 ? m : m - n;
               const long double wave = 2*pi*k/length;
               kernel += std::exp(-wave*wave/6)*std::cos(2*pi*k*(i-j)/n)/n;
            }
            result[i] += kernel*input[j];
         }
      }
      return result;
   }

   /*
   * Unshifted bead propagators of a homopolymer of count beads:
   * forward[j] holds j+1 bead weights and j bonds.
   */
   std::vector<Row> forward(std::vector<double> const & w, int count,
                            double length)
   {
      const int n = (int) w.size();
      std::vector<Row> q(count, Row(n));
      for (int i = 0; i < n; ++i) q[0][i] = std::exp(-(long double)w[i]);
      for (int j = 1; j < count; ++j) {
         q[j] = bond(q[j-1], length);
         for (int i = 0; i < n; ++i) q[j][i] *= std::exp(-(long double)w[i]);
      }
      return q;
   }

   void tearDownModel()
   {  PolymerModel::setModel(PolymerModel::Thread); }

   /*
   * Check the shift of a polymer: the documented one, or, if raised, a
   * larger one for which the shifted partition function is one.
   */
   void checkShift(Polymer<1,T> const & polymer, std::vector<double> const & w,
                   double n, bool raised)
   {
      const double s = (double) rule(w, n);
      if (raised) {
         TEST_ASSERT(polymer.block(0).wShift() > s);
         TEST_ASSERT(std::abs(std::log(polymer.qShifted())) < 1e-6);
      } else {
         TEST_ASSERT(std::abs(polymer.block(0).wShift() - s)
                     < 1e-10*(1.0 + std::abs(s)));
      }
   }

   /*
   * Two thread homopolymers (lengths length and length1) in a given
   * field, compared with reference values for mu, c0 at grid point 0,
   * f and the stress.
   */
   void checkThread(double length, double length1, double ds, int mesh,
                    double cell, std::vector<double> const & w0,
                    std::vector<double> const & w1,
                    double mu, double c0, double f, double stress,
                    bool raised0, bool raised1)
   {
      System<1,T> system;
      makeSystem(system, false, length, ds, mesh, cell, false, length1);
      setFields(system, w0, w1);
      bool threw = false;
      try {
         system.compute();
      } catch (Util::Exception& e) {
         threw = true;
      }
      TEST_ASSERT(!threw);
      Polymer<1,T> const & p0 = system.mixture().polymer(0);
      Polymer<1,T> const & p1 = system.mixture().polymer(1);

      TEST_ASSERT(std::isfinite(p0.mu()));
      TEST_ASSERT(std::abs(p0.mu() - mu) < 1e-9*(1.0 + std::abs(mu)));
      Util::DArray<typename T::Real> host;
      host.allocate(mesh);
      VecOp::eqV(host, system.c().rgrid(0));
      for (int i = 0; i < mesh; ++i) TEST_ASSERT(std::isfinite(host[i]));
      TEST_ASSERT(std::abs(host[0] - c0) < 1e-9);
      system.computeStress();
      TEST_ASSERT(std::isfinite(system.mixture().stress(0)));
      TEST_ASSERT(std::abs(system.mixture().stress(0) - stress)
                  < 1e-8*std::abs(stress));
      system.scft().compute();
      TEST_ASSERT(std::isfinite(system.scft().fHelmholtz()));
      TEST_ASSERT(std::abs(system.scft().fHelmholtz() - f) < 1e-9);

      checkShift(p0, w0, length, raised0);
      checkShift(p1, w1, length1, raised1);
   }

   /*
   * A bead homopolymer in a solvent, compared with the transfer matrix.
   */
   void checkBead(int count, double cell, std::vector<double> const & w0,
                  bool raised)
   {
      const int mesh = (int) w0.size();
      std::vector<double> w1(mesh, 0.0);
      System<1,T> system;
      makeSystem(system, true, (double) count, 0.0, mesh, cell, true);
      setFields(system, w0, w1);
      bool threw = false;
      try {
         system.compute();
      } catch (Util::Exception& e) {
         threw = true;
      }
      TEST_ASSERT(!threw);
      Polymer<1,T> const & polymer = system.mixture().polymer(0);

      // Exact ln q and concentration of the polymer
      std::vector<Row> q = forward(w0, count, cell);
      long double qTot = 0.0L;
      for (int i = 0; i < mesh; ++i) qTot += q[count-1][i]/mesh;
      Row c(mesh, 0.0L);
      for (int j = 0; j < count; ++j) {
         for (int i = 0; i < mesh; ++i) {
            c[i] += q[j][i]*q[count-1-j][i]
                    *std::exp((long double)w0[i])/(count*qTot);
         }
      }
      TEST_ASSERT(std::isfinite(polymer.lnQ()));
      TEST_ASSERT(std::abs(polymer.lnQ() - (double) std::log(qTot))
                  < 1e-10*std::abs((double) std::log(qTot)));
      Util::DArray<typename T::Real> host;
      host.allocate(mesh);
      VecOp::eqV(host, system.c().rgrid(0));
      for (int i = 0; i < mesh; ++i) {
         TEST_ASSERT(std::isfinite(host[i]));
         TEST_ASSERT(std::abs(host[i] - 0.5*(double)c[i]) < 1e-9);
      }

      checkShift(polymer, w0, count, raised);
   }

   /*
   * A thread A-B diblock (lengths lA and lB, phi 1, chi 0, ds 1) in given
   * fields on 16 points (cell 4), compared with reference values for mu,
   * f and the concentrations at grid point 0. The shift of the A block,
   * which is in a uniform field, must not be raised; the B block shift is
   * raised so that Q' = 1.
   */
   void checkDiblock(double lA, double lB, std::vector<double> const & w0,
                     std::vector<double> const & w1, double mu, double f,
                     double c0, double c1)
   {
      std::ostringstream text;
      text << std::setprecision(17)
           << "System{ Mixture{ nMonomer 2 monomers[ 1 1 ] nPolymer 1 "
           << "Polymer{ type linear nBlock 2 blocks[ 0 " << lA << " 1 " << lB
           << " ] phi 1 } ds 1 } Interaction{ chi( 0 1 0 ) } "
           << "Domain{ mesh 16 lattice lamellar } }";
      System<1,T> system;
      std::istringstream in(text.str());
      system.readParam(in);
      Util::FSArray<double,6> cell;
      cell.append(4.0);
      system.setUnitCell(cell);
      setFields(system, w0, w1);
      bool threw = false;
      try {
         system.compute();
      } catch (Util::Exception& e) {
         threw = true;
      }
      TEST_ASSERT(!threw);
      Polymer<1,T> const & polymer = system.mixture().polymer(0);
      TEST_ASSERT(std::isfinite(polymer.mu()));
      TEST_ASSERT(std::abs(polymer.mu() - mu) < 1e-9*(1.0 + std::abs(mu)));
      Util::DArray<typename T::Real> host;
      host.allocate(16);
      VecOp::eqV(host, system.c().rgrid(0));
      TEST_ASSERT(std::abs(host[0] - c0) < 1e-9);
      VecOp::eqV(host, system.c().rgrid(1));
      TEST_ASSERT(std::abs(host[0] - c1) < 1e-9);
      system.scft().compute();
      TEST_ASSERT(std::abs(system.scft().fHelmholtz() - f) < 1e-9);

      // Only the B block is raised
      const double n = lA + lB;
      const double sA = (double) rule(w0, n), sB = (double) rule(w1, n);
      TEST_ASSERT(std::abs(polymer.block(0).wShift() - sA) < 1e-8);
      TEST_ASSERT(polymer.block(1).wShift() > sB + 0.1);
      TEST_ASSERT(std::abs(std::log(polymer.qShifted())) < 1e-6);
   }

   /*
   * Linear chain of 1000 (beads) in a lamellar cell 2048 (mesh 4096) with
   * w = 1000 except w = 0 at one grid point, as one block or as k blocks
   * of equal length (thread ds 1). The k-block chain is the same discrete
   * operator as the single block, so mu, f and c must agree.
   */
   void computeManyBlocks(bool bead, int k, double & mu, double & f,
                          double & c0)
   {
      const int nGrid = 4096;
      std::ostringstream text;
      text << std::setprecision(17) << "System{ ";
      if (bead) text << "polymerModel bead ";
      text << "Mixture{ nMonomer 2 monomers[ 1 1 ] nPolymer 1 "
           << "Polymer{ type linear nBlock " << k << " blocks[ ";
      for (int j = 0; j < k; ++j) {
         if (bead) {
            text << "0 " << 1000/k << " ";
         } else {
            text << "0 " << 1000.0/k << " ";
         }
      }
      text << "] phi 1 } ";
      if (!bead) text << "ds 1 ";
      text << "} Interaction{ chi( 0 1 0 ) } "
           << "Domain{ mesh " << nGrid << " lattice lamellar } }";
      System<1,T> system;
      std::istringstream in(text.str());
      system.readParam(in);
      Util::FSArray<double,6> cell;
      cell.append(2048.0);
      system.setUnitCell(cell);
      std::vector<double> w(nGrid, 1000.0);
      w[0] = 0.0;
      setFields(system, w, w);
      bool threw = false;
      try {
         system.compute();
      } catch (Util::Exception& e) {
         threw = true;
      }
      TEST_ASSERT(!threw);
      mu = system.mixture().polymer(0).mu();
      Util::DArray<typename T::Real> host;
      host.allocate(nGrid);
      VecOp::eqV(host, system.c().rgrid(0));
      c0 = host[0];
      system.scft().compute();
      f = system.scft().fHelmholtz();
      TEST_ASSERT(std::isfinite(mu) && std::isfinite(f) && std::isfinite(c0));
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
      tearDownModel();
      Util::BracketPolicy::set(bracketPolicy_);
   }

   /*
   * Thread chains in strong nonuniform fields, which overflowed with
   * the shift by the average (NaN with exit 0) [v01].
   */
   void checkThreadStrongFields()
   {
      openLogFile("out/ScftFieldShiftThread.log");

      // Step fields 0|D on 8 points, length 1, ds 0.01, cell 4
      double D[2] = {1500.0, 3000.0};
      double mu[2] = {3.5177556860557e-01, 3.5186929242297e-01};
      double c0[2] = {6.8432160663712e-01, 6.8425322869696e-01};
      double f[2] = {-6.4870896522215e-01, -6.4815322200917e-01};
      double stress[2] = {-1.51795723251e-01, -1.51837585834e-01};
      for (int k = 0; k < 2; ++k) {
         std::vector<double> w0(8), w1(8);
         for (int i = 0; i < 8; ++i) {
            w0[i] = i < 4 ? 0.0 : D[k];
            w1[i] = D[k] - w0[i];
         }
         checkThread(1.0, 1.0, 0.01, 8, 4.0, w0, w1, mu[k], c0[k], f[k],
                     stress[k], false, false);
      }

      // Long chains in (0,0,2,2) on 4 points, ds 0.1, cell 50
      double n[2] = {705.0, 1000.0};
      double muL[2] = {9.2581730256342e-01, 1.3137663138651e+00};
      double c0L[2] = {9.9999863499693e-01, 9.9999890995902e-01};
      double fL[2] = {-1.0795369041722e-04, 3.1158623185411e-04};
      double stressL[2] = {-5.24940555301e-05, -5.25160525047e-05};
      std::vector<double> w0 = {0.0, 0.0, 2.0, 2.0};
      std::vector<double> w1 = {2.0, 2.0, 0.0, 0.0};
      for (int k = 0; k < 2; ++k) {
         checkThread(n[k], n[k], 0.1, 4, 50.0, w0, w1, muL[k], c0L[k], fL[k],
                     stressL[k], false, false);
      }
      closeLogFile();
   }

   /*
   * Bead propagator slices follow the documented shift for every chain
   * length; before, they switched from the average to the minimum [v04].
   */
   void checkBeadSlices()
   {
      openLogFile("out/ScftFieldShiftSlices.log");
      int counts[2] = {150, 200};
      double offsets[2] = {0.0, 7.0};
      for (int a = 0; a < 2; ++a) {
         for (int b = 0; b < 2; ++b) {
            const int count = counts[a];
            const double off = offsets[b];
            System<1,T> system;
            makeSystem(system, true, (double) count, 0.0, 4, 20.0, false);
            std::vector<double> w0 = {off, off, 4.0 + off, 4.0 + off};
            std::vector<double> w1 = {4.0 + off, 4.0 + off, off, off};
            setFields(system, w0, w1);
            system.compute();
            Block<1,T> const & block = system.mixture().polymer(0).block(0);
            const long double s = rule(w0, count);
            TEST_ASSERT(std::abs(block.wShift() - (double) s) < 1e-12*(1.0 + off));

            // Slice S holds S bead weights: shifted = exact*exp(S*s)
            const int slice = count/2;
            std::vector<Row> q = forward(w0, slice, 20.0);
            Util::DArray<typename T::Real> host;
            host.allocate(4);
            VecOp::eqV(host, block.propagator(0).q(slice));
            for (int i = 0; i < 4; ++i) {
               TEST_ASSERT(host[i] > 0.0);
               const long double ratio = std::log((long double)host[i])
                                       - std::log(q[slice-1][i]);
               TEST_ASSERT(std::abs((double)(ratio - slice*s)) < 1e-9);
            }
         }
      }
      closeLogFile();
   }

   /*
   * Long chains in fields that vary on the grid scale, for which the
   * documented shift makes the partition function Q' underflow or come
   * close to it: thread length 1000 in a field alternating 0|2 (ln Q'
   * about -810), and thread lengths 640 and 660 in a field alternating
   * 0|2 on half of the cell and 4 on the other half (ln Q' about -679
   * and -700; NaN before, because the average shift then overflowed).
   * The shifts are raised so that Q' = 1; in a diblock, only the block
   * that needs it, and for many blocks confined to a small part of the
   * cell by a common constant. Bead chains: 1000 beads in the
   * alternating field 0|2 (ln Q' about -566, kept; rejected before) and
   * 2000 beads in 0|1 (ln Q' about -760, raised).
   */
   void checkRaisedShift()
   {
      openLogFile("out/ScftFieldShiftRaised.log");

      // Diblock 1000 + 1000 with w_A = 0 and w_B alternating 0|6, and
      // with w_B = 0|6 on half of the cell and 3000 on the other half.
      // A raise common to both blocks multiplied the A propagators by up
      // to exp(780) (an exception, or NaN for the second field). The
      // references come from a numpy RQM4 solver in the unshifted fields
      // that rescales every slice.
      std::vector<double> zero(16, 0.0), wB(16), wH(16);
      for (int i = 0; i < 16; ++i) {
         wB[i] = i % 2 ? 6.0 : 0.0;
         wH[i] = i < 8 ? wB[i] : 3000.0;
      }
      checkDiblock(1000.0, 1000.0, zero, wB, 1.5671960390073e+03,
                   4.3171021096623e-01, 5.0007769668120e-01,
                   8.8287073050384e-01);
      checkDiblock(1000.0, 1000.0, zero, wH, 1.8397430331006e+03,
                   5.6462601056567e-01, 5.0063344331890e-01,
                   8.5934911461644e-01);

      // Second polymer of length 1 in w = 0; references from devel
      std::vector<double> v0(16), v1(16, 0.0);
      for (int i = 0; i < 16; ++i) {
         v0[i] = i < 8 ? (i % 2 ? 2.0 : 0.0) : 4.0;
      }
      double n[2] = {640.0, 660.0};
      double mu[2] = {6.79731866042e+02, 7.00975750196e+02};
      double c0[2] = {6.1186501890893e-01, 6.1185003418158e-01};
      double f[2] = {-6.77531105762e-01, -6.77499223438e-01};
      double stress[2] = {-2.00714953809e-01, -2.00721147560e-01};
      for (int k = 0; k < 2; ++k) {
         checkThread(n[k], 1.0, 1.0, 16, 4.0, v0, v1, mu[k], c0[k], f[k],
                     stress[k], true, false);
      }

      std::vector<double> w0(8), w1(8);
      for (int i = 0; i < 8; ++i) {
         w0[i] = i % 2 ? 2.0 : 0.0;
         w1[i] = 2.0 - w0[i];
      }
      checkThread(1000.0, 1000.0, 1.0, 8, 2.0, w0, w1, 8.10770127762e+02,
                  6.586796404198e-01, 1.27129408601e-01,
                  -6.80273911492e-01, true, true);

      // 100 blocks of 10 confined to one grid point of 4096: calibrating
      // each block from a uniform head over-raised it by up to ln(4096)
      // per junction (an exception); the common raise then makes the
      // chain the same discrete operator as one block.
      for (int b = 0; b < 2; ++b) {
         double mu1, f1, c1, muK, fK, cK;
         computeManyBlocks(b == 1, 1, mu1, f1, c1);
         computeManyBlocks(b == 1, 100, muK, fK, cK);
         TEST_ASSERT(std::abs(muK - mu1) < 1e-10*std::abs(mu1));
         TEST_ASSERT(std::abs(fK - f1) < 1e-10*std::abs(f1));
         TEST_ASSERT(std::abs(cK - c1) < 1e-10*std::abs(c1));
      }
      tearDownModel();

      checkBead(1000, 2.0, w0, false);
      std::vector<double> weak(8);
      for (int i = 0; i < 8; ++i) weak[i] = i % 2 ? 1.0 : 0.0;
      checkBead(2000, 2.0, weak, true);
      closeLogFile();
   }

};

}
}
#endif
