#ifndef RP_POLYMER_TPP
#define RP_POLYMER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include <rp/solvers/Block.h>
#include <rp/solvers/Propagator.h>
#include <pscf/field/RField.h>
#include <pscf/backend/Reduce.h>
#include <pscf/backend/VecOp.h>
#include <pscf/solvers/PolymerTmpl.tpp>
#include <pscf/chem/PolymerModel.h>
#include <pscf/iterator/NanException.h>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <vector>

namespace Pscf {
namespace Rp {

   /*
   * Constructor.
   */
   template <int D, class T>
   Polymer<D,T>::Polymer()
    : stress_(),
      nParam_(0)
   {  ParamComposite::setClassName("Polymer"); }

   /*
   * Destructor.
   */
   template <int D, class T>
   Polymer<D,T>::~Polymer()
   {}

   /*
   * Set the number of unit cell parameters.
   */
   template <int D, class T>
   void Polymer<D,T>::setNParams(int nParam)
   {  nParam_ = nParam; }

   /*
   * Clear all data that depends on unit cell parameters.
   */
   template <int D, class T>
   void Polymer<D,T>::clearUnitCellData()
   {
      for (int j = 0; j < nBlock(); ++j) {
         block(j).clearUnitCellData();
      }
      stress_.clear();
   }

   /*
   * Compute solution to MDE and block concentrations.
   */
   template <int D, class T>
   void Polymer<D,T>::compute(DArray< RField<D,T> > const & wFields,
                              double phiTot)
   {
      // Field shifts (see the documentation of this function). Block j
      // solves the MDE in w - s[j], where w is the field of its monomer
      // type and s[j] = min(w) - ln<exp(-n(w - min(w)))>/n is the free
      // energy per monomer of a whole molecule of n monomers held at one
      // point. Then the propagators stay below about the number of grid
      // points and qShifted() is at most about one, in any field; it is
      // raised below if qShifted() is too small.
      const bool isBead = PolymerModel::isBead();
      const int nb = nBlock();
      const double n = isBead ? (double) nBead() : length();
      std::vector<double> shifts(nb), minima(nb), maxima(nb);
      for (int j = 0; j < nb; ++j) {
         RField<D,T> const & w = wFields[block(j).monomerId()];
         if (wWork_.isAllocated() && wWork_.capacity() != w.capacity()) {
            wWork_.deallocate();
         }
         if (!wWork_.isAllocated()) {
            wWork_.allocate(w.meshDimensions());
         }
         minima[j] = Reduce::min(w);
         if (isBead) maxima[j] = Reduce::max(w);
         VecOp::subVS(wWork_, w, minima[j]);
         VecOp::expVc(wWork_, wWork_, -n);
         shifts[j] = minima[j]
                   - std::log(Reduce::sum(wWork_)/double(w.capacity()))/n;
      }

      const double logLimit = std::log(std::numeric_limits<typename T::Real>::max()) - 32.0;
      // Leave room for Fourier products, contour sums and normalization.
      const double qLimit = std::exp(0.5*logLimit - 16.0)
                           /wFields[block(0).monomerId()].capacity();
      // Smallest accepted qShifted(): products of the shifted propagators
      // are of the order of qShifted() and are normalized by it.
      const double qMin = std::exp(-logLimit);

      // Set up all blocks with given shifts. For bead chains, bound
      // amplification over the entire molecule (including branch sources),
      // not just the exponential for one bead. The direct concentration
      // formula contains two propagators and one reciprocal bead weight.
      auto setupBlocks = [&](std::vector<double> const & s) -> bool {
         double amplification = 0.0, inverseRange = 0.0;
         if (isBead) {
            for (int j = 0; j < nb; ++j) {
               amplification += block(j).nBead()*(s[j] - minima[j]);
               inverseRange = std::max(inverseRange, maxima[j] - s[j]);
            }
         }
         const bool safe = isBead
                        && !(2.0*amplification + inverseRange < logLimit);
         for (int j = 0; j < nb; ++j) {
            block(j).setupSolver(wFields[block(j).monomerId()], s[j], safe);
         }
         return safe;
      };
      bool safeBead = setupBlocks(shifts);

      // Call base class PolymerTmpl solve() function
      // This solves the MDE for all propagators in a precalculated order
      PolymerTmplT::solve(phiTot);

      // Fields that are not finite give a shift that is not finite
      bool finiteFields = true;
      for (int j = 0; j < nb; ++j) {
         if (!std::isfinite(shifts[j])) finiteFields = false;
      }

      // Raise shifts until an average q of shifted propagators is one,
      // given a function that adds a constant to the shifts concerned and
      // returns the new q, and the length (or number of beads) l over which
      // the constant applies: while q < qMin (q may then be zero or an
      // inaccurate subnormal number), in steps of logLimit/l, each of which
      // multiplies q by exp(logLimit) and so cannot make it overflow, and
      // which end before the shifts exceed max(w) (range = max(w) - shift);
      // then by -ln(q)/l.
      auto raiseToOne = [&](double q, double l, double range,
                            std::function<double(double)> const & raise) {
         const double nStep = std::min(std::ceil(l*range/logLimit) + 1.0,
                                       1000.0);
         for (int k = 0; k < nStep && q < qMin; ++k) {
            q = raise(logLimit/l);
         }
         if (q >= qMin && std::isfinite(q)) {
            raise(-std::log(q)/l);
         }
      };

      // If qShifted() is below qMin (a long chain in a field that varies
      // strongly on scales it cannot follow), raise the shifts so that
      // qShifted() becomes one. A raise common to all blocks would also
      // multiply the propagators of blocks that need none (for example a
      // block in a uniform field) by up to 1/qShifted() and could make
      // them overflow, so each block is first calibrated alone: its
      // propagator is solved from a uniform head, and its shift is raised
      // until the average of the tail is one. The polymer is then solved
      // with all shifts raised by a common, usually small, constant so
      // that qShifted() is one.
      if (qShifted() < qMin && finiteFields) {
         if (!isBead) {
            for (int j = 0; j < nb; ++j) {
               maxima[j] = Reduce::max(wFields[block(j).monomerId()]);
            }
         }
         const std::vector<double> documented(shifts);
         VecOp::eqS(wWork_, 1.0);
         for (int j = 0; j < nb; ++j) {
            RField<D,T> const & w = wFields[block(j).monomerId()];
            PropagatorT& p = block(j).propagator(0);
            auto raiseBlock = [&](double delta) -> double {
               shifts[j] += delta;
               block(j).setupSolver(w, shifts[j], true);
               p.solve(wWork_);
               return Reduce::sum(p.q(p.ns()-1))/double(wWork_.capacity());
            };
            const double l = isBead ? (double) block(j).nBead()
                                    : block(j).length();
            raiseToOne(raiseBlock(0.0), l, maxima[j] - shifts[j], raiseBlock);
         }
         auto raisePolymer = [&](double delta) -> double {
            for (int j = 0; j < nb; ++j) shifts[j] += delta;
            safeBead = setupBlocks(shifts);
            PolymerTmplT::solve(phiTot);
            return qShifted();
         };
         auto raiseCommon = [&]() {
            double range = 0.0;
            for (int j = 0; j < nb; ++j) {
               range = std::max(range, maxima[j] - shifts[j]);
            }
            raiseToOne(raisePolymer(0.0), n, range, raisePolymer);
         };
         raiseCommon();

         // A block calibrated from a uniform head is raised too much when
         // its head in the polymer is concentrated where the chain sits,
         // by up to ln(1/f) per junction (f the occupied volume
         // fraction), which can overflow for many blocks. Then raise the
         // documented shifts by a common constant only: equal shifts make
         // such a chain the same discrete operator as a single block.
         if (!(qShifted() >= qMin && std::isfinite(qShifted())
               && std::isfinite(PolymerTmplT::lnQ())
               && std::isfinite(phi())
               && (!isBead || qShifted() <= qLimit))) {
            shifts = documented;
            raiseCommon();
         }
      }

      // Reject results that are still not representable. Thread polymers
      // in fields that are not finite give NaN values, as before, which
      // the iterators report as a NaN error.
      bool valid = qShifted() >= qMin && std::isfinite(qShifted())
                   && std::isfinite(PolymerTmplT::lnQ())
                   && std::isfinite(phi());
      if (isBead && !(qShifted() <= qLimit)) valid = false;
      if (!valid && (isBead || finiteFields)) {
         // A NanException: iterators treat it as a failed solve
         throw NanException(UTIL_FUNC,
               "Polymer normalization exceeds floating-point range",
               __FILE__, __LINE__);
      }

      // Compute block concentration fields, normalizing products of
      // the shifted propagators by their own partition function
      double prefactor;
      if (PolymerModel::isThread()) {
         prefactor = phi() / ( qShifted() * length() );
         for (int i = 0; i < nBlock(); ++i) {
            block(i).computeConcentrationThread(prefactor);
            Block<D,T> const& current = block(i);
            if (finiteFields && !std::isfinite(Reduce::sum(current.cField()))) {
               throw NanException(UTIL_FUNC,
                     "Polymer concentration exceeds floating-point range",
                     __FILE__, __LINE__);
            }
         }
      } else {
         prefactor = phi() / qShifted() / (double)nBead();
         for (int i = 0; i < nBlock(); ++i) {
            if (safeBead) {
               const double logPrefactor = std::log(phi()) - std::log(qShifted())
                                         -std::log((double)nBead());
               block(i).computeConcentrationBeadLog(logPrefactor);
            } else {
               block(i).computeConcentrationBead(prefactor);
            }
            Block<D,T> const& current = block(i);
            if (!std::isfinite(Reduce::sum(current.cField()))) {
               throw NanException(UTIL_FUNC,
                     "Bead polymer concentration exceeds floating-point range",
                     __FILE__, __LINE__);
            }
         }
      }

   }

   /*
   * Compute stress contribution from a polymer species.
   */
   template <int D, class T>
   void Polymer<D,T>::computeStress()
   {
      UTIL_CHECK(nParam_ > 0);

      // Compute stress contributions for all blocks
      double prefactor;
      if (PolymerModel::isThread()) {
         prefactor = phi() / ( qShifted() * length() );
         for (int i = 0; i < nBlock(); ++i) {
            block(i).computeStressThread(prefactor);
         }
      } else {
         prefactor = phi() / qShifted() / (double)nBead();
         if (!std::isfinite(prefactor)) {
            throw NanException(UTIL_FUNC,
               "Bead polymer stress normalization exceeds floating-point range",
               __FILE__, __LINE__);
         }
         for (int i = 0; i < nBlock(); ++i) {
            block(i).computeStressBead(prefactor);
         }
      }

      // Initialize all stress_ elements to zero
      stress_.clear();
      for (int i = 0; i < nParam_; ++i) {
         stress_.append(0.0);
      }

      // Sum over all block stress contributions
      for (int i = 0; i < nBlock(); ++i) {
         for (int j = 0; j < nParam_; ++j){
            stress_[j] += block(i).stress(j);
            if (PolymerModel::isBead() && !std::isfinite(stress_[j])) {
               throw NanException(UTIL_FUNC,
                     "Bead polymer stress exceeds floating-point range",
                     __FILE__, __LINE__);
            }
         }
      }

   }

}
}
#endif
