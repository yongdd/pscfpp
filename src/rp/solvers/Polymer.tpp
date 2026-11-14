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
#include <pscf/solvers/PolymerTmpl.tpp>
#include <pscf/chem/PolymerModel.h>
#include <pscf/iterator/NanException.h>
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
      // Retain the mean shift for ordinary fields. For bead chains, bound
      // amplification over the entire molecule (including branch sources),
      // not just the exponential for one bead. The direct concentration
      // formula contains two propagators and one reciprocal bead weight.
      std::vector<double> means(nBlock()), minima(nBlock());
      double amplification = 0.0, inverseRange = 0.0, largestStep = 0.0;
      for (int j = 0; j < nBlock(); ++j) {
         const int monomerId = block(j).monomerId();
         RField<D,T> const & w = wFields[monomerId];
         means[j] = Reduce::sum(w) / double(w.capacity());
         if (PolymerModel::isBead()) {
            minima[j] = Reduce::min(w);
            const double step = means[j] - minima[j];
            amplification += block(j).nBead()*step;
            largestStep = std::max(largestStep, step);
            inverseRange = std::max(inverseRange, Reduce::max(w)-means[j]);
         }
      }
      const double logLimit = std::log(std::numeric_limits<typename T::Real>::max()) - 32.0;
      // Leave room for Fourier products, contour sums and normalization.
      const double qLimit = std::exp(0.5*logLimit - 16.0)
                           /wFields[block(0).monomerId()].capacity();
      bool safeBead = PolymerModel::isBead()
                  && !(2.0*amplification + inverseRange < logLimit);
      bool useMinimum = safeBead && !(largestStep < logLimit);
      for (int j = 0; j < nBlock(); ++j) {
         block(j).setupSolver(wFields[block(j).monomerId()],
                              useMinimum ? minima[j] : means[j], safeBead);
      }

      // Call base class PolymerTmpl solve() function
      // This solves the MDE for all propagators in a precalculated order
      PolymerTmplT::solve(phiTot);

      // A mean-shifted Q can overflow even when the physical Q and all
      // normalized observables are representable. Retry with bounded local
      // weights before publishing concentrations. The raw q-slice convention
      // is unchanged: these are still propagators in w - block.wShift().
      if (PolymerModel::isBead()) {
         if (!(qShifted() > 0.0) || !std::isfinite(qShifted())
             || qShifted() > qLimit) {
            if (!useMinimum) {
               safeBead = true;
               for (int j = 0; j < nBlock(); ++j) {
                  block(j).setupSolver(wFields[block(j).monomerId()], minima[j], true);
               }
               PolymerTmplT::solve(phiTot);
            }
         }
         if (!(qShifted() >= 1.0/qLimit && qShifted() <= qLimit)
             || !std::isfinite(PolymerTmplT::lnQ()) || !std::isfinite(phi())) {
            // A NanException: iterators treat it as a failed solve
            throw NanException(UTIL_FUNC,
                  "Bead polymer normalization exceeds floating-point range",
                  __FILE__, __LINE__);
         }
      }

      // Compute block concentration fields, normalizing products of
      // the shifted propagators by their own partition function
      double prefactor;
      if (PolymerModel::isThread()) {
         prefactor = phi() / ( qShifted() * length() );
         for (int i = 0; i < nBlock(); ++i) {
            block(i).computeConcentrationThread(prefactor);
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
