/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include "Block.h"
#include "Propagator.h"
#include <r1/field/Domain.h>
#include <pscf/solvers/PolymerTmpl.tpp>

#include <pscf/iterator/NanException.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

// Explicit instantiation of base class
namespace Pscf {
   template class PolymerTmpl<R1::Block, R1::Propagator, double>;
}

namespace Pscf {
namespace R1 { 

   /*
   * Constructor
   */
   Polymer::Polymer()
   {  setClassName("Polymer"); }

   /*
   * Destructor
   */
   Polymer::~Polymer()
   {}

   /*
   * Compute solution to MDE and concentrations.
   */ 
   void Polymer::compute(DArray< DArray<double> > const & wFields)
   {

      // Setup solvers for all blocks. Each block solves the MDE in the
      // field w of its monomer type minus the shift
      // min(w) - ln<exp(-N(w - min(w)))>/N, N = length() (see Polymer.h),
      // which PolymerTmpl::solve removes from ln q analytically. This
      // keeps the propagators representable, makes the results
      // independent of the constant of w, and keeps the shifted field
      // small where the molecules are.
      const double n = length();
      const int nb = nBlock();
      std::vector<double> shifts(nb), maxima(nb);
      for (int j = 0; j < nb; ++j) {
         DArray<double> const & w = wFields[block(j).monomerId()];
         const int nx = w.capacity();
         if (wWork_.isAllocated() && wWork_.capacity() != nx) {
            wWork_.deallocate();
         }
         if (!wWork_.isAllocated()) {
            wWork_.allocate(nx);
         }
         double wMin = w[0], wMax = w[0];
         for (int i = 1; i < nx; ++i) {
            if (w[i] < wMin) wMin = w[i];
            if (w[i] > wMax) wMax = w[i];
         }
         for (int i = 0; i < nx; ++i) {
            wWork_[i] = std::exp(-n*(w[i] - wMin));
         }
         shifts[j] = wMin - std::log(block(j).domain().spatialAverage(wWork_))/n;
         maxima[j] = wMax;

         // A field that is not finite gives a shift that is not finite
         // (NaN or -inf at any point) or an infinite maximum (+inf).
         // The tridiagonal solver rejects such coefficients, so report
         // it as a NanException, which iterators treat as a failed solve.
         if (!std::isfinite(shifts[j]) || !std::isfinite(wMax)) {
            throw NanException(UTIL_FUNC, "Polymer field w is not finite",
                               __FILE__, __LINE__);
         }
         block(j).setupSolver(w, shifts[j]);
      }

      // Solve MDE for all propagators
      solve();

      // Raise shifts until an average q of shifted propagators is one,
      // given a function that adds a constant to the shifts concerned and
      // returns the new q, and the length l over which the constant
      // applies: while q < qMin (q may then be zero or an inaccurate
      // subnormal number), in steps of logLimit/l, each of which
      // multiplies q by about exp(logLimit) and so cannot make it
      // overflow, and which end before the shifts exceed max(w) (range =
      // max(w) - shift); then by -ln(q)/l.
      const double logLimit
                  = std::log(std::numeric_limits<double>::max()) - 32.0;
      const double qMin = std::exp(-logLimit);
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
      // until the spatial average of the tail is one. The polymer is then
      // solved with all shifts raised by a common, usually small,
      // constant so that qShifted() is one.
      if (qShifted() < qMin) {
         const std::vector<double> documented(shifts);
         for (int i = 0; i < wWork_.capacity(); ++i) {
            wWork_[i] = 1.0;
         }
         for (int j = 0; j < nb; ++j) {
            DArray<double> const & w = wFields[block(j).monomerId()];
            Propagator& p = block(j).propagator(0);
            auto raiseBlock = [&](double delta) -> double {
               shifts[j] += delta;
               block(j).setupSolver(w, shifts[j]);
               p.solve(wWork_);
               return block(j).domain().spatialAverage(p.tail());
            };
            raiseToOne(raiseBlock(0.0), block(j).length(),
                       maxima[j] - shifts[j], raiseBlock);
         }
         auto raisePolymer = [&](double delta) -> double {
            for (int j = 0; j < nb; ++j) {
               shifts[j] += delta;
               block(j).setupSolver(wFields[block(j).monomerId()],
                                    shifts[j]);
            }
            solve();
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
               && std::isfinite(lnQ()) && std::isfinite(phi()))) {
            shifts = documented;
            raiseCommon();
         }
      }

      // Reject results that are still not representable (a NanException,
      // which iterators treat as a failed solve).
      if (!(qShifted() >= qMin) || !std::isfinite(qShifted())
          || !std::isfinite(lnQ()) || !std::isfinite(phi())) {
         throw NanException(UTIL_FUNC,
               "Polymer normalization exceeds floating-point range",
               __FILE__, __LINE__);
      }

      // Compute block concentration fields (thread model), normalizing
      // products of the shifted propagators by their own partition
      // function
      double prefactor;
      prefactor = phi() / ( qShifted() * length() );
      for (int i = 0; i < nb; ++i) {
         block(i).computeConcentration(prefactor);
         Block const & current = block(i);
         DArray<double> const & c = current.cField();
         for (int k = 0; k < c.capacity(); ++k) {
            if (!std::isfinite(c[k])) {
               throw NanException(UTIL_FUNC,
                     "Polymer concentration exceeds floating-point range",
                     __FILE__, __LINE__);
            }
         }
      }

   }

}
}
