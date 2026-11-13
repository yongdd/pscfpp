#ifndef PSCF_WAVE_NORM_H
#define PSCF_WAVE_NORM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/Sort.h>
#include <util/global.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace Pscf {

   /** Relative precision used to compare squared wavevector norms. */
   template <typename Real>
   inline double waveNormTolerance()
   {
      return std::max(1.0E-12,
                      64.0*double(std::numeric_limits<Real>::epsilon()));
   }

   /**
   * Compare nonnegative squared norms independently of length units.
   *
   * Zero belongs only to the zero shell. Comparing to the larger norm
   * avoids a dimensionful absolute floor that would merge small modes.
   */
   inline bool equalWaveNorm(double a, double b,
                            double tolerance = waveNormTolerance<double>())
   {
      if (a == 0.0 || b == 0.0) return a == b;
      return std::abs(a-b) <= tolerance*std::max(a,b);
   }

   /** Group sorted wavevector norms using the same relative criterion. */
   inline void findWaveBunches(std::vector< Sort::Item<double> > const & items,
                              Util::GArray< Sort::Bunch >& bunches,
                              double tolerance)
   {
      UTIL_CHECK(!items.empty());
      bunches.clear();
      Sort::Bunch bunch;
      bunch[0] = 0;
      for (int i = 0; i < int(items.size()); ++i) {
         UTIL_CHECK(std::isfinite(items[i].value) && items[i].value >= 0.0);
         if (i == 0) continue;
         UTIL_CHECK(items[i].value >= items[i-1].value);
         if (!equalWaveNorm(items[i].value, items[bunch[0]].value,
                            tolerance)) {
            bunch[1] = i;
            bunches.append(bunch);
            bunch[0] = i;
         }
      }
      bunch[1] = items.size();
      bunches.append(bunch);
   }

}
#endif
