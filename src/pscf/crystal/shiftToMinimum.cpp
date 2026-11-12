/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "shiftToMinimum.h"
#include <pscf/math/RealVec.h>
#include <util/containers/FArray.h>
#include <util/global.h>

namespace Pscf {

   using namespace Util;

   namespace {

      /*
      * Scalar product of the reciprocal lattice vectors with integer
      * coordinates a and b.
      */
      template <int D>
      double kDot(IntVec<D> const & a, IntVec<D> const & b,
                  UnitCell<D> const & cell)
      {
         RealVec<D> ga(0.0);
         RealVec<D> gb(0.0);
         RealVec<D> p;
         for (int i = 0; i < D; ++i) {
            p.multiply(cell.kBasis(i), a[i]);
            ga += p;
            p.multiply(cell.kBasis(i), b[i]);
            gb += p;
         }
         double value = 0.0;
         for (int i = 0; i < D; ++i) {
            value += ga[i]*gb[i];
         }
         return value;
      }

      /*
      * Image shift neighbors of the last mesh and cell used by this
      * thread (shiftToMinimum is called once per wavevector).
      */
      template <int D>
      struct NeighborCache
      {
         IntVec<D> d;
         double kBasis[D][D];
         FSArray<IntVec<D>, 14> neighbors;
         bool isValid;

         NeighborCache()
          : isValid(false)
         {}
      };

      template <int D>
      FSArray<IntVec<D>, 14> const &
      cachedNeighbors(IntVec<D> const & d, UnitCell<D> const & cell)
      {
         thread_local NeighborCache<D> cache;
         bool same = cache.isValid && (cache.d == d);
         for (int i = 0; same && i < D; ++i) {
            for (int j = 0; j < D; ++j) {
               if (cache.kBasis[i][j] != cell.kBasis(i)[j]) same = false;
            }
         }
         if (!same) {
            imageShiftNeighbors<D>(d, cell, cache.neighbors);
            cache.d = d;
            for (int i = 0; i < D; ++i) {
               for (int j = 0; j < D; ++j) {
                  cache.kBasis[i][j] = cell.kBasis(i)[j];
               }
            }
            cache.isValid = true;
         }
         return cache.neighbors;
      }

      /*
      * Move image r to a minimum image.
      *
      * Moves r by the image shift neighbor that lowers |k|^2 most, for
      * as long as one lowers it by more than epsilon. The image is then
      * moved to the largest (in the order of IntVec) of the neighboring
      * images with equal |k|^2, repeatedly. The neighbors contain every
      * difference between tied minimum images (they are the Voronoi-
      * relevant vectors), so this reaches the largest tied image from
      * any starting image: the result depends only on the class of v,
      * not on the representative passed in (Basis relies on this when
      * it compares the images of v and -v). Where the window search
      * already found the largest tied image, nothing changes.
      * The CUDA WaveList kernels implement the same steps.
      */
      template <int D>
      void descendToMinimum(IntVec<D>& r, IntVec<D> const & d,
                            UnitCell<D> const & cell, double epsilon)
      {
         FSArray<IntVec<D>, 14> const & nb = cachedNeighbors<D>(d, cell);
         const int n = nb.size();
         double Gsq = cell.ksq(r);
         IntVec<D> u;
         for (;;) {
            int jBest = -1;
            double GsqBest = Gsq - epsilon;
            for (int j = 0; j < n; ++j) {
               u.add(r, nb[j]);
               double g = cell.ksq(u);
               if (g < GsqBest) {
                  jBest = j;
                  GsqBest = g;
               }
            }
            if (jBest < 0) break;
            r += nb[jBest];
            Gsq = GsqBest;
         }
         bool changed = true;
         while (changed) {
            changed = false;
            for (int j = 0; j < n; ++j) {
               u.add(r, nb[j]);
               if ((cell.ksq(u) < Gsq + epsilon) && (r < u)) {
                  r = u;
                  changed = true;
               }
            }
         }
      }

   }

   template <>
   IntVec<1> shiftToMinimum(IntVec<1> const & v,
                            IntVec<1> const & d,
                            UnitCell<1> const & cell)
   {
      IntVec<1> u;
      if (v[0] > d[0]/2) {
         u[0] = v[0] - d[0];
      } else
      if (v[0] < -(d[0]/2) ) {
         u[0] = v[0] + d[0];
      } else {
         u[0] = v[0];
      }
      UTIL_ASSERT(abs(u[0]) <= d[0]/2);
      return u;
   }

   template <>
   IntVec<2> shiftToMinimum(IntVec<2> const & v,
                            IntVec<2> const & d,
                            UnitCell<2> const & cell)
   {
      // Initialize minimum to input value
      const double epsilon = 1.0E-8;
      double Gsq = cell.ksq(v);
      double Gsq_min_lo = Gsq - epsilon;
      double Gsq_min_hi = Gsq + epsilon;

      // Loop over neighboring images
      IntVec<2> r = v;               // Minimum image of v
      IntVec<2> u;                   // Vector used in loop
      int s0, s1;
      const int q = 2;
      for (s0 = -q; s0 < q+1; ++s0) {
         u[0] = v[0] + s0*d[0];
         for (s1 = -q; s1 < q+1; ++s1) {
            u[1] = v[1] + s1*d[1];
            Gsq = cell.ksq(u);
            if (Gsq < Gsq_min_hi) {
               if ((Gsq < Gsq_min_lo) || (r < u)) {
                  Gsq_min_lo = Gsq - epsilon;
                  Gsq_min_hi = Gsq + epsilon;
                  r = u;
               }
            }
         }
      }

      // Continue from r if the window missed the minimum image
      descendToMinimum<2>(r, d, cell, epsilon);
      return r;
   }

   template <>
   IntVec<3> shiftToMinimum(IntVec<3> const & v,
                            IntVec<3> const & d,
                            UnitCell<3> const & cell)
   {
      // Initialize minimum to input value
      const double epsilon = 1.0E-8;
      double Gsq = cell.ksq(v);
      double Gsq_min_lo = Gsq - epsilon;
      double Gsq_min_hi = Gsq + epsilon;

      // Loop over neighboring images
      IntVec<3> r = v;               // Minimum image
      IntVec<3> u;                   // Vector used in loop
      int s0, s1, s2;
      const int q = 2;
      for (s0 = -q; s0 < q + 1; ++s0) {
         u[0] = v[0] + s0*d[0];
         for (s1 = -q; s1 < q + 1; ++s1) {
            u[1] = v[1] + s1*d[1];
            for (s2 = -q; s2 < q + 1; ++s2) {
               u[2] = v[2] + s2*d[2];
               Gsq = cell.ksq(u);
               if (Gsq < Gsq_min_hi) {
                  if ((Gsq < Gsq_min_lo) || (r < u)) {
                     Gsq_min_lo = Gsq - epsilon;
                     Gsq_min_hi = Gsq + epsilon;
                     r = u;
                  }
               }
            }
         }
      }

      // Continue from r if the window missed the minimum image
      descendToMinimum<3>(r, d, cell, epsilon);
      return r;
   }

   /*
   * Voronoi-relevant vectors of the lattice of image shifts.
   */
   template <int D>
   void imageShiftNeighbors(IntVec<D> const & d,
                            UnitCell<D> const & cell,
                            FSArray<IntVec<D>, 14>& neighbors)
   {
      neighbors.clear();
      if (D == 1) {
         neighbors.append(d);
         IntVec<D> w;
         w.negate(d);
         neighbors.append(w);
         return;
      }

      // Superbase: v[i+1] = d[i]*e_i, v[0] = -(v[1] + ... + v[D])
      FArray<IntVec<D>, D+1> v;
      for (int i = 0; i < D; ++i) {
         v[i+1] = IntVec<D>(0);
         v[i+1][i] = d[i];
      }
      v[0] = IntVec<D>(0);
      for (int i = 1; i <= D; ++i) {
         v[0] -= v[i];
      }

      // Selling reduction: while some v_i.v_j > 0, apply the reduction
      // step, which lowers the sum of |v_k|^2 by 4 v_i.v_j (2D) or by
      // 2 v_i.v_j (3D) and keeps a superbase of the same lattice.
      int iter = 0;
      for (;;) {
         double norm = 0.0;
         for (int i = 0; i <= D; ++i) {
            norm += kDot<D>(v[i], v[i], cell);
         }
         double pMax = 1.0E-12*norm;
         int iMax = -1;
         int jMax = -1;
         for (int i = 0; i <= D; ++i) {
            for (int j = i + 1; j <= D; ++j) {
               double p = kDot<D>(v[i], v[j], cell);
               if (p > pMax) {
                  pMax = p;
                  iMax = i;
                  jMax = j;
               }
            }
         }
         if (iMax < 0) break;
         for (int k = 0; k <= D; ++k) {
            if (k == iMax || k == jMax) continue;
            v[k] += v[iMax];
            if (D == 2) v[k] += v[iMax];
         }
         v[iMax].negate();
         ++iter;
         if (iter > 10000) {
            UTIL_THROW("Selling reduction of image shifts did not converge");
         }
      }

      // Neighbors: +-v_i, and in 3D +-(v_0 + v_i)
      IntVec<D> w;
      for (int i = 0; i <= D; ++i) {
         neighbors.append(v[i]);
         w.negate(v[i]);
         neighbors.append(w);
      }
      if (D == 3) {
         for (int i = 1; i <= D; ++i) {
            w.add(v[0], v[i]);
            neighbors.append(w);
            w.negate();
            neighbors.append(w);
         }
      }
   }

   template void imageShiftNeighbors<1>(IntVec<1> const &,
                                        UnitCell<1> const &,
                                        FSArray<IntVec<1>, 14>&);
   template void imageShiftNeighbors<2>(IntVec<2> const &,
                                        UnitCell<2> const &,
                                        FSArray<IntVec<2>, 14>&);
   template void imageShiftNeighbors<3>(IntVec<3> const &,
                                        UnitCell<3> const &,
                                        FSArray<IntVec<3>, 14>&);

}
