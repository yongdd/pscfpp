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
#include <algorithm>
#include <cmath>

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
         double minNormSq; // smallest |k|^2 of the neighbors
         bool isValid;

         NeighborCache()
          : isValid(false)
         {}
      };

      template <int D>
      NeighborCache<D> const &
      neighborCache(IntVec<D> const & d, UnitCell<D> const & cell)
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
            cache.minNormSq = cell.ksq(cache.neighbors[0]);
            for (int j = 1; j < cache.neighbors.size(); ++j) {
               double g = cell.ksq(cache.neighbors[j]);
               if (g < cache.minNormSq) cache.minNormSq = g;
            }
            cache.d = d;
            for (int i = 0; i < D; ++i) {
               for (int j = 0; j < D; ++j) {
                  cache.kBasis[i][j] = cell.kBasis(i)[j];
               }
            }
            cache.isValid = true;
         }
         return cache;
      }

      template <int D>
      FSArray<IntVec<D>, 14> const &
      cachedNeighbors(IntVec<D> const & d, UnitCell<D> const & cell)
      {  return neighborCache<D>(d, cell).neighbors; }

      /*
      * Move image r to a minimum image.
      *
      * Moves r by the image shift neighbor that lowers |k|^2 most, for
      * as long as one lowers it beyond the relative tolerance epsilon.
      * The image is then moved to the largest (in the order of IntVec) of the neighboring
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
      void descendWith(IntVec<D>& r, FSArray<IntVec<D>, 14> const & nb,
                       UnitCell<D> const & cell, double epsilon)
      {
         const int n = nb.size();
         double Gsq = cell.ksq(r);
         IntVec<D> u;
         for (;;) {
            int jBest = -1;
            double GsqBest = Gsq;
            for (int j = 0; j < n; ++j) {
               u.add(r, nb[j]);
               double g = cell.ksq(u);
               double tol = epsilon*std::max(g, GsqBest);
               if (g < GsqBest - tol) {
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
               double g = cell.ksq(u);
               double tol = epsilon*std::max(g, Gsq);
               if ((std::abs(g - Gsq) <= tol) && (r < u)) {
                  r = u;
                  changed = true;
               }
            }
         }
      }

      /*
      * Move image r to a minimum image, by the image shift neighbors
      * of the mesh (see descendWith).
      */
      template <int D>
      void descendToMinimum(IntVec<D>& r, IntVec<D> const & d,
                            UnitCell<D> const & cell, double epsilon)
      {
         FSArray<IntVec<D>, 14> const & nb = cachedNeighbors<D>(d, cell);
         descendWith<D>(r, nb, cell, epsilon);
      }

      /*
      * Is the wave index vi on a Nyquist plane of a mesh dimension di?
      */
      inline bool isNyquistIndex(int vi, int di)
      {
         if (di % 2 != 0) return false;
         int r = vi % di;
         if (r < 0) r += di;
         return (2*r == di);
      }

      /*
      * Voronoi-relevant vectors of the 2D lattice of image shifts in
      * directions j and k: +-v_a for the obtuse superbase v_0, v_1, v_2
      * obtained by Selling reduction of d[j]*e_j, d[k]*e_k and
      * -(d[j]*e_j + d[k]*e_k), as in imageShiftNeighbors.
      */
      template <int D>
      void planeNeighbors(IntVec<D> const & d, int j, int k,
                          UnitCell<D> const & cell,
                          FSArray<IntVec<D>, 14>& nb)
      {
         // Superbase of the plane lattice of image shifts
         FArray<IntVec<D>, 3> v;
         v[1] = IntVec<D>(0);
         v[1][j] = d[j];
         v[2] = IntVec<D>(0);
         v[2][k] = d[k];
         v[0].negate(v[1]);
         v[0] -= v[2];

         // Selling reduction (2D step, as in imageShiftNeighbors)
         int iter = 0;
         for (;;) {
            double norm = 0.0;
            for (int a = 0; a < 3; ++a) {
               norm += kDot<D>(v[a], v[a], cell);
            }
            double pMax = 1.0E-12*norm;
            int aMax = -1;
            int bMax = -1;
            for (int a = 0; a < 3; ++a) {
               for (int b = a + 1; b < 3; ++b) {
                  double p = kDot<D>(v[a], v[b], cell);
                  if (p > pMax) {
                     pMax = p;
                     aMax = a;
                     bMax = b;
                  }
               }
            }
            if (aMax < 0) break;
            int c = 3 - aMax - bMax;
            v[c] += v[aMax];
            v[c] += v[aMax];
            v[aMax].negate();
            ++iter;
            if (iter > 10000) {
               UTIL_THROW("Selling reduction of image shifts did not converge");
            }
         }

         // Voronoi-relevant vectors +-v_a
         nb.clear();
         IntVec<D> w;
         for (int a = 0; a < 3; ++a) {
            nb.append(v[a]);
            w.negate(v[a]);
            nb.append(w);
         }
      }

      /*
      * Plane neighbors (see planeNeighbors) of the last mesh and cell
      * used by this thread, for the plane normal to each direction i
      * (shiftToMinimumNyquist is called once per wavevector, and the
      * waves of the three Nyquist planes alternate).
      */
      template <int D>
      struct PlaneNeighborCache
      {
         IntVec<D> d;
         double kBasis[D][D];
         FSArray<IntVec<D>, 14> neighbors[D];
         bool hasPlane[D];
         bool isValid;

         PlaneNeighborCache()
          : isValid(false)
         {}
      };

      /*
      * Plane neighbors for the directions j and k normal to i (D = 3),
      * computed once per mesh, cell and plane.
      */
      template <int D>
      FSArray<IntVec<D>, 14> const &
      cachedPlaneNeighbors(IntVec<D> const & d, int i, int j, int k,
                           UnitCell<D> const & cell)
      {
         thread_local PlaneNeighborCache<D> cache;
         bool same = cache.isValid && (cache.d == d);
         for (int m = 0; same && m < D; ++m) {
            for (int n = 0; n < D; ++n) {
               if (cache.kBasis[m][n] != cell.kBasis(m)[n]) same = false;
            }
         }
         if (!same) {
            cache.d = d;
            for (int m = 0; m < D; ++m) {
               for (int n = 0; n < D; ++n) {
                  cache.kBasis[m][n] = cell.kBasis(m)[n];
               }
               cache.hasPlane[m] = false;
            }
            cache.isValid = true;
         }
         if (!cache.hasPlane[i]) {
            planeNeighbors<D>(d, j, k, cell, cache.neighbors[i]);
            cache.hasPlane[i] = true;
         }
         return cache.neighbors[i];
      }

      /*
      * Minimum image of r among the images r + s[j]*d[j]*e_j +
      * s[k]*d[k]*e_k (s[j], s[k] integers), with r[i] = 0 for the
      * third direction i.
      *
      * The image is moved by the Voronoi-relevant vectors of the 2D
      * lattice of image shifts in directions j and k (see
      * planeNeighbors), as in descendToMinimum, so the result is the
      * largest minimum image.
      */
      template <int D>
      void planeMinimum(IntVec<D>& r, IntVec<D> const & d,
                        int i, int j, int k,
                        UnitCell<D> const & cell, double epsilon)
      {
         FSArray<IntVec<D>, 14> const & nb
                             = cachedPlaneNeighbors<D>(d, i, j, k, cell);
         descendWith<D>(r, nb, cell, epsilon);
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
      // Compare squared norms relatively: changing length units must not
      // change the selected integer image. Zero is never tied with nonzero.
      // Keep the band at rounding scale; a 1e-8 relative band treats small
      // physical angle perturbations as ties and can break inversion pairs.
      const double epsilon = 1.0E-12;
      double GsqMin = cell.ksq(v);

      // Loop over neighboring images
      IntVec<2> r = v;               // Minimum image of v
      IntVec<2> u;                   // Vector used in loop
      int s0, s1;
      const int q = 2;
      for (s0 = -q; s0 < q+1; ++s0) {
         u[0] = v[0] + s0*d[0];
         for (s1 = -q; s1 < q+1; ++s1) {
            u[1] = v[1] + s1*d[1];
            double Gsq = cell.ksq(u);
            double tol = epsilon*std::max(Gsq, GsqMin);
            if (Gsq < GsqMin - tol ||
                (std::abs(Gsq - GsqMin) <= tol && r < u)) {
               GsqMin = Gsq;
               r = u;
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
      // Compare squared norms relatively: changing length units must not
      // change the selected integer image. Zero is never tied with nonzero.
      // Keep the band at rounding scale; a 1e-8 relative band treats small
      // physical angle perturbations as ties and can break inversion pairs.
      const double epsilon = 1.0E-12;
      double GsqMin = cell.ksq(v);

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
               double Gsq = cell.ksq(u);
               double tol = epsilon*std::max(Gsq, GsqMin);
               if (Gsq < GsqMin - tol ||
                   (std::abs(Gsq - GsqMin) <= tol && r < u)) {
                  GsqMin = Gsq;
                  r = u;
               }
            }
         }
      }

      // Continue from r if the window missed the minimum image
      descendToMinimum<3>(r, d, cell, epsilon);
      return r;
   }

   /*
   * Does any component of v lie on a Nyquist plane?
   */
   template <int D>
   bool hasNyquistIndex(IntVec<D> const & v, IntVec<D> const & d)
   {
      for (int i = 0; i < D; ++i) {
         if (isNyquistIndex(v[i], d[i])) return true;
      }
      return false;
   }

   namespace {

      /*
      * Images tied with r (r first), reached from r by steps of the
      * vectors nb through images with |k|^2 equal to that of r within
      * the relative tolerance TiedImageTolerance (see tiedMinimumImages).
      */
      template <int D>
      void collectTies(IntVec<D> const & r,
                       FSArray<IntVec<D>, 14> const & nb,
                       UnitCell<D> const & cell,
                       FSArray<IntVec<D>, 16>& tied)
      {
         tied.clear();
         tied.append(r);
         double g0 = cell.ksq(r);
         IntVec<D> u;
         for (int t = 0; t < tied.size(); ++t) {
            for (int j = 0; j < nb.size(); ++j) {
               u.add(tied[t], nb[j]);
               double g = cell.ksq(u);
               if (std::abs(g - g0) > TiedImageTolerance*std::max(g, g0)) {
                  continue;
               }
               bool found = false;
               for (int m = 0; m < tied.size(); ++m) {
                  if (tied[m] == u) found = true;
               }
               if (!found) {
                  UTIL_CHECK(tied.size() < tied.capacity());
                  tied.append(u);
               }
            }
         }
      }

      /*
      * Derivatives of |k|^2 with respect to all unit cell parameters,
      * averaged over the tied images (see tiedImageDKSq).
      */
      template <int D>
      void averageDKSq(FSArray<IntVec<D>, 16> const & tied,
                       UnitCell<D> const & cell, double* dksq)
      {
         for (int n = 0; n < cell.nParameter(); ++n) {
            dksq[n] = tiedImageDKSq<D>(tied, cell, n);
         }
      }

      /*
      * Smooth step of the weight of the Nyquist sign average, as a
      * function of an absolute cosine t (see nyquistDecouplingWeight):
      * 1 for t <= NyquistAverageCosineMin, 0 for t >= 
      * NyquistAverageCosineMax, 1 - x^3 (10 - 15x + 6x^2) in between. 
      * Sets ds to the derivative with respect to t.
      */
      inline double nyquistStep(double t, double& ds)
      {
         const double t1 = NyquistAverageCosineMin;
         const double t2 = NyquistAverageCosineMax;
         ds = 0.0;
         if (t <= t1) return 1.0;
         if (t >= t2) return 0.0;
         double x = (t - t1)/(t2 - t1);
         ds = -30.0*x*x*(1.0 - x)*(1.0 - x)/(t2 - t1);
         return 1.0 - x*x*x*(10.0 - 15.0*x + 6.0*x*x);
      }

      /*
      * Smooth steps s[i][j] of the absolute cosines of the angles 
      * between the reciprocal basis vectors i and j != i (see 
      * nyquistDecouplingWeight). If ds is not null, also sets ds[i][j][n]
      * to the derivative of s[i][j] with respect to cell parameter n.
      */
      template <int D>
      void pairSteps(UnitCell<D> const & cell, double (*s)[D],
                     double (*ds)[D][6])
      {
         // Reciprocal metric g[i][j] = b_i.b_j
         double g[D][D];
         for (int i = 0; i < D; ++i) {
            for (int j = 0; j < D; ++j) {
               g[i][j] = 0.0;
               for (int k = 0; k < D; ++k) {
                  g[i][j] += cell.kBasis(i)[k]*cell.kBasis(j)[k];
               }
            }
         }
         int nParam = cell.nParameter();
         for (int i = 0; i < D; ++i) {
            s[i][i] = 1.0;
            for (int j = i + 1; j < D; ++j) {
               double c = g[i][j]/std::sqrt(g[i][i]*g[j][j]);
               double dsdt;
               s[i][j] = nyquistStep(std::abs(c), dsdt);
               s[j][i] = s[i][j];
               if (!ds) continue;
               for (int n = 0; n < nParam; ++n) {
                  double value = 0.0;
                  if (dsdt != 0.0) {
                     double dc = cell.dkkBasis(n, i, j)
                                 / std::sqrt(g[i][i]*g[j][j])
                               - 0.5*c*( cell.dkkBasis(n, i, i)/g[i][i]
                                       + cell.dkkBasis(n, j, j)/g[j][j] );
                     value = dsdt*( c < 0.0 ? -dc : dc );
                  }
                  ds[i][j][n] = value;
                  ds[j][i][n] = value;
               }
            }
         }
      }

      /*
      * Components of v reduced to (-d[i]/2, d[i]/2]. Sets isNyquist[i]
      * and returns the number of Nyquist indices.
      */
      template <int D>
      int reduceComponents(IntVec<D> const & v, IntVec<D> const & d,
                           IntVec<D>& u, bool* isNyquist)
      {
         int nNyquist = 0;
         for (int i = 0; i < D; ++i) {
            int r = v[i] % d[i];
            if (r < 0) r += d[i];
            if (2*r > d[i]) r -= d[i];
            u[i] = r;
            isNyquist[i] = (2*r == d[i]);
            if (isNyquist[i]) ++nNyquist;
         }
         return nNyquist;
      }

      /*
      * Image r of the wavevector with reduced components u (see 
      * reduceComponents) that minimizes |k|^2 of its components in the
      * directions of the block (bit i of block set for direction i),
      * over the shifts of these directions; the other components of r
      * are zero. Returns this minimum, and sets dksq[n] to its 
      * derivatives if dksq is not null, averaged over tied images (see
      * tiedImageDKSq).
      */
      template <int D>
      double blockMinimum(IntVec<D> const & u, IntVec<D> const & d,
                          int block, UnitCell<D> const & cell,
                          IntVec<D>& r, double* dksq)
      {
         int nB = 0;
         FSArray<IntVec<D>, 16> tied;
         r = u;
         for (int i = 0; i < D; ++i) {
            if (block & (1 << i)) {
               ++nB;
            } else {
               r[i] = 0;
            }
         }
         if (nB == D) {
            // Minimum image: descent by the Voronoi-relevant vectors of
            // the lattice of image shifts from the reduced components,
            // which reaches the (largest) minimum image (see descendWith)
            descendToMinimum<D>(r, d, cell, 1.0E-12);
            if (dksq) {
               collectTies<D>(r, cachedNeighbors<D>(d, cell), cell, tied);
            }
         } else if (nB == 2) {
            // D = 3: minimum in the plane of the two directions j, k
            int i = -1;
            int j = -1;
            int k = -1;
            for (int m = 0; m < D; ++m) {
               if (!(block & (1 << m))) {
                  i = m;
               } else if (j < 0) {
                  j = m;
               } else {
                  k = m;
               }
            }
            UTIL_CHECK(i >= 0 && j >= 0 && k >= 0);
            planeMinimum<D>(r, d, i, j, k, cell, 1.0E-12);
            if (dksq) {
               collectTies<D>(r, cachedPlaneNeighbors<D>(d, i, j, k, cell),
                              cell, tied);
            }
         } else if (dksq) {
            // One direction: the reduced component is the minimum
            tied.clear();
            tied.append(r);
         }
         if (dksq) {
            averageDKSq<D>(tied, cell, dksq);
         }
         return cell.ksq(r);
      }

      /*
      * ksqNyquist, its derivatives (if dksq is not null) and the 
      * representative image (if image is not null, see 
      * shiftToMinimumNyquist), for a wavevector with a Nyquist index.
      */
      template <int D>
      double ksqNyquistImpl(IntVec<D> const & v, IntVec<D> const & d,
                            UnitCell<D> const & cell, double* dksq,
                            IntVec<D>* image)
      {
         IntVec<D> u;
         bool isNyquist[D];
         int nNyquist = reduceComponents<D>(v, d, u, isNyquist);
         UTIL_CHECK(nNyquist > 0);
         int nParam = cell.nParameter();

         // Nyquist directions and directions of nonzero indices (masks)
         int nyquistMask = 0;
         int nonzeroMask = 0;
         int nonzero[D];
         int nNonzero = 0;
         for (int i = 0; i < D; ++i) {
            if (isNyquist[i]) nyquistMask |= (1 << i);
            if (u[i] != 0) {
               nonzeroMask |= (1 << i);
               nonzero[nNonzero] = i;
               ++nNonzero;
            }
         }

         // Pairs of directions of nonzero indices, decoupling weights
         double s[D][D];
         double ds[D][D][6];
         pairSteps<D>(cell, s, dksq ? ds : nullptr);
         int pairI[3], pairJ[3];
         int nPair = 0;
         for (int a = 0; a < nNonzero; ++a) {
            for (int b = a + 1; b < nNonzero; ++b) {
               pairI[nPair] = nonzero[a];
               pairJ[nPair] = nonzero[b];
               ++nPair;
            }
         }

         // Minima of blocks of directions (bit masks), each computed 
         // once (blocks recur in different configurations)
         double blockK[1 << D];
         double blockDK[1 << D][6];
         IntVec<D> blockR[1 << D];
         bool hasBlock[1 << D];
         for (int b = 0; b < (1 << D); ++b) hasBlock[b] = false;

         // Sum over configurations of coupled pairs (bits of c)
         double value = 0.0;
         if (dksq) {
            for (int n = 0; n < nParam; ++n) dksq[n] = 0.0;
         }
         double pMax = -1.0;
         double dK[6];
         for (int c = 0; c < (1 << nPair); ++c) {

            // Probability of the configuration, factors of the product
            double f[3];
            double prob = 1.0;
            for (int e = 0; e < nPair; ++e) {
               double se = s[pairI[e]][pairJ[e]];
               f[e] = ((c >> e) & 1) ? 1.0 - se : se;
               prob *= f[e];
            }
            if (prob == 0.0) continue;

            // Connected components of the nonzero directions (masks)
            int comp[D];
            for (int i = 0; i < D; ++i) comp[i] = (1 << i);
            for (int e = 0; e < nPair; ++e) {
               if (!((c >> e) & 1)) continue;
               int merged = comp[pairI[e]] | comp[pairJ[e]];
               for (int i = 0; i < D; ++i) {
                  if (merged & (1 << i)) comp[i] = merged;
               }
            }

            // Blocks: components of Nyquist indices only (averaged over
            // their sign), and all other directions (rest). A component
            // of all nonzero indices is not averaged: reversing all of
            // them is the isometry k -> -k.
            int rest = (1 << D) - 1;
            int blocks[D];
            int nBlock = 0;
            for (int a = 0; a < nNonzero; ++a) {
               int m = comp[nonzero[a]];
               if ((m & nyquistMask) != m) continue;
               if (m == nonzeroMask) continue;
               if (!(rest & m)) continue; // already listed
               blocks[nBlock] = m;
               ++nBlock;
               rest &= ~m;
            }

            // K of this configuration: minimum over the rest, plus the
            // minimum of each averaged component
            int list[D + 1];
            int nList = 0;
            if (rest) {
               list[nList] = rest;
               ++nList;
            }
            for (int b = 0; b < nBlock; ++b) {
               list[nList] = blocks[b];
               ++nList;
            }
            double K = 0.0;
            if (dksq) {
               for (int n = 0; n < nParam; ++n) dK[n] = 0.0;
            }
            IntVec<D> img(0);
            for (int b = 0; b < nList; ++b) {
               int m = list[b];
               if (!hasBlock[m]) {
                  blockK[m] = blockMinimum<D>(u, d, m, cell, blockR[m],
                                              dksq ? blockDK[m] : nullptr);
                  hasBlock[m] = true;
               }
               K += blockK[m];
               img += blockR[m];
               if (dksq) {
                  for (int n = 0; n < nParam; ++n) dK[n] += blockDK[m][n];
               }
            }
            value += prob*K;
            if (image && prob > pMax) {
               pMax = prob;
               *image = img;
            }

            // Derivatives: prob*dK + K*d(prob)
            if (dksq) {
               for (int n = 0; n < nParam; ++n) {
                  double dProb = 0.0;
                  for (int e = 0; e < nPair; ++e) {
                     double term = ds[pairI[e]][pairJ[e]][n];
                     if (term == 0.0) continue;
                     if ((c >> e) & 1) term = -term;
                     for (int g = 0; g < nPair; ++g) {
                        if (g != e) term *= f[g];
                     }
                     dProb += term;
                  }
                  dksq[n] += prob*dK[n] + K*dProb;
               }
            }
         }
         return value;
      }

      /*
      * May ksqNyquist differ from the minimum-image |k|^2 for this 
      * wavevector, i.e. may some configuration of coupled pairs with a
      * positive probability have an averaged component? False if all
      * pairs of nonzero indices have a decoupling weight of 0 (then all
      * nonzero indices form one component, which is not averaged), in
      * particular if there is only one nonzero index.
      */
      template <int D>
      bool usesAverage(IntVec<D> const & u, UnitCell<D> const & cell)
      {
         double s[D][D];
         pairSteps<D>(cell, s, nullptr);
         for (int i = 0; i < D; ++i) {
            if (u[i] == 0) continue;
            for (int j = i + 1; j < D; ++j) {
               if (u[j] != 0 && s[i][j] > 0.0) return true;
            }
         }
         return false;
      }

   }

   /*
   * Decoupling weight of directions i and j.
   */
   template <int D>
   double nyquistDecouplingWeight(int i, int j, UnitCell<D> const & cell)
   {
      UTIL_CHECK(i >= 0 && i < D);
      UTIL_CHECK(j >= 0 && j < D);
      double s[D][D];
      pairSteps<D>(cell, s, nullptr);
      return s[i][j];
   }

   /*
   * Does ksqNyquist average |k|^2 for this wavevector?
   */
   template <int D>
   bool usesNyquistAverage(IntVec<D> const & v, IntVec<D> const & d,
                           UnitCell<D> const & cell)
   {
      IntVec<D> u;
      bool isNyquist[D];
      int nNyquist = reduceComponents<D>(v, d, u, isNyquist);
      if (nNyquist == 0) return false;
      return usesAverage<D>(u, cell);
   }

   /*
   * Square magnitude of a wavevector with Nyquist indices.
   */
   template <int D>
   double ksqNyquist(IntVec<D> const & v, IntVec<D> const & d,
                     UnitCell<D> const & cell)
   {
      if (!hasNyquistIndex(v, d)) {
         return cell.ksq(shiftToMinimum(v, d, cell));
      }
      return ksqNyquistImpl<D>(v, d, cell, nullptr, nullptr);
   }

   /*
   * Square magnitude of a wavevector with Nyquist indices, derivatives.
   */
   template <int D>
   double ksqNyquist(IntVec<D> const & v, IntVec<D> const & d,
                     UnitCell<D> const & cell, double* dksq)
   {
      UTIL_CHECK(dksq);
      if (!hasNyquistIndex(v, d)) {
         IntVec<D> r = shiftToMinimum(v, d, cell);
         FSArray<IntVec<D>, 16> tied;
         collectTies<D>(r, cachedNeighbors<D>(d, cell), cell, tied);
         averageDKSq<D>(tied, cell, dksq);
         return cell.ksq(r);
      }
      return ksqNyquistImpl<D>(v, d, cell, dksq, nullptr);
   }

   /*
   * ksqNyquist, its derivatives and the representative image.
   */
   template <int D>
   double ksqNyquist(IntVec<D> const & v, IntVec<D> const & d,
                     UnitCell<D> const & cell, double* dksq,
                     IntVec<D>& image)
   {
      UTIL_CHECK(hasNyquistIndex(v, d));
      return ksqNyquistImpl<D>(v, d, cell, dksq, &image);
   }

   /*
   * Derivative of ksqNyquist with respect to cell parameter n.
   */
   template <int D>
   double dksqNyquist(IntVec<D> const & v, IntVec<D> const & d,
                      UnitCell<D> const & cell, int n)
   {
      UTIL_CHECK(n >= 0 && n < cell.nParameter());
      double dksq[6];
      ksqNyquist<D>(v, d, cell, dksq);
      return dksq[n];
   }

   /*
   * Representative image of v for ksqNyquist.
   */
   template <int D>
   IntVec<D> shiftToMinimumNyquist(IntVec<D> const & v,
                                   IntVec<D> const & d,
                                   UnitCell<D> const & cell)
   {
      if (!hasNyquistIndex(v, d)) return shiftToMinimum(v, d, cell);
      IntVec<D> image;
      ksqNyquistImpl<D>(v, d, cell, nullptr, &image);
      return image;
   }

   /*
   * Minimum images tied with the minimum image r.
   */
   template <int D>
   void tiedMinimumImages(IntVec<D> const & r, IntVec<D> const & d,
                          UnitCell<D> const & cell,
                          FSArray<IntVec<D>, 16>& tied)
   {
      NeighborCache<D> const & cache = neighborCache<D>(d, cell);

      // An image tied with r + w (w a nonzero lattice vector) has |k|^2
      // of at least |w|^2/4: no search for images inside the sphere 
      // inscribed in the Voronoi cell (most images)
      double limit = (1.0 - 4.0*TiedImageTolerance)*cache.minNormSq;
      if (4.0*cell.ksq(r) < limit) {
         tied.clear();
         tied.append(r);
         return;
      }
      collectTies<D>(r, cache.neighbors, cell, tied);
   }

   /*
   * Derivative of |k|^2 averaged over tied images.
   */
   template <int D>
   double tiedImageDKSq(FSArray<IntVec<D>, 16> const & tied,
                        UnitCell<D> const & cell, int n)
   {
      int nTied = tied.size();
      UTIL_CHECK(nTied > 0);
      double value0 = cell.dksq(tied[0], n);
      if (nTied == 1) return value0;
      double sum = value0;
      double maxDiff = 0.0;
      double scale = std::abs(value0);
      for (int t = 1; t < nTied; ++t) {
         double value = cell.dksq(tied[t], n);
         sum += value;
         maxDiff = std::max(maxDiff, std::abs(value - value0));
         scale = std::max(scale, std::abs(value));
      }
      if (maxDiff <= TiedImageTolerance*scale) return value0;
      return sum/double(nTied);
   }

   template void tiedMinimumImages<1>(IntVec<1> const &, IntVec<1> const &,
                                      UnitCell<1> const &,
                                      FSArray<IntVec<1>, 16>&);
   template void tiedMinimumImages<2>(IntVec<2> const &, IntVec<2> const &,
                                      UnitCell<2> const &,
                                      FSArray<IntVec<2>, 16>&);
   template void tiedMinimumImages<3>(IntVec<3> const &, IntVec<3> const &,
                                      UnitCell<3> const &,
                                      FSArray<IntVec<3>, 16>&);

   template double tiedImageDKSq<1>(FSArray<IntVec<1>, 16> const &,
                                    UnitCell<1> const &, int);
   template double tiedImageDKSq<2>(FSArray<IntVec<2>, 16> const &,
                                    UnitCell<2> const &, int);
   template double tiedImageDKSq<3>(FSArray<IntVec<3>, 16> const &,
                                    UnitCell<3> const &, int);

   template double nyquistDecouplingWeight<1>(int, int, 
                                              UnitCell<1> const &);
   template double nyquistDecouplingWeight<2>(int, int, 
                                              UnitCell<2> const &);
   template double nyquistDecouplingWeight<3>(int, int, 
                                              UnitCell<3> const &);

   template bool usesNyquistAverage<1>(IntVec<1> const &, 
                                       IntVec<1> const &,
                                       UnitCell<1> const &);
   template bool usesNyquistAverage<2>(IntVec<2> const &, 
                                       IntVec<2> const &,
                                       UnitCell<2> const &);
   template bool usesNyquistAverage<3>(IntVec<3> const &, 
                                       IntVec<3> const &,
                                       UnitCell<3> const &);

   template double ksqNyquist<1>(IntVec<1> const &, IntVec<1> const &,
                                 UnitCell<1> const &, double*, 
                                 IntVec<1>&);
   template double ksqNyquist<2>(IntVec<2> const &, IntVec<2> const &,
                                 UnitCell<2> const &, double*, 
                                 IntVec<2>&);
   template double ksqNyquist<3>(IntVec<3> const &, IntVec<3> const &,
                                 UnitCell<3> const &, double*, 
                                 IntVec<3>&);

   template double ksqNyquist<1>(IntVec<1> const &, IntVec<1> const &,
                                 UnitCell<1> const &, double*);
   template double ksqNyquist<2>(IntVec<2> const &, IntVec<2> const &,
                                 UnitCell<2> const &, double*);
   template double ksqNyquist<3>(IntVec<3> const &, IntVec<3> const &,
                                 UnitCell<3> const &, double*);

   template bool hasNyquistIndex<1>(IntVec<1> const &, IntVec<1> const &);
   template bool hasNyquistIndex<2>(IntVec<2> const &, IntVec<2> const &);
   template bool hasNyquistIndex<3>(IntVec<3> const &, IntVec<3> const &);

   template double ksqNyquist<1>(IntVec<1> const &, IntVec<1> const &,
                                 UnitCell<1> const &);
   template double ksqNyquist<2>(IntVec<2> const &, IntVec<2> const &,
                                 UnitCell<2> const &);
   template double ksqNyquist<3>(IntVec<3> const &, IntVec<3> const &,
                                 UnitCell<3> const &);

   template double dksqNyquist<1>(IntVec<1> const &, IntVec<1> const &,
                                  UnitCell<1> const &, int);
   template double dksqNyquist<2>(IntVec<2> const &, IntVec<2> const &,
                                  UnitCell<2> const &, int);
   template double dksqNyquist<3>(IntVec<3> const &, IntVec<3> const &,
                                  UnitCell<3> const &, int);

   template IntVec<1> shiftToMinimumNyquist<1>(IntVec<1> const &,
                                              IntVec<1> const &,
                                              UnitCell<1> const &);
   template IntVec<2> shiftToMinimumNyquist<2>(IntVec<2> const &,
                                              IntVec<2> const &,
                                              UnitCell<2> const &);
   template IntVec<3> shiftToMinimumNyquist<3>(IntVec<3> const &,
                                              IntVec<3> const &,
                                              UnitCell<3> const &);

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
