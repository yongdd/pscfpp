#ifndef PSCF_SHIFT_MINIMUM_H
#define PSCF_SHIFT_MINIMUM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/crystal/UnitCell.h>
#include <pscf/math/IntVec.h>
#include <util/containers/FSArray.h>
#include <iostream>

namespace Pscf {

   using namespace Util;

   /**
   * Returns minimum magnitude image of DFT wavevector.
   *
   * The images of v differ by vectors of the lattice of image shifts,
   * spanned by d[i]*e_i in reciprocal basis coordinates. The function
   * first searches the shifts s[i]*d[i] with |s[i]| <= 2 and keeps,
   * among images of equal |k|^2, the one that is largest in the order
   * of IntVec. For strongly anisotropic meshes on skewed cells, this
   * window can miss the minimum image. The function therefore then
   * moves the image by the Voronoi-relevant vectors of the shift
   * lattice (see imageShiftNeighbors) for as long as one of them lowers
   * |k|^2, and then to the largest of the neighboring images with equal
   * |k|^2. The result is the largest minimum image of the class of v,
   * whatever representative v is. It differs from the result of the
   * window search only where that search missed the minimum or a larger
   * tied image.
   *
   * \ingroup Pscf_Crystal_Module
   *
   * \param v  IntVec<D> containing integer indices of a wavevector
   * \param d  dimensions of the discrete Fourier transform grid
   * \param cell  UnitCell
   * \return  integer indices of minimum image of v
   */
   template <int D>
   IntVec<D> shiftToMinimum(IntVec<D> const & v, 
                            IntVec<D> const & d, 
                            UnitCell<D> const & cell);

   // Explicit specializations for D=1, 2 and 3

   template <> 
   IntVec<1> shiftToMinimum(IntVec<1> const & v, 
                            IntVec<1> const & d, 
                            UnitCell<1> const & cell);

   template <> 
   IntVec<2> shiftToMinimum(IntVec<2> const & v, 
                            IntVec<2> const & d, 
                            UnitCell<2> const & cell);

   template <>
   IntVec<3> shiftToMinimum(IntVec<3> const & v,
                            IntVec<3> const & d,
                            UnitCell<3> const & cell);

   /**
   * Does any index of a DFT wavevector lie on a Nyquist plane?
   *
   * Index v[i] is a Nyquist index if d[i] is even and v[i] is congruent
   * to d[i]/2 modulo d[i]. The images of such a wavevector with v[i] 
   * and -v[i] are then equivalent (they differ by v[i] - (-v[i]), a 
   * multiple of d[i]).
   *
   * \ingroup Pscf_Crystal_Module
   *
   * \param v  integer indices of a wavevector
   * \param d  dimensions of the discrete Fourier transform grid
   */
   template <int D>
   bool hasNyquistIndex(IntVec<D> const & v, IntVec<D> const & d);

   /**
   * Absolute cosine below which two directions are fully decoupled.
   *
   * See nyquistDecouplingWeight.
   *
   * \ingroup Pscf_Crystal_Module
   */
   const double NyquistAverageCosineMin = 0.1;

   /**
   * Absolute cosine above which two directions are fully coupled.
   *
   * See nyquistDecouplingWeight.
   *
   * \ingroup Pscf_Crystal_Module
   */
   const double NyquistAverageCosineMax = 0.25;

   /**
   * Decoupling weight of two directions for the Nyquist sign average.
   *
   * Returns a smooth step function of t = |b_i.b_j|/(|b_i||b_j|), where
   * b_i are the reciprocal basis vectors: 1 for t <= 
   * NyquistAverageCosineMin (an angle within about 5.7 degrees of 90
   * degrees), 0 for t >= NyquistAverageCosineMax (about 14.5 degrees), 
   * and 1 - x^3 (10 - 15x + 6x^2) in between, with x the fractional 
   * position of t between these bounds, so that it has continuous first
   * and second derivatives (see ksqNyquist).
   *
   * \ingroup Pscf_Crystal_Module
   *
   * \param i  index of a direction (0 <= i < D)
   * \param j  index of another direction (0 <= j < D, j != i)
   * \param cell  UnitCell
   */
   template <int D>
   double nyquistDecouplingWeight(int i, int j, UnitCell<D> const & cell);

   /**
   * Square magnitude |k|^2 of a wavevector with Nyquist indices.
   *
   * Reversing the signs of a set F of the Nyquist indices of a 
   * wavevector v (see hasNyquistIndex) changes |k|^2 only through the 
   * cross terms 2 v[i] v[j] b_i.b_j with i in F and j not in F, where
   * b_i are the reciprocal basis vectors. Where all these b_i.b_j with
   * v[j] != 0 vanish, the images that differ by such a reversal have
   * equal |k|^2, and where the angles pass through 90 degrees they 
   * exchange as minimum image, so that the minimum |k|^2 has a kink 
   * there. Near such geometries, kSq is therefore averaged over these
   * reversals, which is smooth there, and far from them the minimum
   * image is used: the average is invariant only under symmetry 
   * operations that map Nyquist indices onto Nyquist indices (signed
   * permutations of the basis vectors), and breaks the symmetry of 
   * cells whose metric has other operations, e.g. of a rhombic cell at
   * 120 degrees (hexagonal lattice) or a rhombohedral cell at 60 or 
   * 109.47 degrees (fcc or bcc lattice), whereas the minimum image has
   * the symmetry of every lattice.
   *
   * Precisely: each pair of directions i, j with nonzero indices of v
   * (modulo d) is coupled with probability 1 - w_ij, where w_ij = 
   * nyquistDecouplingWeight(i, j, cell), and decoupled otherwise. For
   * each configuration of coupled pairs, the connected components of
   * the directions with nonzero indices whose indices are all Nyquist
   * indices are averaged over their sign, except a component of all 
   * nonzero indices (whose reversal is the isometry k -> -k): their 
   * cross terms with the other directions are dropped, |k|^2 of each 
   * such component is minimized over the shifts of its directions, and
   * |k|^2 of all other directions (including those with a zero index)
   * over the shifts of these. This function returns the sum of these 
   * values weighted by the probabilities of the configurations 
   * (products of w_ij and 1 - w_ij). Near orthogonal geometries (all 
   * weights 1), this is the average of |k|^2 over the signs of all 
   * Nyquist indices (if v has another nonzero index); where the weights
   * of all pairs of nonzero indices are 0, and for a wavevector with
   * only one nonzero index, it is the minimum-image |k|^2. The 
   * result depends only on the class of v (v modulo d), is continuous 
   * with continuous first and second derivatives where the weights 
   * change, and is smooth where an angle between directions with 
   * nonzero indices passes through 90 degrees whatever the other angles
   * are. For a wavevector v with no Nyquist index, this returns |k|^2 
   * of its minimum image.
   *
   * The cosines of the symmetric geometries named above (1/2 for 
   * hexagonal and bcc, 1/3 for fcc) lie in the range where the weights
   * are 0, so that kSq has the symmetry of these cells, as it has for
   * operations that are signed permutations. A lattice whose symmetry
   * operations include others that are not signed permutations at 
   * angles closer to 90 degrees (e.g. a centered rectangular lattice in
   * a primitive oblique cell with very unequal sides) is not treated
   * symmetrically by the average.
   *
   * \ingroup Pscf_Crystal_Module
   *
   * \param v  integer indices of a wavevector (any image)
   * \param d  dimensions of the discrete Fourier transform grid
   * \param cell  UnitCell
   */
   template <int D>
   double ksqNyquist(IntVec<D> const & v, IntVec<D> const & d,
                     UnitCell<D> const & cell);

   /**
   * Does ksqNyquist average |k|^2 for this wavevector?
   *
   * Returns true if v has a Nyquist index and some pair of directions
   * with nonzero indices has a positive decoupling weight, so that some
   * configuration of coupled pairs with a positive probability may have
   * an averaged component (see ksqNyquist). If it returns false, 
   * ksqNyquist is the minimum-image |k|^2.
   *
   * \ingroup Pscf_Crystal_Module
   *
   * \param v  integer indices of a wavevector
   * \param d  dimensions of the discrete Fourier transform grid
   * \param cell  UnitCell
   */
   template <int D>
   bool usesNyquistAverage(IntVec<D> const & v, IntVec<D> const & d,
                           UnitCell<D> const & cell);

   /**
   * Square magnitude of a wavevector with Nyquist indices, derivatives.
   *
   * Returns ksqNyquist(v, d, cell), and sets dksq[n] to its derivative
   * with respect to unit cell parameter n, for all parameters.
   *
   * \ingroup Pscf_Crystal_Module
   *
   * \param v  integer indices of a wavevector (any image)
   * \param d  dimensions of the discrete Fourier transform grid
   * \param cell  UnitCell
   * \param dksq  output, array of cell.nParameter() derivatives
   */
   template <int D>
   double ksqNyquist(IntVec<D> const & v, IntVec<D> const & d,
                     UnitCell<D> const & cell, double* dksq);

   /**
   * Square magnitude, derivatives and image of a Nyquist wavevector.
   *
   * Returns ksqNyquist(v, d, cell), sets dksq[n] to its derivatives if
   * dksq is not null, and sets image to shiftToMinimumNyquist(v, d, 
   * cell), all in one evaluation.
   *
   * \pre v has a Nyquist index (see hasNyquistIndex)
   *
   * \ingroup Pscf_Crystal_Module
   *
   * \param v  integer indices of a wavevector (any image)
   * \param d  dimensions of the discrete Fourier transform grid
   * \param cell  UnitCell
   * \param dksq  output, array of cell.nParameter() derivatives, or null
   * \param image  output, the representative image
   */
   template <int D>
   double ksqNyquist(IntVec<D> const & v, IntVec<D> const & d,
                     UnitCell<D> const & cell, double* dksq,
                     IntVec<D>& image);

   /**
   * Derivative of ksqNyquist with respect to a unit cell parameter.
   *
   * \ingroup Pscf_Crystal_Module
   *
   * \param v  integer indices of a wavevector (any image)
   * \param d  dimensions of the discrete Fourier transform grid
   * \param cell  UnitCell
   * \param n  index of the unit cell parameter
   */
   template <int D>
   double dksqNyquist(IntVec<D> const & v, IntVec<D> const & d,
                      UnitCell<D> const & cell, int n);

   /**
   * Returns a representative image of a DFT wavevector for ksqNyquist.
   *
   * If v has no Nyquist index, this is shiftToMinimum(v, d, cell). 
   * Otherwise it is the image of the configuration of coupled pairs of
   * the largest probability (see ksqNyquist): the minimizing images of
   * its averaged components (Nyquist indices at +d[i]/2 for a component
   * of one direction) and of the other directions. Where all weights 
   * are 1, ksqNyquist is the average of |k|^2 of this image over the 
   * signs of its Nyquist indices; where ksqNyquist is the minimum-image
   * |k|^2, it is the (largest) minimum image.
   *
   * \ingroup Pscf_Crystal_Module
   *
   * \param v  integer indices of a wavevector
   * \param d  dimensions of the discrete Fourier transform grid
   * \param cell  UnitCell
   * \return  integer indices of the image
   */
   template <int D>
   IntVec<D> shiftToMinimumNyquist(IntVec<D> const & v,
                                   IntVec<D> const & d,
                                   UnitCell<D> const & cell);

   /**
   * Voronoi-relevant vectors of the lattice of image shifts.
   *
   * The lattice of image shifts is spanned by the vectors d[i]*e_i, in
   * reciprocal basis coordinates, with the metric of the reciprocal
   * lattice of the cell. Selling reduction turns this basis into an
   * obtuse superbase v_0, ..., v_D (sum zero, all v_i.v_j <= 0). The
   * Voronoi-relevant vectors are then among the sums of the elements of
   * proper nonempty subsets of the superbase (Conway and Sloane, Proc.
   * R. Soc. Lond. A 436, 55 (1992)): +-v_i in 2D (6 vectors), and +-v_i
   * and +-(v_0+v_i) in 3D (14 vectors). An image k is a minimum image
   * if and only if no k + w, with w one of these vectors, has smaller
   * |k|^2. In 1D the vectors are +-d[0].
   *
   * \ingroup Pscf_Crystal_Module
   *
   * \param d  dimensions of the discrete Fourier transform grid
   * \param cell  UnitCell
   * \param neighbors  output, the vectors (in a fixed order)
   */
   template <int D>
   void imageShiftNeighbors(IntVec<D> const & d,
                            UnitCell<D> const & cell,
                            FSArray<IntVec<D>, 14>& neighbors);

}
#endif
