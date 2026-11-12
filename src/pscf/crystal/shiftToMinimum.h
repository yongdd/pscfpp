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
