#ifndef R1_POLYMER_H
#define R1_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/PolymerTmpl.h>   // base class template
#include <util/containers/DArray.h>      // member template

// Forward declarations
namespace Pscf {
   namespace R1 {
      class Block;
      class Propagator;
   }
}

namespace Pscf {
namespace R1 {

   using namespace Util;

   /**
   * Descriptor and solver for a block polymer species.
   *
   * A Polymer contains an array of Block objects, each of which contains
   * two propagator objects. After Polymer::compute is called, each block
   * stores the block concentration (volume fraction) fields computed for 
   * that block, and each Propagator contains the solution of the modified
   * diffusion equation for that propagator. Functions that return values
   * for molecular volume fraction phi, molecular chemical potential mu, 
   * and molecular partition function q are indirectly inherited from 
   * Pscf::PolymerSpecies.
   * 
   * \ref user_param_polymer_sec "Parameter File Format"
   * \ingroup R1_Solver_Module
   */
   class Polymer : public PolymerTmpl<Block, Propagator, double>
   {

   public:

      /**
      * Default constructor.
      */
      Polymer();

      /**
      * Destructor.
      */
      ~Polymer();

      /**
      * Compute solution to modified diffusion equation and concentrations.
      *
      * Upon return, propagator solutions and block concentration fields
      * have been computed for all propagators and blocks, along with the
      * molecular partition function q and phi or mu.
      *
      * Each block solves the MDE in the field w of its monomer type
      * minus a constant shift (see Block::setupSolver), so the stored
      * propagators are those of the shifted fields. Their partition
      * function is qShifted(); ln q, mu, phi and the block
      * concentrations do not depend on the shifts. The shift is
      *
      *    min(w) - ln < exp(-N (w - min(w))) > / N,
      *
      * where < > is the spatial (volume) average and N is length():
      * the free energy per monomer of a whole molecule held at one
      * point. It lies between the minimum and the average of w, so
      * that the propagators cannot overflow and qShifted() is at most
      * about one. It follows a constant added to w exactly, and it is
      * close to the field where the molecules are, which keeps the
      * Crank-Nicolson contour error small. If qShifted() is
      * nevertheless below about 1e-294 (a long chain in a field that
      * varies strongly on scales it cannot follow), the shifts are
      * raised so that qShifted() is one. First each block is calibrated
      * alone: its propagator is solved from a uniform head, and its
      * shift is raised until the spatial average of the tail is one
      * (in steps of about 678/l, for a block of length l, while that
      * average is below 1e-294, and then by minus its logarithm over
      * l). Then all shifts are raised by a common constant so that
      * qShifted() is one. Raising only the blocks that need it keeps
      * the propagators of the other blocks bounded. If that result is
      * not representable (a block whose head in the polymer is
      * concentrated where the chain sits is raised by up to the
      * logarithm of the reciprocal occupied volume fraction too much,
      * which adds up over many junctions), the documented shifts are
      * instead raised by a common constant only. A raised shift can
      * exceed the field by more than 2/ds where the field is low, which
      * makes the Crank-Nicolson factor negative there; this only
      * happens for fields that vary by hundreds of kT per unit length
      * on scales a chain cannot follow. A NanException is thrown if the
      * partition function or the concentrations are still not
      * representable, except for fields that are not finite (a shift
      * that is not finite), which give NaN values as before.
      *
      * \param wFields  array of chemica potential fields.
      */
      void compute(DArray< DArray<double> > const & wFields);

   private:

      /// Work array for the field shifts (allocated by compute).
      DArray<double> wWork_;

      // Restrict access
      using PolymerTmpl<Block, Propagator, double>::solve;

   };

} // namespace R1
} // namespace Pscf

// Explicit instantiation declaration
namespace Pscf {
   extern template class PolymerTmpl<R1::Block, R1::Propagator, double>;
}
#endif
