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
      // field of its monomer type minus the spatial average of that
      // field, which is removed from ln q analytically by
      // PolymerTmpl::solve. This keeps the propagators representable
      // for long chains and makes the results independent of the
      // constant of w.
      int monomerId;
      double wShift;
      for (int j = 0; j < nBlock(); ++j) {
         monomerId = block(j).monomerId();
         DArray<double> const & w = wFields[monomerId];
         wShift = block(j).domain().spatialAverage(w);
         block(j).setupSolver(w, wShift);
      }

      // Solve MDE for all propagators
      solve();

      // Compute block concentration fields (thread model), normalizing
      // products of the shifted propagators by their own partition
      // function
      double prefactor;
      prefactor = phi() / ( qShifted() * length() );
      for (int i = 0; i < nBlock(); ++i) {
         block(i).computeConcentration(prefactor);
      }

   }

}
}
