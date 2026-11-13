/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.h"
#include <r1/field/Domain.h>

#include <cmath>

namespace Pscf {
namespace R1 { 

   /*
   * Constructor.
   */
   Solvent::Solvent() 
    : SolventSpecies()  
   {  setClassName("Solvent"); }

   /*
   * Destructor.
   */
   Solvent::~Solvent()
   {}

   /*
   * Set association with Domain and allocate concentration field.
   */
   void Solvent::setDiscretization(Domain const & domain)
   {
      domainPtr_ = &domain;
      int nx = domain.nx();
      if (nx > 0) {
         cField_.allocate(nx);
      }
   }

   /*
   * Compute concentration, q, and phi or mu.
   */ 
   void Solvent::compute(FieldT const & wField)
   {
      UTIL_CHECK(cField_.isAllocated());

      // Shift by the minimum field so that every Boltzmann weight is
      // at most one and at least one is exactly one. Subtracting the
      // average instead can overflow for finite nonuniform fields.
      int nx = domain().nx();
      UTIL_CHECK(wField.capacity() == nx);
      double s = size();
      double wShift = wField[0];
      for (int i = 1; i < nx; ++i) {
         if (wField[i] < wShift) wShift = wField[i];
      }
      for (int i = 0; i < nx; ++i) {
          cField_[i] = exp(-s*(wField[i] - wShift));
      }

      // Compute spatial average Q of the shifted weight
      double Q = domain().spatialAverage(cField_);

      // Set ln q = ln Q - s*wShift, compute mu or phi (by ensemble)
      Species::setLnQ(std::log(Q) - s*wShift);

      // Normalize concentration (the shift cancels in phi/Q)
      double prefactor = phi()/Q;
      for (int i = 0; i < nx; ++i) {
          cField_[i] *= prefactor;
      }
    
   }

}
}
