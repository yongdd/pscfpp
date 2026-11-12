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

      // Evaluate unnormalized concentration, the Boltzmann weight of
      // the field shifted by its spatial average wShift, which keeps
      // the weight representable for large solvent molecules
      int nx = domain().nx();
      double s = size();
      double wShift = domain().spatialAverage(wField);
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
