/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.tpp"

namespace Pscf {

   using namespace Util;

   // Explicit class instantiation definitions
   template class FFT<1,CPT>;
   template class FFT<2,CPT>;
   template class FFT<3,CPT>;

} // namespace Pscf
