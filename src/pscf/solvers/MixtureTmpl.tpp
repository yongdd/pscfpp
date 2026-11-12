#ifndef PSCF_MIXTURE_TMPL_TPP
#define PSCF_MIXTURE_TMPL_TPP

/*
* PSCF - Mixture Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "MixtureTmpl.h"
#include <util/misc/Log.h>
#include <cmath>
#include <complex>

namespace Pscf
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <class PT, class ST, typename WT>
   MixtureTmpl<PT,ST,WT>::MixtureTmpl()
    : Composition<WT>(),
      ParamComposite(),
      polymers_(),
      solvents_()
   {}

   /*
   * Get a PolymerSpecies descriptor by non-const reference.
   */
   template <class PT, class ST, typename WT>
   PolymerSpecies<WT> const & 
   MixtureTmpl<PT,ST,WT>::polymerSpecies(int id) const
   {
      UTIL_CHECK(id >= 0 && id < nPolymer_);
      return polymers_[id];
   }

   /*
   * Get a SolventSpecies descriptor by const reference.
   */
   template <class PT, class ST, typename WT>
   SolventSpecies<WT> const & 
   MixtureTmpl<PT,ST,WT>::solventSpecies(int id) const
   {
      UTIL_CHECK(id >= 0 && id < nSolvent_);
      return solvents_[id];
   }

   /*
   * Read all parameters and initialize.
   */
   template <class PT, class ST, typename WT>
   void MixtureTmpl<PT,ST,WT>::readParameters(std::istream& in)
   {
      // Read nMonomer and monomers array
      read<int>(in, "nMonomer", nMonomer_);
      monomers_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         monomers_[i].setId(i);
      }
      readDArray< Monomer >(in, "monomers", monomers_, nMonomer_);

      /*
      * The input format for a single monomer is defined in the istream
      * extraction operation (operator >>) for a Pscf::Monomer, in file
      * pscf/chem/Monomer.cpp. The text representation contains only the
      * value for the monomer statistical segment Monomer::kuhn.
      */

      // Read nPolymer (required parameter, must be > 0)
      read<int>(in, "nPolymer", nPolymer_);
      UTIL_CHECK(nPolymer_ > 0);

      // Optionally read nSolvent, with nSolvent=0 by default
      nSolvent_ = 0;
      readOptional<int>(in, "nSolvent", nSolvent_);

      // Read polymers and compute nBlock
      nBlock_ = 0;
      if (nPolymer_ > 0) {

         polymers_.allocate(nPolymer_);
         for (int i = 0; i < nPolymer_; ++i) {
            readParamComposite(in, polymer(i));
            nBlock_ = nBlock_ + polymer(i).nBlock();
         }

         // Set statistical segment lengths for all blocks
         double kuhn;
         int monomerId;
         for (int i = 0; i < nPolymer_; ++i) {
            for (int j = 0; j < polymer(i).nBlock(); ++j) {
               monomerId = polymer(i).block(j).monomerId();
               kuhn = Composition<WT>::monomer(monomerId).kuhn();
               polymer(i).block(j).setKuhn(kuhn);
            }
         }

      }

      // Read solvents
      if (nSolvent_ > 0) {

         solvents_.allocate(nSolvent_);
         for (int i = 0; i < nSolvent_; ++i) {
            readParamComposite(in, solvent(i));
            // The monomerId is used as an array index by the solvers
            int monomerId = solvent(i).monomerId();
            if (monomerId < 0 || monomerId >= nMonomer_) {
               Log::file() << "Solvent " << i << ": monomerId = "
                           << monomerId << ", nMonomer = " << nMonomer_
                           << "\n";
               UTIL_THROW("Solvent monomerId out of range");
            }
         }

      }

      // Optionally read monomer reference value
      vMonomer_ = 1.0; // Default value
      readOptional(in, "vMonomer", vMonomer_);
      if (vMonomer_ <= 0.0) {
         Log::file() << "vMonomer = " << vMonomer_ << "\n";
         UTIL_THROW("vMonomer must be positive");
      }

      // If every species is closed, the volume fractions must add up
      // to 1 (see normalizePhi)
      if (Composition<WT>::isCanonical()) {
         WT sum = normalizePhi();
         if (std::abs(sum - 1.0) > 1.0E-12) {
            Log::file() << "Sum of species volume fractions = " << sum
                        << ", rescaled to 1\n";
         }
      }

   }

   /*
   * Rescale the volume fractions of a closed mixture to add up to 1.
   */
   template <class PT, class ST, typename WT>
   WT MixtureTmpl<PT,ST,WT>::normalizePhi()
   {
      // If every species is closed, the volume fractions must add up
      // to 1. The SCFT iterators remove the homogeneous part of the
      // residual in this case, so they never test incompressibility at
      // k = 0 and would converge to a state with a total concentration
      // equal to the sum of the volume fractions.
      UTIL_CHECK(Composition<WT>::isCanonical());
      WT sum = 0.0;
      for (int i = 0; i < nPolymer_; ++i) {
         sum += polymer(i).phi();
      }
      for (int i = 0; i < nSolvent_; ++i) {
         sum += solvent(i).phi();
      }
      if (std::abs(sum - 1.0) > 1.0E-5) {
         Log::file() << "Sum of species volume fractions = " << sum
                     << "\n";
         UTIL_THROW("Volume fractions of a closed mixture must add up to 1");
      }

      // Within the tolerance (e.g., 1/3 written as 0.333333), rescale
      // to a sum of exactly 1: the FTS compressors keep the k = 0
      // residual, which would otherwise stay at |sum - 1|
      if (std::abs(sum - 1.0) > 1.0E-12) {
         const double s = std::real(sum);
         for (int i = 0; i < nPolymer_; ++i) {
            polymer(i).setPhi(std::real(polymer(i).phi())/s);
         }
         for (int i = 0; i < nSolvent_; ++i) {
            solvent(i).setPhi(std::real(solvent(i).phi())/s);
         }
      }
      return sum;
   }

}
#endif
