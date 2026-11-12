#ifndef RP_INTRACORRELATION_TPP
#define RP_INTRACORRELATION_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraCorrelation.h"

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Solvent.h>
#include <rp/field/Domain.h>
#include <pscf/chem/Ensemble.h>

#include <pscf/field/FFT.h>
#include <pscf/field/RField.h>
#include <pscf/field/WaveList.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/correlation/Mixture.h>

#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   IntraCorrelation<D,T>::IntraCorrelation(
                                     System<D,T> const & system)
    : systemPtr_(&system),
      correlationMixturePtr_(nullptr),
      kSize_(-1)
   {
      correlationMixturePtr_
          = new Pscf::Correlation::Mixture<RealT>(system.mixture());
   }

   /*
   * Destructor.
   */
   template <int D, class T>
   IntraCorrelation<D,T>::~IntraCorrelation()
   {  delete correlationMixturePtr_; }

   /*
   * Compute k-space array of intramolecular correlation functions (on device).
   */
   template<int D, class T>
   void
   IntraCorrelation<D,T>::computeOmegaTotal(RField<D,T>& correlations)
   {
      getMeshDimensions();
      int nk = kSize();
      UTIL_CHECK(correlations.isAllocated());
      UTIL_CHECK(correlations.capacity() == nk);

      correlations_h_.associate(correlations);
      computeOmegaTotalArray(correlations_h_);
      correlations = correlations_h_;
      // Release the association, so that the correlations can be
      // computed again (e.g., after a ramp step changes parameters)
      correlations_h_.dissociate();
   }

   /*
   * Compute k-space array of intramolecular correlation functions (on host).
   */
   template <int D, class T>
   void
   IntraCorrelation<D,T>::computeOmegaTotalArray(Array<RealT>& correlations)
   {
      getMeshDimensions();
      UTIL_CHECK(correlations.capacity() == kSize_);
      computeGsq();
      if (!correlationMixturePtr_->isAllocated()) {
         correlationMixturePtr_->allocate();
      }
      correlationMixturePtr_->setup();
      correlationMixturePtr_->computeOmegaTotal(Gsq_, correlations);

      // The LR compressors use these values as the response of the
      // total concentration to a change of the pressure-like field. In a
      // mixed ensemble (some species open, some closed) the amount of
      // each closed species is fixed, so at k = 0 only the open species
      // respond: use the sum of their k = 0 terms, phi*N/v each. (In a
      // canonical system the k = 0 residual vanishes, and in a
      // grand-canonical one every species counts, so the value is kept.)
      Mixture<D,T> const & mixture = system().mixture();
      const double vMonomer = mixture.vMonomer();
      bool hasOpen = false;
      bool hasClosed = false;
      double omegaOpen = 0.0;
      for (int i = 0; i < mixture.nPolymer(); ++i) {
         if (mixture.polymer(i).ensemble() == Ensemble::Open) {
            Correlation::Polymer<RealT> const & polymer
                                     = correlationMixturePtr_->polymer(i);
            omegaOpen += polymer.phi()*polymer.totalLength()/vMonomer;
            hasOpen = true;
         } else {
            hasClosed = true;
         }
      }
      for (int i = 0; i < mixture.nSolvent(); ++i) {
         if (mixture.solvent(i).ensemble() == Ensemble::Open) {
            omegaOpen += mixture.solvent(i).phi()*mixture.solvent(i).size()
                         / vMonomer;
            hasOpen = true;
         } else {
            hasClosed = true;
         }
      }
      if (hasOpen && hasClosed) {
         for (int j = 0; j < kSize_; ++j) {
            if (Gsq_[j] == 0.0) {
               correlations[j] = omegaOpen;
            }
         }
      }
   }

   /*
   * Compute k-grid mesh dimensions and Gsq.
   */
   template <int D, class T>
   void IntraCorrelation<D,T>::getMeshDimensions()
   {
      meshDimensions_ = system().domain().mesh().dimensions();

      // Compute k-space mesh dimensions kMeshDimensions_ and size Size_
      FFT<D,T>::computeKMesh(meshDimensions_, kMeshDimensions_, kSize_);
   }

   /*
   * Construct array of squared wavevector values (on host).
   */
   template <int D, class T>
   void IntraCorrelation<D,T>::computeGsq()
   {
      // Get mesh dimensions if not done previously
      if (kSize_ <= 0) {
         getMeshDimensions();
      }

      // Check allocation of Gsq_ (k-space array of square wavenumbers)
      if (!Gsq_.isAllocated()) {
         Gsq_.allocate(kSize_);
      }
      UTIL_CHECK(Gsq_.capacity() == kSize_);

      // Compute Gsq_, equal to kSq of the WaveList used by the MDE 
      // solver: for a wavevector with a Nyquist index, averaged over 
      // the signs of the Nyquist indices if the WaveList does so (see
      // WaveList::allowNyquistAverage), otherwise |k|^2 of the minimum
      // image. 
      UnitCell<D> const & unitCell = system().domain().unitCell();
      bool const averageNyquist 
                       = system().domain().waveList().hasNyquistAverage();
      IntVec<D> G, Gmin;
      MeshIterator<D> iter;
      iter.setDimensions(kMeshDimensions_);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         G = iter.position();
         if (averageNyquist && hasNyquistIndex(G, meshDimensions_)) {
            Gmin = shiftToMinimumNyquist(G, meshDimensions_, unitCell);
            Gsq_[iter.rank()] = ksqNyquist(Gmin, meshDimensions_, unitCell);
         } else {
            Gmin = shiftToMinimum(G, meshDimensions_, unitCell);
            Gsq_[iter.rank()] = unitCell.ksq(Gmin);
         }
      }
   }

}
}
#endif
