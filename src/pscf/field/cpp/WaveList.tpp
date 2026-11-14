#ifndef PSCF_WAVE_LIST_CP_TPP
#define PSCF_WAVE_LIST_CP_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "WaveList.h"

#include <pscf/field/cpp/FFT.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/crystal/hasVariableAngle.h>
#include <pscf/crystal/waveNorm.h>

#include <pscf/mesh/MeshIterator.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/Sort.h>

namespace Pscf {

   /*
   * Constructor.
   */
   template <int D>
   WaveList<D,CPT>::WaveList(bool isRealField)
    : kSize_(0),
      nBunch_(0),
      isAllocated_(false),
      hasMinImages_(false),
      hasKSq_(false),
      hasdKSq_(false),
      isSorted_(false),
      isRealField_(isRealField),
      averageNyquist_(false),
      unitCellPtr_(nullptr),
      meshPtr_(nullptr)
   {}

   /*
   * Destructor.
   */
   template <int D>
   WaveList<D,CPT>::~WaveList() 
   {}

   /*
   * Allocate memory used by WaveList.
   */
   template <int D>
   void WaveList<D,CPT>::allocate(Mesh<D> const & m, UnitCell<D> const & c) 
   {
      UTIL_CHECK(m.size() > 0);
      UTIL_CHECK(c.nParameter() > 0);
      UTIL_CHECK(!isAllocated_);

      // Create permanent associations with mesh and unit cell
      unitCellPtr_ = &c;
      meshPtr_ = &m;

      // Local copies of properties
      int nParams = unitCell().nParameter();
      IntVec<D> const & meshDimensions = mesh().dimensions();

      // Average kSq over the signs of Nyquist indices if the lattice
      // system has a variable angle (see class documentation)
      averageNyquist_ = hasVariableAngle<D>(unitCell().lattice());

      // Compute kMeshDimensions_ and kSize_
      if (isRealField_) {
         FFT<D,CPT>::computeKMesh(meshDimensions, kMeshDimensions_, kSize_);
      } else {
         kMeshDimensions_ = meshDimensions;
         kSize_ = mesh().size();
      }

      // Allocate memory
      minImages_.allocate(kSize_);
      kSq_.allocate(kMeshDimensions_);
      dKSq_.allocate(nParams);
      for (int i = 0; i < nParams; i++) {
         dKSq_[i].allocate(kMeshDimensions_);
      }

      // Allocate and set up implicitInverse_ array if isRealField_ == true
      // (only depends on mesh dimensions, only used for real fields)
      if (isRealField_) {
         implicitInverse_.allocate(kSize_);

         MeshIterator<D> kItr(kMeshDimensions_);
         int rank;

         for (kItr.begin(); !kItr.atEnd(); ++kItr) {
            rank = kItr.rank();
            implicitInverse_[rank] = 
               FFT<D,CPT>::hasImplicitInverse(kItr.position(), meshDimensions);
         }
      }

      clearUnitCellData();
      isAllocated_ = true;
   }

   /*
   * Clear all data that depends on unit cell parameters.
   */
   template <int D>
   void WaveList<D,CPT>::clearUnitCellData()
   {
      hasKSq_ = false;
      hasdKSq_ = false;
      if (unitCell().nParameter() > 1) {
         isSorted_ = false;
         sortedBunches_.clear();
         nBunch_ = 0;
      }
      if (hasVariableAngle<D>(unitCell().lattice())) {
         hasMinImages_ = false;
      }
   }

   /*
   * Allow or forbid averaging kSq over the signs of Nyquist indices.
   */
   template <int D>
   void WaveList<D,CPT>::allowNyquistAverage(bool allow)
   {
      UTIL_CHECK(isAllocated_);
      bool average = allow && hasVariableAngle<D>(unitCell().lattice());
      if (average != averageNyquist_) {
         averageNyquist_ = average;
         hasMinImages_ = false;
         hasKSq_ = false;
         hasdKSq_ = false;
         isSorted_ = false;
         sortedBunches_.clear();
         nBunch_ = 0;
      }
   }

   /*
   * Compute minimum image vectors for all wavevectors.
   */
   template <int D>
   void WaveList<D,CPT>::computeMinimumImages() 
   {
      if (hasMinImages_) return; // min images already calculated

      // Precondition
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().isInitialized());
      UTIL_CHECK(minImages_.capacity() == kSize_);

      IntVec<D> const & d = mesh().dimensions();
      int const nParam = unitCell().nParameter();

      // Wavevectors with a Nyquist index for which ksqNyquist averages
      // |k|^2 (see usesNyquistAverage): kSq, image and dKSq
      // from ksqNyquist, in one evaluation; dKSq is kept for 
      // computedKSq. Otherwise ksqNyquist is the minimum-image |k|^2, 
      // and these wavevectors are treated as all others.
      nyquistRanks_.clear();
      nyquistDKSq_.clear();
      double dk[6];
      IntVec<D> image;

      // Wavevectors with tied minimum images (see tiedRanks_)
      bool const variableAngle = hasVariableAngle<D>(unitCell().lattice());
      tiedRanks_.clear();
      tiedImages_.clear();
      FSArray<IntVec<D>, 16> tied;

      MeshIterator<D> kItr(kMeshDimensions_);
      int rank;
      for (kItr.begin(); !kItr.atEnd(); ++kItr) {
         rank = kItr.rank();
         IntVec<D> const & v = kItr.position();
         if (averageNyquist_ && hasNyquistIndex(v, d)) {
            if (usesNyquistAverage(v, d, unitCell())) {
               kSq_[rank] = ksqNyquist(v, d, unitCell(), dk, image);
               minImages_[rank] = image;
               nyquistRanks_.push_back(rank);
               for (int n = 0; n < nParam; ++n) {
                  nyquistDKSq_.push_back(dk[n]);
               }
               continue;
            }
         }
         minImages_[rank] = shiftToMinimum(v, d, unitCell());
         kSq_[rank] = unitCell().ksq(minImages_[rank]);

         // Record tied minimum images (variable angles only: for other
         // lattice systems, tied images have equal dKSq)
         if (variableAngle) {
            tiedMinimumImages(minImages_[rank], d, unitCell(), tied);
            if (tied.size() > 1) {
               tiedRanks_.push_back(rank);
               tiedImages_.push_back(tied);
            }
         }
      }

      hasMinImages_ = true;
      hasKSq_ = true;
   }

   /*
   * Compute array of value of |k|^2
   */
   template <int D>
   void WaveList<D,CPT>::computeKSq() 
   {
      // If kSq_ is valid, return immediately without recomputing
      if (hasKSq_) return; 

      // If necessary, compute minimum images.
      if (!hasMinImages_) {
         computeMinimumImages(); // computes both min images and kSq
         return;
      }

      // Preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(unitCell().nParameter() > 0);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().isInitialized());

      // Compute kSq_ (ksqNyquist for the wavevectors in nyquistRanks_)
      IntVec<D> const & d = mesh().dimensions();
      for (int rank = 0; rank < kSize_; ++rank) {
         kSq_[rank] = unitCell().ksq(minImages_[rank]);
      }
      int nNyquist = nyquistRanks_.size();
      for (int k = 0; k < nNyquist; ++k) {
         int rank = nyquistRanks_[k];
         kSq_[rank] = ksqNyquist(minImages_[rank], d, unitCell());
      }
      
      hasKSq_ = true;
   }

   /*
   * Compute derivatives of |k|^2 w/ respect to unit cell parameters.
   */
   template <int D>
   void WaveList<D,CPT>::computedKSq()
   {
      if (hasdKSq_) return; // dKSq already calculated

      // Compute minimum images if needed
      if (!hasMinImages_) {
         computeMinimumImages(); 
      }

      // Preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(unitCell().nParameter() > 0);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().isInitialized());

      int const nParam = unitCell().nParameter();
      MeshIterator<D> kItr(kMeshDimensions_);
      int i, rank;
      for (i = 0 ; i < nParam; ++i) {
         RField<D,CPT>& dksq = dKSq_[i];
         for (kItr.begin(); !kItr.atEnd(); ++kItr) {
            rank = kItr.rank();
            dksq[rank] = unitCell().dksq(minImages_[rank], i);
            if (isRealField_) {
               if (implicitInverse_[rank]) {
                  dksq[rank] *= 2.0;
               }
            }
         }
      }

      // Wavevectors with tied minimum images: average over these images
      int nTiedWaves = tiedRanks_.size();
      for (int t = 0; t < nTiedWaves; ++t) {
         rank = tiedRanks_[t];
         FSArray<IntVec<D>, 16> const & tied = tiedImages_[t];
         for (i = 0; i < nParam; ++i) {
            double value = tiedImageDKSq(tied, unitCell(), i);
            if (isRealField_) {
               if (implicitInverse_[rank]) {
                  value *= 2.0;
               }
            }
            dKSq_[i][rank] = value;
         }
      }

      // Wavevectors in nyquistRanks_: derivatives of ksqNyquist, from
      // computeMinimumImages (for the same cell)
      int nNyquist = nyquistRanks_.size();
      UTIL_CHECK((int) nyquistDKSq_.size() == nNyquist*nParam);
      for (int k = 0; k < nNyquist; ++k) {
         rank = nyquistRanks_[k];
         for (i = 0; i < nParam; ++i) {
            dKSq_[i][rank] = nyquistDKSq_[k*nParam + i];
            if (isRealField_) {
               if (implicitInverse_[rank]) {
                  dKSq_[i][rank] *= 2.0;
               }
            }
         }
      }

      hasdKSq_ = true;
   }

   /*
   * Sort waves by magnitude.
   */
   template <int D>
   void WaveList<D,CPT>::sortWaves()
   {
      if (isSorted_) return;

      // Compute wavenumbers if needed
      UTIL_CHECK(isAllocated_);
      if (!hasKSq_) {
         computeKSq();
      }

      // Construct sorted array of items 
      std::vector< Sort::Item<double> > items;
      Sort::Item<double> item;
      for (int i = 0; i < kSize_; ++i) {
         item.value = kSq_[i];
         item.id = i;
         items.push_back(item);
      }
      UTIL_CHECK((int)items.size() == kSize_);
      Sort::sort(items);

      // Fill sortedIds_ array with ids
      if (!sortedIds_.isAllocated()) {
         sortedIds_.allocate(kSize_);
      }
      UTIL_CHECK(sortedIds_.capacity() == kSize_);
      for (int i = 0; i < kSize_; ++i) {
         sortedIds_[i] = items[i].id;
      }

      // Construct sortedBunches_ array and set nBunch_
      sortedBunches_.clear();
      findWaveBunches(items, sortedBunches_, waveNormTolerance<double>());
      nBunch_ = sortedBunches_.size();

      // Fill bunchIds_ array
      UTIL_CHECK(nBunch_ > 0);
      if (!bunchIds_.isAllocated()) {
         bunchIds_.allocate(kSize_);
      }
      UTIL_CHECK(bunchIds_.capacity() == kSize_);
      int begin = 0;
      int end = 0;
      int ib, iw;
      for (ib = 0; ib < nBunch_; ++ib) {
         begin = sortedBunches_[ib][0];
         end = sortedBunches_[ib][1];
         UTIL_CHECK(end > begin);
         for (iw = begin; iw < end; ++iw) {
            bunchIds_[sortedIds_[iw]] = ib;
         }
      }
      UTIL_CHECK(end == kSize_);

      isSorted_ = true;
   }

}
#endif
