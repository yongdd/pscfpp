#ifndef PSCF_FFT_CP_TPP
#define PSCF_FFT_CP_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.h"
#include "RField.h"
#include "CField.h"

/*
* A note about const_casts:
*
* The FFTW library is used in this file to perform discrete Fourier
* transforms. FFTW's complex-to-real inverse transform overwrites its
* input array, but all other out-of-place transforms leave the input
* array unaltered. However, all transforms in the FFTW library require
* non-const pointers to the input array, even though they do not alter
* the array.
*
* In order to maintain const-correctness in PSCF, this FFT class accepts
* const input arrays for its methods that perform a Fourier transform,
* unless the transform is expected to modify / overwrite its input (as
* is the case for complex-to-real inverse transforms). This denotes to
* the caller of the method that the input array will not be altered,
* which is an accurate representation of the expected behavior.
*
* However, the const-correctness of this FFT class creates a conflict
* with the FFTW library's requirement of non-const inputs. This conflict
* is resolved using a const_cast, in which the const pointer to the input
* array is made non-const when passed into FFTW functions. The use of
* const_cast is reserved only for these few select cases in which we are
* confident that the input array will not be modified.
*
* For more information about the relevant FFTW methods, see the FFTW
* documentation at https://www.fftw.org/fftw3_doc/index.html. In Section
* 4.3.2, it is stated that, by default, "an out-of-place transform must
* not change its input array," except for complex-to-real transforms, in
* which case "no input-preserving algorithms are implemented." Finally,
* we note that the unit tests for this FFT class check that the input
* array is unaltered, allowing developers to continually ensure that
* the FFTW functions do not modify their input unexpectedly.
*/

namespace Pscf {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   FFT<D,CPT>::FFT()
    : meshDimensions_(0),
      rSize_(0),
      kSize_(0),
      rcfPlan_(0),
      criPlan_(0),
      ccfPlan_(0),
      cciPlan_(0),
      rcfUnalignedPlan_(0),
      criUnalignedPlan_(0),
      ccfUnalignedPlan_(0),
      cciUnalignedPlan_(0),
      rAlignment_(0),
      kAlignment_(0),
      cInAlignment_(0),
      cOutAlignment_(0),
      isSetup_(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   FFT<D,CPT>::~FFT()
   {
      if (rcfPlan_) {
         fftw_destroy_plan(rcfPlan_);
      }
      if (criPlan_) {
         fftw_destroy_plan(criPlan_);
      }
      if (ccfPlan_) {
         fftw_destroy_plan(ccfPlan_);
      }
      if (cciPlan_) {
         fftw_destroy_plan(cciPlan_);
      }
      if (rcfUnalignedPlan_) fftw_destroy_plan(rcfUnalignedPlan_);
      if (criUnalignedPlan_) fftw_destroy_plan(criUnalignedPlan_);
      if (ccfUnalignedPlan_) fftw_destroy_plan(ccfUnalignedPlan_);
      if (cciUnalignedPlan_) fftw_destroy_plan(cciUnalignedPlan_);
   }

   /*
   * Setup mesh dimensions, work memory and FFT plans.
   */
   template <int D>
   void FFT<D,CPT>::setup(IntVec<D> const & meshDimensions)
   {
      // Precondition
      UTIL_CHECK(!isSetup_);

      // Set and check mesh dimensions
      rSize_ = 1;
      kSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         rSize_ *= meshDimensions[i];
         if (i < D - 1) {
            kSize_ *= meshDimensions[i];
         } else {
            kSize_ *= (meshDimensions[i]/2 + 1);
         }
      }

      // Allocate kFieldCopy_ array if necessary
      if (!kFieldCopy_.isAllocated()) {
          kFieldCopy_.allocate(meshDimensions);
      } else {
          if (kFieldCopy_.capacity() != kSize_) {
             kFieldCopy_.deallocate();
             kFieldCopy_.allocate(meshDimensions);
          }
      }
      UTIL_CHECK(meshDimensions == kFieldCopy_.meshDimensions());
      UTIL_CHECK(kFieldCopy_.capacity() == kSize_);

      // Create temporary RField and CField objects used for plans
      RField<D,CPT> rField;
      rField.allocate(meshDimensions);
      UTIL_CHECK(meshDimensions == rField.meshDimensions());
      UTIL_CHECK(rField.capacity() == rSize_);

      CField<D,CPT> cFieldIn;
      cFieldIn.allocate(meshDimensions);
      UTIL_CHECK(meshDimensions == cFieldIn.meshDimensions());
      UTIL_CHECK(cFieldIn.capacity() == rSize_);

      CField<D,CPT> cFieldOut;
      cFieldOut.allocate(meshDimensions);
      UTIL_CHECK(meshDimensions == cFieldOut.meshDimensions());
      UTIL_CHECK(cFieldOut.capacity() == rSize_);

      // Make aligned plans and fallback plans for associated slices
      makePlans(rField, kFieldCopy_, cFieldIn, cFieldOut);

      isSetup_ = true;
   }

   /*
   * Make normal plans and alignment-independent fallback plans.
   *
   * FFTW new-array execution requires the planning arrays' alignment.
   * An associated field slice need not share it. Keep optimized plans for
   * the usual allocated arrays, selecting FFTW_UNALIGNED only when needed.
   * Both kinds are created during setup, not inside transform execution.
   */
   template <int D>
   void FFT<D,CPT>::makePlans(RField<D,CPT>& rField,
                              RFieldDft<D,CPT>& kField,
                              CField<D,CPT>& cFieldIn,
                              CField<D,CPT>& cFieldOut)
   {
      int dims[D];
      for (int i = 0; i < D; ++i) dims[i] = meshDimensions_[i];
      rAlignment_ = fftw_alignment_of(&rField[0]);
      kAlignment_ = fftw_alignment_of(&kField[0][0]);
      cInAlignment_ = fftw_alignment_of(&cFieldIn[0][0]);
      cOutAlignment_ = fftw_alignment_of(&cFieldOut[0][0]);

      unsigned int flags = FFTW_ESTIMATE;
      rcfPlan_ = fftw_plan_dft_r2c(D, dims, &rField[0], &kField[0], flags);
      criPlan_ = fftw_plan_dft_c2r(D, dims, &kField[0], &rField[0], flags);
      ccfPlan_ = fftw_plan_dft(D, dims, &cFieldIn[0], &cFieldOut[0],
                             FFTW_FORWARD, flags);
      cciPlan_ = fftw_plan_dft(D, dims, &cFieldOut[0], &cFieldIn[0],
                             FFTW_BACKWARD, flags);
      flags |= FFTW_UNALIGNED;
      rcfUnalignedPlan_ = fftw_plan_dft_r2c(D, dims, &rField[0],
                                          &kField[0], flags);
      criUnalignedPlan_ = fftw_plan_dft_c2r(D, dims, &kField[0],
                                          &rField[0], flags);
      ccfUnalignedPlan_ = fftw_plan_dft(D, dims, &cFieldIn[0],
                                       &cFieldOut[0], FFTW_FORWARD, flags);
      cciUnalignedPlan_ = fftw_plan_dft(D, dims, &cFieldOut[0],
                                       &cFieldIn[0], FFTW_BACKWARD, flags);
      UTIL_CHECK(rcfPlan_ && criPlan_ && ccfPlan_ && cciPlan_);
      UTIL_CHECK(rcfUnalignedPlan_ && criUnalignedPlan_ &&
                 ccfUnalignedPlan_ && cciUnalignedPlan_);
   }

   // Real <-> Complex Transforms

   /*
   * Execute real-to-complex forward transform.
   */
   template <int D>
   void FFT<D,CPT>::forwardTransform(RField<D,CPT> const & rField,
                                 RFieldDft<D,CPT>& kField)
   const
   {
      UTIL_CHECK(isSetup_)
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == kSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);

      // Execute preplanned forward transform
      // (See note at top of file explaining this use of const_cast)
      double* in = const_cast<double*>(&rField[0]);
      bool aligned = fftw_alignment_of(in) == rAlignment_ &&
                     fftw_alignment_of(&kField[0][0]) == kAlignment_;
      fftw_execute_dft_r2c(aligned ? rcfPlan_ : rcfUnalignedPlan_,
                           in, &kField[0]);

      // Rescale the resulting array
      double scale = 1.0/double(rSize_);
      for (int i = 0; i < kSize_; ++i) {
         kField[i][0] *= scale;
         kField[i][1] *= scale;
      }
   }

   /*
   * Execute inverse (complex-to-real) transform.
   */
   template <int D>
   void
   FFT<D,CPT>::inverseTransformUnsafe(
                              RFieldDft<D,CPT> & kField, 
                              RField<D,CPT>& rField)
   const
   {
      UTIL_CHECK(isSetup_)
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == kSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);

      // Execute preplanned inverse transform
      bool aligned = fftw_alignment_of(&kField[0][0]) == kAlignment_ &&
                     fftw_alignment_of(&rField[0]) == rAlignment_;
      fftw_execute_dft_c2r(aligned ? criPlan_ : criUnalignedPlan_,
                           &kField[0], &rField[0]);
   }

   /*
   * Execute inverse (complex-to-real) transform without destroying input.
   */
   template <int D>
   void
   FFT<D,CPT>::inverseTransformSafe(
                              RFieldDft<D,CPT> const & kField,
                              RField<D,CPT>& rField) const
   {
      UTIL_CHECK(kFieldCopy_.capacity() == kField.capacity());

      kFieldCopy_ = kField;
      inverseTransformUnsafe(kFieldCopy_, rField);
   }

   // Complex <-> Complex Transforms

   /*
   * Execute complex-to-complex forward transform.
   */
   template <int D>
   void
   FFT<D,CPT>::forwardTransform(
                           CField<D,CPT> const & rField, 
                           CField<D,CPT>& kField) const
   {
      UTIL_CHECK(isSetup_)
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.capacity() == rSize_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);

      // Execute preplanned forward transform
      // (See note at top of file explaining this use of const_cast)
      fftw_complex* in = const_cast<fftw_complex*>(&rField[0]);
      bool aligned = fftw_alignment_of(&in[0][0]) == cInAlignment_ &&
                     fftw_alignment_of(&kField[0][0]) == cOutAlignment_;
      fftw_execute_dft(aligned ? ccfPlan_ : ccfUnalignedPlan_,
                       in, &kField[0]);

      // Rescale the resulting array
      double scale = 1.0/double(rSize_);
      for (int i = 0; i < rSize_; ++i) {
         kField[i][0] *= scale;
         kField[i][1] *= scale;
      }
   }

   /*
   * Execute inverse (complex-to-complex) transform.
   */
   template <int D>
   void
   FFT<D,CPT>::inverseTransform(
                           CField<D,CPT> const & kField, 
                           CField<D,CPT>& rField)
   const
   {
      UTIL_CHECK(isSetup_)
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == rSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);

      // Execute preplanned inverse transform
      // (See note at top of file explaining this use of const_cast)
      fftw_complex* in = const_cast<fftw_complex*>(&kField[0]);
      bool aligned = fftw_alignment_of(&in[0][0]) == cOutAlignment_ &&
                     fftw_alignment_of(&rField[0][0]) == cInAlignment_;
      fftw_execute_dft(aligned ? cciPlan_ : cciUnalignedPlan_,
                       in, &rField[0]);
   }

   // Static function

   /*
   * Compute dimensions and size of the k-size mesh (static)
   */
   template <int D>
   void FFT<D,CPT>::computeKMesh(
                             IntVec<D> const & rMeshDimensions,
                             IntVec<D> & kMeshDimensions,
                             int & kSize )
   {
      kSize = 1;
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions[i] = rMeshDimensions[i];
         } else {
            kMeshDimensions[i] = rMeshDimensions[i]/2 + 1;
         }
         kSize *= kMeshDimensions[i];
      }
   }

} // namespace Pscf
#endif
