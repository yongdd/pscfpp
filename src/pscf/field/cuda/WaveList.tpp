#ifndef PSCF_WAVE_LIST_CU_TPP
#define PSCF_WAVE_LIST_CU_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "WaveList.h"

#include <pscf/field/cuda/resources.h>
#include <pscf/field/cuda/FFT.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/crystal/hasVariableAngle.h>
#include <pscf/crystal/waveNorm.h>
#include <pscf/backend/cuda/HostRecvArray.h>
#include <pscf/backend/cuda/HostSendArray.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/Sort.h>
#include <pscf/math/IntVec.h>

namespace Pscf {

   // CUDA kernels:
   // (defined in an anonymous namespace, used only in this file)
   namespace {

      /*
      * Voronoi-relevant vectors of the lattice of image shifts (see
      * imageShiftNeighbors in shiftToMinimum.h), passed to the kernels
      * by value. Vector j has components v[D*j], ..., v[D*j+D-1].
      */
      struct ImageShifts
      {
         int n;
         int v[14*3];
      };

      /*
      * Square magnitude of the wavevector with integer coordinates img.
      */
      template <int D>
      __device__ cudaReal _kSqImage(int const * img,
                                    cudaReal const * kBasis)
      {
         cudaReal kSq = 0.0;
         for (int c = 0; c < D; ++c) {
            cudaReal kVec = 0.0;
            for (int k = 0; k < D; ++k) {
               kVec += kBasis[D*k + c] * img[k];
            }
            kSq += kVec*kVec;
         }
         return kSq;
      }

      /*
      * Store the images and kSq of wavevectors with a Nyquist index,
      * computed on the host (see ksqNyquist).
      *
      * Wavevector i < n of the list has rank ranks[i], its image has
      * components images[i + j*stride], j = 0, ..., D-1, and its kSq 
      * is kSqValues[i].
      */
      template <int D> __global__
      void _setNyquistImages(int* minImages,
                             cudaReal* kSq,
                             int const * ranks,
                             int const * images,
                             cudaReal const * kSqValues,
                             int const n,
                             int const stride,
                             int const kSize)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            int rank = ranks[i];
            for (int j = 0; j < D; ++j) {
               minImages[rank + j*kSize] = images[i + j*stride];
            }
            kSq[rank] = kSqValues[i];
         }
      }

      /*
      * Store derivatives dKSq computed on the host for wavevectors with
      * a Nyquist index: out[ranks[i] + p*kSize] = values[i + p*stride]
      * for i < n and p < nSlice (the number of unit cell parameters).
      */
      __global__
      void _setNyquistValues(cudaReal* out,
                             int const * ranks,
                             cudaReal const * values,
                             int const n,
                             int const stride,
                             int const nSlice,
                             int const kSize)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            int rank = ranks[i];
            for (int p = 0; p < nSlice; ++p) {
               out[rank + p*kSize] = values[i + p*stride];
            }
         }
      }

      /*
      * Move image r to a minimum image, with the steps of the function
      * descendToMinimum in shiftToMinimum.cpp: move by the image shift
      * that lowers kSq most beyond the relative tolerance epsilon,
      * then move to the largest neighboring image of equal kSq,
      * repeatedly. Returns true if r changed, and then sets kSq.
      */
      template <int D>
      __device__ bool _descendToMinimum(int* r, cudaReal& kSq,
                                        cudaReal const * kBasis,
                                        ImageShifts const & shifts,
                                        cudaReal const epsilon)
      {
         int u[D];
         bool moved = false;
         for (;;) {
            int jBest = -1;
            cudaReal kSqBest = kSq;
            for (int j = 0; j < shifts.n; ++j) {
               for (int i = 0; i < D; ++i) {
                  u[i] = r[i] + shifts.v[D*j + i];
               }
               cudaReal g = _kSqImage<D>(u, kBasis);
               cudaReal tol = epsilon * (g > kSqBest ? g : kSqBest);
               if (g < kSqBest - tol) {
                  jBest = j;
                  kSqBest = g;
               }
            }
            if (jBest < 0) break;
            for (int i = 0; i < D; ++i) {
               r[i] += shifts.v[D*jBest + i];
            }
            kSq = kSqBest;
            moved = true;
         }
         bool changed = true;
         while (changed) {
            changed = false;
            for (int j = 0; j < shifts.n; ++j) {
               for (int i = 0; i < D; ++i) {
                  u[i] = r[i] + shifts.v[D*j + i];
               }
               cudaReal g = _kSqImage<D>(u, kBasis);
               cudaReal tol = epsilon * (g > kSq ? g : kSq);
               if (fabs(g - kSq) <= tol) {
                  // Is r < u in the order of IntVec?
                  bool less = false;
                  for (int i = 0; i < D; ++i) {
                     if (r[i] > u[i]) break;
                     if (r[i] < u[i]) {
                        less = true;
                        break;
                     }
                  }
                  if (less) {
                     for (int i = 0; i < D; ++i) {
                        r[i] = u[i];
                     }
                     changed = true;
                     moved = true;
                  }
               }
            }
         }
         if (moved) {
            kSq = _kSqImage<D>(r, kBasis);
         }
         return moved;
      }

      /*
      * Compute minimum images of each wavevector, template declaration.
      */
      template <int D> __global__
      void _computeMinimumImages(int* minImages,
                                 cudaReal* kSq,
                                 cudaReal const * kBasis,
                                 int const * meshDims,
                                 int const kSize,
                                 int const kBegin,
                                 ImageShifts const shifts);

      /*
      * Compute minimum images of each wavevector, 1D specialization.
      *
      * Kernel must be launched with kSize threads. (One thread calculates
      * one minimum image.) Wavevectors with index below kBegin are not
      * evaluated.
      *
      * When launched, this kernel requires no dynamic memory allocation.
      */
      template <> __global__
      void _computeMinimumImages<1>(int* minImages,
                                    cudaReal* kSq,
                                    cudaReal const * kBasis,
                                    int const * meshDims,
                                    int const kSize,
                                    int const kBegin,
                                    ImageShifts const shifts)
      {
         unsigned int nThreads = blockDim.x * gridDim.x;
         unsigned int startID = blockIdx.x * blockDim.x + threadIdx.x;

         // Load kBasis and meshDims
         cudaReal kBasis_ = kBasis[0];
         int meshDims_ = meshDims[0];

         for (int i = kBegin + startID; i < kSize; i += nThreads) {
            int img = minImages[i];

            while (img > (meshDims_>>1)) { // note: x>>1 is same as x/2
               img -= meshDims_;
            }
            while (img < -1*(meshDims_>>1)) {
               img += meshDims_;
            }

            minImages[i] = img;

            // Compute kSq in floating point (img*img overflows int
            // for |img| > 46340, i.e. for meshes of 92682 or more points)
            cudaReal kVec = kBasis_ * img;
            kSq[i] = kVec * kVec;
         }
      }

      /*
      * Compute minimum images of each wavevector, 2D specialization.
      *
      * This kernel should be launched with >=32*kSize threads.
      *
      * 32 threads are used to calculate one minimum image. In the CPU code,
      * we compare 25 different images of a wave to determine the minimum
      * image, and here we choose 32 instead because it is a power of 2,
      * allowing us to use a reduction algorithm to compare the images.
      * The 7 spare threads evaluate the unshifted image again, so that
      * the images compared are the same as in the CPU code. As in the
      * CPU code, the image is then moved by the shifts in shifts while
      * one lowers kSq (needed where the window misses the minimum).
      *
      * When launched, this kernel requires dynamic memory allocation of
      * (2*sizeof(int) + sizeof(cudaReal)) * nThreadsPerBlock.
      */
      template <> __global__
      void _computeMinimumImages<2>(int* minImages,
                                    cudaReal* kSq,
                                    cudaReal const * kBasis,
                                    int const * meshDims,
                                    int const kSize,
                                    int const kBegin,
                                    ImageShifts const shifts)
      {
         unsigned int startID = blockIdx.x * blockDim.x + threadIdx.x;

         // Determine which lattice parameter and which image to evaluate
         // (this launch treats wavevectors kBegin, kBegin+1, ...)
         unsigned int localID = (startID >> 5); // equivalent to startID / 32
         unsigned int paramID = kBegin + localID;
         unsigned int imageID = startID - (localID * 32); // value in [0, 31]

         // Determine the image that will be evaluated by this thread
         // (uses integer division, not ideal for speed)
         int s0, s1;
         if (imageID < 25) {
            s0 = imageID / 5;
            s1 = imageID - (s0 * 5);
            s0 -= 2; // shift values to go from -2 to 2, not 0 to 4
            s1 -= 2;
         } else {
            s0 = 0;  // spare thread: unshifted image
            s1 = 0;
         }

         // Relative squared-norm tolerance, independent of length units
         #ifdef SINGLE_PRECISION
         cudaReal epsilon = 1.0E-4;
         #else
         cudaReal epsilon = 1.0E-12;
         #endif

         // Declare dynamically allocated arrays in shared memory
         // (allocated space is shared between cudaReal array and int array)
         extern __shared__ cudaReal kSqVals_[];
         int* images_ = (int*)&kSqVals_[blockDim.x];

         if (paramID < kSize) { // only evaluate if on the k-grid

            // Load image from global memory and shift based on s0 and s1
            images_[2*threadIdx.x] = minImages[paramID] + (s0*meshDims[0]);
            images_[2*threadIdx.x+1] = minImages[kSize + paramID]
                                     + (s1 * meshDims[1]);

            // Calculate kSq for this wave
            cudaReal kVec0 = (kBasis[0] * images_[2*threadIdx.x]) +
                             (kBasis[2] * images_[2*threadIdx.x+1]);
            cudaReal kVec1 = (kBasis[1] * images_[2*threadIdx.x]) +
                             (kBasis[3] * images_[2*threadIdx.x+1]);

            kSqVals_[threadIdx.x] = (kVec0*kVec0) + (kVec1*kVec1);
         }
         __syncthreads(); // wait for all threads to finish

         // Perform a parallel reduction on 32 threads to find min image
         // (note: stride >>= 1 is equivalent to stride /= 2)
         for (int stride = 16; stride > 0; stride >>= 1) {
            if (paramID < kSize) { // only evaluate if on the k-grid
               if (imageID < stride) {
                  bool swap = false;
                  cudaReal a = kSqVals_[threadIdx.x];
                  cudaReal b = kSqVals_[threadIdx.x+stride];
                  cudaReal tol = epsilon * (a > b ? a : b);
                  if (b < a - tol) {
                     swap = true;
                  } else if (fabs(b - a) <= tol) {
                     // kSq values effectively equal.
                     // Determine whether to swap based on hkl indices
                     for (int i = 0; i < 2; ++i) {
                        if (images_[2*threadIdx.x+i] >
                                       images_[2*(threadIdx.x+stride)+i]) {
                           break;
                        } else
                        if (images_[2*threadIdx.x+i] <
                                       images_[2*(threadIdx.x+stride)+i]) {
                           swap = true;
                           break;
                        }
                     }
                  }
                  if (swap) {
                     images_[2*threadIdx.x] =
                                       images_[2*(threadIdx.x+stride)];
                     images_[2*threadIdx.x+1] =
                                       images_[2*(threadIdx.x+stride)+1];
                     kSqVals_[threadIdx.x] = kSqVals_[threadIdx.x+stride];
                  }
               }
            }
            // The 32 threads of this reduction form one warp. Since Volta
            // (independent thread scheduling), the warp must still be
            // synchronized between steps, so that each step reads the
            // values written in the previous one.
            __syncwarp();
         }

         // At this point, for any thread with imageID == 0, the
         // corresponding entries in kSqVals_ and images_ should contain
         // the kSq value and the minimum image, respectively. Store
         // values and exit
         if ((imageID == 0) && (paramID < kSize)) {
            int r[2];
            r[0] = images_[2*threadIdx.x];
            r[1] = images_[2*threadIdx.x+1];
            cudaReal kSqMin = kSqVals_[threadIdx.x];
            _descendToMinimum<2>(r, kSqMin, kBasis, shifts, epsilon);
            kSq[paramID] = kSqMin;
            minImages[paramID] = r[0];
            minImages[paramID + kSize] = r[1];
         }
      }

      /*
      * Compute minimum images of each wavevector, 3D specialization.
      *
      * This kernel should be launched with >=128*kSize threads.
      *
      * 128 threads are used to calculate one minimum image. In the CPU
      * code, we compare 125 different images of a wave to determine the
      * minimum image, and here we choose 128 instead because it is a
      * power of 2, allowing us to use a reduction algorithm to compare
      * the images. The 3 spare threads evaluate the unshifted image
      * again, so that the images compared are the same as in the CPU
      * code. As in the CPU code, the image is then moved by the shifts
      * in shifts while one lowers kSq.
      *
      * When launched, this kernel requires dynamic memory allocation of
      * (3*sizeof(int) + sizeof(cudaReal)) * nThreadsPerBlock.
      */
      template <> __global__
      void _computeMinimumImages<3>(int* minImages,
                                    cudaReal* kSq,
                                    cudaReal const * kBasis,
                                    int const * meshDims,
                                    int const kSize,
                                    int const kBegin,
                                    ImageShifts const shifts)
      {
         unsigned int startID = blockIdx.x * blockDim.x + threadIdx.x;

         // Determine which lattice parameter and which image to evaluate
         // (this launch treats wavevectors kBegin, kBegin+1, ...)
         unsigned int localID = (startID >> 7); // equivalent to startID / 128
         unsigned int paramID = kBegin + localID;
         unsigned int imageID = startID - (localID * 128);

         // Determine the image that will be evaluated by this thread
         // (uses integer division, not ideal for speed)
         int s0, s1, s2;
         if (imageID < 125) {
            s0 = imageID / 25;
            s1 = (imageID - (s0 * 25)) / 5;
            s2 = imageID - (s0 * 25) - (s1 * 5);
            s0 -= 2; // shift values to go from -2 to 2
            s1 -= 2;
            s2 -= 2;
         } else {
            s0 = 0;  // spare thread: unshifted image
            s1 = 0;
            s2 = 0;
         }

         // Relative squared-norm tolerance, independent of length units
         #ifdef SINGLE_PRECISION
         cudaReal epsilon = 1.0E-4;
         #else
         cudaReal epsilon = 1.0E-12;
         #endif

         // Declare dynamically allocated arrays in shared memory
         // (allocated space is shared between cudaReal array and int array)
         extern __shared__ cudaReal kSqVals_[];
         int* images_ = (int*)&kSqVals_[blockDim.x];

         if (paramID < kSize) { // only evaluate if on the k-grid

            // Load image data from global memory
            images_[3*threadIdx.x] = minImages[paramID] + (s0 * meshDims[0]);
            images_[3*threadIdx.x+1] = minImages[kSize + paramID] +
                                       (s1 * meshDims[1]);
            images_[3*threadIdx.x+2] = minImages[kSize + kSize + paramID] +
                                       (s2 * meshDims[2]);

            // Calculate kSq for this wave
            cudaReal kVec0(0.0), kVec1(0.0), kVec2(0.0);
            for (int k = 0; k < 3; ++k) {
               kVec0 += kBasis[3*k]   * images_[3*threadIdx.x + k];
               kVec1 += kBasis[3*k+1] * images_[3*threadIdx.x + k];
               kVec2 += kBasis[3*k+2] * images_[3*threadIdx.x + k];
            }

            kSqVals_[threadIdx.x] = (kVec0*kVec0) + (kVec1*kVec1) +
                                    (kVec2*kVec2);
         }
         __syncthreads(); // wait for all threads to finish

         // Perform a parallel reduction on 128 threads to find min image
         // (note: stride >>= 1 is equivalent to stride /= 2)
         for (int stride = 64; stride > 0; stride >>= 1) {
            if (paramID < kSize) { // only evaluate if on the k-grid
               if (imageID < stride) {
                  bool swap = false;
                  cudaReal a = kSqVals_[threadIdx.x];
                  cudaReal b = kSqVals_[threadIdx.x+stride];
                  cudaReal tol = epsilon * (a > b ? a : b);
                  if (b < a - tol) {
                     swap = true;
                  } else if (fabs(b - a) <= tol) {
                     // kSq values effectively equal.
                     // Determine whether to swap based on hkl indices
                     for (int i = 0; i < 3; ++i) {
                        if (images_[3*threadIdx.x+i] >
                                     images_[3*(threadIdx.x+stride)+i]) {
                           break;
                        } else
                        if (images_[3*threadIdx.x+i] <
                                     images_[3*(threadIdx.x+stride)+i]) {
                           swap = true;
                           break;
                        }
                     }
                  }
                  if (swap) {
                     images_[3*threadIdx.x] =
                                         images_[3*(threadIdx.x+stride)];
                     images_[3*threadIdx.x+1] =
                                         images_[3*(threadIdx.x+stride)+1];
                     images_[3*threadIdx.x+2] =
                                         images_[3*(threadIdx.x+stride)+2];
                     kSqVals_[threadIdx.x] = kSqVals_[threadIdx.x+stride];
                  }
               }
            }
            __syncthreads(); // wait for all threads to finish
         }

         // At this point, for any thread with imageID == 0, the 
         // corresponding entries in kSqVals_ and images_ should 
         // contain the kSq value and the minimum image, respectively. 
         // Store values and exit
         if ((imageID == 0) && (paramID < kSize)) {
            int r[3];
            r[0] = images_[3*threadIdx.x];
            r[1] = images_[3*threadIdx.x+1];
            r[2] = images_[3*threadIdx.x+2];
            cudaReal kSqMin = kSqVals_[threadIdx.x];
            _descendToMinimum<3>(r, kSqMin, kBasis, shifts, epsilon);
            kSq[paramID] = kSqMin;
            minImages[paramID] = r[0];
            minImages[paramID + kSize] = r[1];
            minImages[paramID + kSize + kSize] = r[2];
         }
      }

      /*
      * Compute the kSq array on the GPU.
      */
      template <int D> __global__ 
      void _computeKSq(cudaReal* kSq, 
                       int const * waveBz,
                       cudaReal const * kBasis,
                       int const nParams, 
                       int const kSize)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;

         // Load kBasis into shared memory for fast access
         __shared__ cudaReal kBasis_[D*D];
         if (threadIdx.x < D * D) {
            kBasis_[threadIdx.x] = kBasis[threadIdx.x];
         }
         __syncthreads(); // wait for all threads to finish

         // Variables to be used in the loop
         int i, j, k, waveBz_;
         cudaReal kVec[D], kSqVal;
         // Note: usually local arrays are very slow in a CUDA kernel, but
         // the compiler should ideally be able to unroll the loops below
         // and store the kVec array in the register, because the full
         // structure of the loops below is known at compile time.

         // Loop through array
         for (i = startID; i < kSize; i += nThreads) {

            // initialize kVec to 0
            for (j = 0; j < D; ++j) {
               kVec[j] = 0.0;
            }

            // Calculate kVec
            for (j = 0; j < D; ++j) {
               waveBz_ = waveBz[i + (j * kSize)];
               for (k = 0; k < D; ++k) {
                  kVec[k] += kBasis_[k + (D*j)] * waveBz_;
               }
            }

            // Compute kSq
            kSqVal = 0.0;
            for (j = 0; j < D; ++j) {
               kSqVal += kVec[j] * kVec[j];
            }

            // Store value in global memory
            kSq[i] = kSqVal;

         } // kSize
      }

      /*
      * Derivative of |k|^2 of the image img with respect to parameter
      * param (dkkBasis as in _computedKSq), in floating point, because
      * the product of two wave indices can overflow int.
      */
      template <int D>
      __device__ cudaReal _dKSqImage(int const * img,
                                     cudaReal const * dkkBasis,
                                     int const param)
      {
         cudaReal value = 0.0;
         for (int j = 0; j < D; ++j) {
            for (int k = 0; k < D; ++k) {
               value += cudaReal(img[j]) * cudaReal(img[k])
                        * dkkBasis[k + (j * D) + (param * D * D)];
            }
         }
         return value;
      }

      /*
      * Compute the dKSq array on the GPU.
      *
      * The implicitInverse array is only used if isRealField is true.
      * If isRealField is false, implicitInverse may be a null pointer.
      * If averageTies is true, dKSq of a wavevector whose minimum image
      * is tied with other images, within the relative tolerance 
      * tieTolerance, is the average over these images, unless their 
      * derivatives agree within this tolerance (see tiedMinimumImages
      * and tiedImageDKSq in shiftToMinimum.h): they are found from the
      * stored image by steps of the Voronoi-relevant vectors in shifts,
      * using kBasis.
      *
      * When launched, this kernel requires dynamic memory allocation of
      * sizeof(cudaReal) * (nParams + 1) * D * D.
      */
      template <int D> __global__ 
      void _computedKSq(cudaReal* dKSq, 
                        int const * waveBz,
                        cudaReal const * dkkBasis,
                        bool const * implicitInverse,
                        int const nParams, 
                        int const kSize,
                        bool isRealField,
                        cudaReal const * kBasis,
                        ImageShifts const shifts,
                        bool const averageTies,
                        cudaReal const tieTolerance)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;

         // Load dkkBasis and kBasis into shared memory for fast access
         // (max size of dkkBasis is 54 elements for triclinic unit cell)
         extern __shared__ cudaReal dkkBasis_[];
         cudaReal* kBasis_ = &dkkBasis_[nParams * D * D];
         if (threadIdx.x < nParams * D * D) {
            dkkBasis_[threadIdx.x] = dkkBasis[threadIdx.x];
         }
         if (threadIdx.x < D * D) {
            kBasis_[threadIdx.x] = kBasis[threadIdx.x];
         }
         __syncthreads(); // wait for all threads to finish

         // Tied images (image t has components tied[D*t], ...)
         const int maxTied = 16;
         int tied[maxTied * D];
         int u[D];

         // Loop through array
         for (int i = startID; i < kSize; i += nThreads) {

            // Load waveBz to local memory
            for (int j = 0; j < D; ++j) {
               tied[j] = waveBz[i + (j * kSize)];
            }

            // Find tied images, as tiedMinimumImages
            int nTied = 1;
            if (averageTies) {
               cudaReal g0 = _kSqImage<D>(tied, kBasis_);
               for (int t = 0; t < nTied; ++t) {
                  for (int s = 0; s < shifts.n; ++s) {
                     for (int j = 0; j < D; ++j) {
                        u[j] = tied[D*t + j] + shifts.v[D*s + j];
                     }
                     cudaReal g = _kSqImage<D>(u, kBasis_);
                     cudaReal tol = tieTolerance * (g > g0 ? g : g0);
                     if (fabs(g - g0) > tol) continue;
                     bool found = false;
                     for (int m = 0; m < nTied; ++m) {
                        bool same = true;
                        for (int j = 0; j < D; ++j) {
                           if (tied[D*m + j] != u[j]) same = false;
                        }
                        if (same) found = true;
                     }
                     if (!found && nTied < maxTied) {
                        for (int j = 0; j < D; ++j) {
                           tied[D*nTied + j] = u[j];
                        }
                        ++nTied;
                     }
                  }
               }
            }

            for (int param = 0; param < nParams; ++param) {
               cudaReal dKSqVal = _dKSqImage<D>(tied, dkkBasis_, param);
               if (nTied > 1) {
                  // Average, unless all derivatives agree (as in 
                  // tiedImageDKSq in shiftToMinimum.cpp)
                  cudaReal sum = dKSqVal;
                  cudaReal maxDiff = 0.0;
                  cudaReal scale = fabs(dKSqVal);
                  for (int t = 1; t < nTied; ++t) {
                     cudaReal value 
                             = _dKSqImage<D>(&tied[D*t], dkkBasis_, param);
                     sum += value;
                     cudaReal diff = fabs(value - dKSqVal);
                     if (diff > maxDiff) maxDiff = diff;
                     if (fabs(value) > scale) scale = fabs(value);
                  }
                  if (maxDiff > tieTolerance*scale) {
                     dKSqVal = sum/cudaReal(nTied);
                  }
               }
               if (isRealField) {
                  if (implicitInverse[i]) {
                     // if element i's inverse is implicit
                     dKSqVal *= 2;
                  }
               }
               dKSq[(param * kSize) + i] = dKSqVal;
            } // nParams

         } // kSize
      }

   } // anonymous namespace Pscf::Prdc::(unnamed)

   // Member function definitions

   /*
   * Constructor.
   */
   template <int D>
   WaveList<D,CUT>::WaveList(bool isRealField)
    : kSize_(0),
      nBunch_(0),
      isAllocated_(false),
      hasMinImages_(false),
      hasMinImages_h_(false),
      hasKSq_(false),
      hasdKSq_(false),
      isSorted_(false),
      isRealField_(isRealField),
      nNyquist_(0),
      nAverage_(0),
      averageNyquist_(false),
      unitCellPtr_(nullptr),
      meshPtr_(nullptr)
   {}

   /*
   * Destructor.
   */
   template <int D>
   WaveList<D,CUT>::~WaveList()
   {

      if (dKSqSlices_.isAllocated()) {
         int nParams = dKSqSlices_.capacity();
         for (int i = 0; i < nParams; i++) {
            if (dKSqSlices_[i].isAssociated()) {
               dKSqSlices_[i].dissociate();
            }
         }
      }

      /*
      * The above loop dissociates elements of dKSqSlices_ from memory
      * owned by dKSq_. Because the destructor body is called before
      * destructors for members, this operation will occur before either
      * container is destroyed.  Strictly speaking, this is a redundant
      * safety measure as long as dKSqSlices_ appears in the declaration
      * of members after dKSq_, because this order causes dKSqSlices_
      * to to be destroyed before dKSq_, and the destructor for each
      * element of dKSqSlices_ would also call dissociate() if the
      * association still existed at that point.
      */

   }

   /*
   * Allocate memory.
   *
   * This function allocates and construct implicitInverse_d_ if and
   * only if isRealField_ == true.
   */
   template <int D>
   void WaveList<D,CUT>::allocate(Mesh<D> const & m, UnitCell<D> const & c)
   {
      UTIL_CHECK(m.size() > 0);
      UTIL_CHECK(c.nParameter() > 0);
      UTIL_CHECK(!isAllocated_);

      // Create permanent associations with mesh and unit cell
      unitCellPtr_ = &c;
      meshPtr_ = &m;

      int nParams = unitCell().nParameter();
      IntVec<D> const & meshDimensions = mesh().dimensions();

      // Compute kMeshDimensions_ and kSize_
      if (isRealField_) {
         FFT<D,CUT>::computeKMesh(meshDimensions, kMeshDimensions_, kSize_);
      } else {
         kMeshDimensions_ = meshDimensions;
         kSize_ = mesh().size();
      }
      UTIL_CHECK(kSize_ > 0);

      minImages_.allocate(kSize_ * D);
      kSq_.allocate(kMeshDimensions_);
      dKSq_.allocate(kSize_ * nParams);
      dKSqSlices_.allocate(nParams);
      for (int i = 0; i < nParams; i++) {
         dKSqSlices_[i].associate(dKSq_, i*kSize_, kMeshDimensions_);
      }

      // Allocate and fill implicitInverse arrays if isRealField_ == true
      // (only depends on mesh dimensions, only used for real fields)
      if (isRealField_) {
         implicitInverse_d_.allocate(kSize_);
         implicitInverse_h_.allocate(kSize_);
         HostSendArray<bool,CUT> implicitTemp;
         implicitTemp.allocate(kSize_);
         MeshIterator<D> kItr(kMeshDimensions_);
         int rank;
         for (kItr.begin(); !kItr.atEnd(); ++kItr) {
            rank = kItr.rank();
            implicitTemp[rank] = 
               FFT<D,CUT>::hasImplicitInverse(kItr.position(), 
                                              meshDimensions);
            implicitInverse_h_[rank] = implicitTemp[rank];
         }
         implicitInverse_d_ = implicitTemp; // copy to device memory
      }

      // For a lattice system with a variable angle, list wavevectors
      // that have a Nyquist index (see class documentation)
      averageNyquist_ = hasVariableAngle<D>(unitCell().lattice());
      nNyquist_ = 0;
      if (averageNyquist_) {
         MeshIterator<D> kItr(kMeshDimensions_);
         for (kItr.begin(); !kItr.atEnd(); ++kItr) {
            if (hasNyquistIndex(kItr.position(), meshDimensions)) {
               ++nNyquist_;
            }
         }
         if (nNyquist_ > 0) {
            nyquistWaves_h_.allocate(nNyquist_);
            nyquistRanks_h_.allocate(nNyquist_);
            int i = 0;
            for (kItr.begin(); !kItr.atEnd(); ++kItr) {
               if (hasNyquistIndex(kItr.position(), meshDimensions)) {
                  nyquistWaves_h_[i] = kItr.position();
                  nyquistRanks_h_[i] = kItr.rank();
                  ++i;
               }
            }
            UTIL_CHECK(i == nNyquist_);
            averageRanks_h_.allocate(nNyquist_);
            averageRanks_d_.allocate(nNyquist_);
            averageImages_h_.allocate(D*nNyquist_);
            averageImages_d_.allocate(D*nNyquist_);
            averageKSq_h_.allocate(nNyquist_);
            averageKSq_d_.allocate(nNyquist_);
            averageDKSq_h_.allocate(nParams*nNyquist_);
            averageDKSq_d_.allocate(nParams*nNyquist_);
         }
      }

      clearUnitCellData();
      isAllocated_ = true;
   }

   /*
   * Clear data that depends on unit cell parameters.
   */
   template <int D>
   void WaveList<D,CUT>::clearUnitCellData()
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
         hasMinImages_h_ = false;
      }
   }

   /*
   * Compute minimum image vectors and kSq.
   */
   template <int D>
   void WaveList<D,CUT>::computeMinimumImages()
   {
      // If min images are valid, return immediately
      if (hasMinImages_) return;

      // Precondition
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().isInitialized());
      UTIL_CHECK(minImages_.capacity() == kSize_ * D);

      // Set initial array of images to contain the k-grid points
      HostSendArray<int,CUT> imagesTmp(D*kSize_);
      MeshIterator<D> kItr(kMeshDimensions_);
      for (int i = 0; i < D; i++) {
         for (kItr.begin(); !kItr.atEnd(); ++kItr) {
            imagesTmp[kItr.rank() + (i*kSize_)] = kItr.position(i);
         }
      }
      minImages_ = imagesTmp; // copy to device

      // Get kBasis and meshDims and store on device
      HostSendArray<cudaReal,CUT> kBasis_h(D*D);
      HostSendArray<int,CUT> meshDims_h(D);
      int idx = 0;
      for (int j = 0; j < D; ++j) {
         for (int k = 0; k < D; ++k) {
            kBasis_h[idx] = unitCell().kBasis(j)[k];
            idx++;
         }
         meshDims_h[j] = mesh().dimension(j);
      }

      BackendArray<cudaReal,CUT> kBasis(D*D);
      BackendArray<int,CUT> meshDims(D);
      kBasis = kBasis_h;
      meshDims = meshDims_h;

      // Set number of threads per gridpoint (depends on D)
      int threadsPerGP;
      if (D == 3) {
         threadsPerGP = 128;
      } else if (D == 2) {
         threadsPerGP = 32;
      } else if (D == 1) {
         threadsPerGP = 1;
      }

      // The kernel uses threadsPerGP threads per wavevector. For a large
      // mesh, kSize_*threadsPerGP exceeds the range of int (e.g., in 3D
      // for kSize_ >= 2^24). The kernel is therefore launched for chunks
      // of at most kChunk wavevectors, with at most 2^30 threads each.
      const int kChunkMax = (1 << 30) / threadsPerGP;
      const int kChunk = (kSize_ < kChunkMax) ? kSize_ : kChunkMax;

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(kChunk*threadsPerGP, nBlocks, nThreads);

      if ((D == 3) && (nThreads < 128)) {
         // Thread blocks too small. Manually set nThreads to 128
         ThreadArray::setThreadsPerBlock(128);
         ThreadArray::setThreadsLogical(kChunk*threadsPerGP, nBlocks, nThreads);

         // If the above was successful, print warning
         Log::file() << "\nWarning: "
                     << "nThreads too small for computeMinimumImages.\n"
                     << "Setting nThreads equal to 128." << std::endl;
      }

      UTIL_CHECK(kSq_.isAllocated());
      UTIL_CHECK(minImages_.isAllocated());

      // Voronoi-relevant vectors of the lattice of image shifts
      FSArray<IntVec<D>, 14> neighbors;
      imageShiftNeighbors<D>(mesh().dimensions(), unitCell(), neighbors);
      ImageShifts shifts;
      shifts.n = neighbors.size();
      for (int j = 0; j < neighbors.size(); ++j) {
         for (int i = 0; i < D; ++i) {
            shifts.v[D*j + i] = neighbors[j][i];
         }
      }

      // Launch kernel, once per chunk of wavevectors
      size_t sz = (D * sizeof(int) + sizeof(cudaReal)) * nThreads;
      for (int kBegin = 0; kBegin < kSize_; kBegin += kChunk) {
         _computeMinimumImages<D><<<nBlocks, nThreads, sz>>>
            (minImages_.cArray(), kSq_.cArray(), kBasis.cArray(),
             meshDims.cArray(), kSize_, kBegin, shifts);
         cudaErrorCheck( cudaGetLastError() );
      }

      // Wavevectors with a Nyquist index for which ksqNyquist averages
      // |k|^2 (see usesNyquistAverage): kSq, image and dKSq
      // from ksqNyquist on the host, in one evaluation, copied to the
      // device; dKSq is kept for computedKSq. Otherwise ksqNyquist is
      // the minimum-image |k|^2 computed by the kernels.
      nAverage_ = 0;
      if (averageNyquist_ && nNyquist_ > 0) {
         IntVec<D> const & d = mesh().dimensions();
         int nParams = unitCell().nParameter();
         double dk[6];
         IntVec<D> image;
         for (int i = 0; i < nNyquist_; ++i) {
            IntVec<D> const & v = nyquistWaves_h_[i];
            if (!usesNyquistAverage(v, d, unitCell())) continue;
            int k = nAverage_;
            int rank = nyquistRanks_h_[i];
            averageRanks_h_[k] = rank;
            averageKSq_h_[k] = ksqNyquist(v, d, unitCell(), dk, image);
            for (int j = 0; j < D; ++j) {
               averageImages_h_[k + j*nNyquist_] = image[j];
            }
            double factor = 1.0;
            if (isRealField_ && implicitInverse_h_[rank]) factor = 2.0;
            for (int n = 0; n < nParams; ++n) {
               averageDKSq_h_[k + n*nNyquist_] = factor*dk[n];
            }
            ++nAverage_;
         }
         if (nAverage_ > 0) {
            averageRanks_d_ = averageRanks_h_; // copy to device memory
            averageImages_d_ = averageImages_h_;
            averageKSq_d_ = averageKSq_h_;
            averageDKSq_d_ = averageDKSq_h_;
            int nBlocksN, nThreadsN;
            ThreadArray::setThreadsLogical(nAverage_, nBlocksN, nThreadsN);
            _setNyquistImages<D><<<nBlocksN, nThreadsN>>>
               (minImages_.cArray(), kSq_.cArray(),
                averageRanks_d_.cArray(), averageImages_d_.cArray(),
                averageKSq_d_.cArray(), nAverage_, nNyquist_, kSize_);
            cudaErrorCheck( cudaGetLastError() );
         }
      }

      hasMinImages_ = true;
      hasKSq_ = true;
   }

   /*
   * Compute values of k^2, using existing minImages if possible.
   */
   template <int D>
   void WaveList<D,CUT>::computeKSq()
   {
      // If kSq values are valid, return immediately without recomputing
      if (hasKSq_) return;

      // If necessary, compute minimum images
      if (!hasMinImages_) {
         computeMinimumImages(); // compute both min images and kSq
         return;
      }

      // If this point is reached, calculate kSq using _computeKSq kernel

      // Preconditions
      UTIL_CHECK(unitCell().nParameter() > 0);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().isInitialized());
      UTIL_CHECK(isAllocated_);

      // Wavevectors whose kSq is given by ksqNyquist exist only for
      // lattice systems with a variable angle, for which
      // clearUnitCellData also clears the minimum images: their kSq is
      // computed by computeMinimumImages
      UTIL_CHECK(nAverage_ == 0);

      // Get kBasis and store on device
      HostSendArray<cudaReal,CUT> kBasis_h(D*D);
      BackendArray<cudaReal,CUT> kBasis(D*D);
      int idx = 0;
      for (int j = 0; j < D; ++j) {
         for (int k = 0; k < D; ++k) {
            kBasis_h[idx] = unitCell().kBasis(j)[k];
            idx++;
         }
      }
      kBasis = kBasis_h;

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(kSize_, nBlocks, nThreads);

      // Launch kernel to calculate kSq on device
      UTIL_CHECK(kSq_.isAllocated());
      UTIL_CHECK(minImages_.isAllocated());
      _computeKSq<D><<<nBlocks, nThreads>>>
            (kSq_.cArray(), minImages_.cArray(), kBasis.cArray(),
             unitCell().nParameter(), kSize_);
      cudaErrorCheck( cudaGetLastError() );

      hasKSq_ = true;
   }

   /*
   * Compute derivatives of k^2 with respect to unit cell parameters.
   */
   template <int D>
   void WaveList<D,CUT>::computedKSq()
   {
      if (hasdKSq_) return; // dKSq already calculated

      // Compute minimum images if needed
      if (!hasMinImages_) {
         computeMinimumImages();
      }

      // Preconditions
      UTIL_CHECK(unitCell().nParameter() > 0);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().isInitialized());
      UTIL_CHECK(isAllocated_);

      // Calculate dkkBasis and store on device
      int idx;
      HostSendArray<cudaReal,CUT> dkkBasis_h(unitCell().nParameter() * D * D);
      BackendArray<cudaReal,CUT> dkkBasis;
      for (int i = 0 ; i < unitCell().nParameter(); ++i) {
         for (int j = 0; j < D; ++j) {
            for (int k = 0; k < D; ++k) {
               idx = k + (j * D) + (i * D * D);
               dkkBasis_h[idx] = unitCell().dkkBasis(i, j, k);
            }
         }
      }
      dkkBasis = dkkBasis_h;

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(kSize_, nBlocks, nThreads);

      // Kernel requires block size to be >= the size of dkkBasis.
      // Max size of dkkBasis is 54, so this should always be satisfied
      UTIL_CHECK(nThreads > dkkBasis.capacity());

      // Preconditions for kernel
      UTIL_CHECK(dKSq_.isAllocated());
      UTIL_CHECK(minImages_.isAllocated());
      bool const * implicitInversePtr = nullptr;
      if (isRealField_) {
         UTIL_CHECK(implicitInverse_d_.isAllocated());
         implicitInversePtr = implicitInverse_d_.cArray();
      }

      // Data to average dKSq over tied minimum images (variable angles
      // only: for other lattice systems, tied images have equal dKSq)
      bool averageTies = hasVariableAngle<D>(unitCell().lattice());
      HostSendArray<cudaReal,CUT> kBasis_h(D*D);
      BackendArray<cudaReal,CUT> kBasis(D*D);
      for (int j = 0; j < D; ++j) {
         for (int k = 0; k < D; ++k) {
            kBasis_h[D*j + k] = unitCell().kBasis(j)[k];
         }
      }
      kBasis = kBasis_h;
      ImageShifts shifts;
      shifts.n = 0;
      if (averageTies) {
         FSArray<IntVec<D>, 14> neighbors;
         imageShiftNeighbors<D>(mesh().dimensions(), unitCell(), neighbors);
         shifts.n = neighbors.size();
         for (int j = 0; j < neighbors.size(); ++j) {
            for (int i = 0; i < D; ++i) {
               shifts.v[D*j + i] = neighbors[j][i];
            }
         }
      }
      #ifdef SINGLE_PRECISION
      cudaReal tieTolerance = 1.0E-4;
      #else
      cudaReal tieTolerance = TiedImageTolerance;
      #endif

      // Launch kernel to calculate dKSq on device
      size_t sz = sizeof(cudaReal)*(dkkBasis.capacity() + D*D);
      _computedKSq<D><<<nBlocks, nThreads, sz>>>
         (dKSq_.cArray(), minImages_.cArray(), dkkBasis.cArray(),
          implicitInversePtr, unitCell().nParameter(), kSize_,
          isRealField_, kBasis.cArray(), shifts, averageTies, 
          tieTolerance);
      cudaErrorCheck( cudaGetLastError() );

      // Wavevectors with a Nyquist index whose kSq is given by 
      // ksqNyquist: derivatives from computeMinimumImages (same cell)
      if (nAverage_ > 0) {
         int nBlocksN, nThreadsN;
         ThreadArray::setThreadsLogical(nAverage_, nBlocksN, nThreadsN);
         _setNyquistValues<<<nBlocksN, nThreadsN>>>
            (dKSq_.cArray(), averageRanks_d_.cArray(), 
             averageDKSq_d_.cArray(), nAverage_, nNyquist_, 
             unitCell().nParameter(), kSize_);
         cudaErrorCheck( cudaGetLastError() );
      }

      hasdKSq_ = true;
   }

   /*
   * Sort waves by magnitude.
   */
   template <int D>
   void WaveList<D,CUT>::sortWaves()
   {
      // If waves are already sorted, return immediately
      if (isSorted_) return;

      // Compute wavenumbers if necessary
      UTIL_CHECK(isAllocated_);
      if (!hasKSq_) {
         computeKSq();
      }

      // Copy values of kSq to host
      HostRecvArray<cudaReal,CUT> kSq_h(kSq_);

      // Construct Sort::Item objects with value = kSq, id = wave id
      std::vector< Sort::Item<double> > items;
      Sort::Item<double> item;
      for (int i = 0; i < kSize_; ++i) {
         item.value = kSq_h[i];
         item.id = i;
         items.push_back(item);
      }
      UTIL_CHECK((int)items.size() == kSize_);

      // Sort array of items
      Sort::sort(items);

      // Fill sortedIds_ member array with sorted ids
      if (!sortedIds_.isAllocated()) {
         sortedIds_.allocate(kSize_);
      }
      UTIL_CHECK(sortedIds_.capacity() == kSize_);
      for (int i = 0; i < kSize_; ++i) {
         sortedIds_[i] = items[i].id;
      }

      // Construct sortedBunches_ array and set nBunch_
      sortedBunches_.clear();
      findWaveBunches(items, sortedBunches_, waveNormTolerance<cudaReal>());
      nBunch_ = sortedBunches_.size();

      // Construct bunchIds_ array
      UTIL_CHECK(nBunch_ > 0);
      if (!bunchIds_.isAllocated()) {
         bunchIds_.allocate(kSize_);
      }
      UTIL_CHECK(bunchIds_.capacity() == kSize_);
      int begin, end, ib, iw;
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

   /*
   * Get minimum images on the host by reference.
   */
   template <int D>
   DArray< IntVec<D> > const & 
   WaveList<D,CUT>::minImages_h() const
   {
      UTIL_CHECK(hasMinImages_);
      if (!hasMinImages_h_) {
         if (!minImages_h_.isAllocated()) {
            minImages_h_.allocate(kSize_);
         }
         HostRecvArray<int,CUT> minImages_temp;
         minImages_temp = minImages_;
         int offset;
         for (int j = 0; j < D; ++j) {
            offset = j*kSize_;
            for (int i = 0; i < kSize_; ++i) {
               minImages_h_[i][j] = minImages_temp[i + offset];
            }
         }
         hasMinImages_h_ = true;
      }
      return minImages_h_;
   }

} // namespace Pscf
#endif
