/*
 * PSCF - Polymer Self-Consistent Field
 *
 * Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
 * Distributed under the terms of the GNU General Public License.
 */

#include "BinaryRelaxIterator.h"
#include <r1/system/System.h>
#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <r1/solvers/Solvent.h>
#include <pscf/chem/Ensemble.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/iterator/NanException.h>
#include <util/misc/Log.h>
#include <util/misc/Timer.h>
#include <math.h>
#include <limits>

namespace Pscf {
namespace R1
{
   using namespace Util;

   // Constructor
   BinaryRelaxIterator::BinaryRelaxIterator(System& system)
    : Iterator(system),
      epsilon_(0.0),
      lambdaPlus_(0.0),
      lambdaMinus_(0.0),
      lambdaPlusWork_(0.0),
      lambdaMinusWork_(0.0),
      maxItr_(200),
      isAllocated_(false),
      isCanonical_(true)
   {  setClassName("BinaryRelaxIterator"); }

   // Destructor
   BinaryRelaxIterator::~BinaryRelaxIterator()
   {}

   // Read parameter file block and allocate memory
   void BinaryRelaxIterator::readParameters(std::istream& in)
   {
      UTIL_CHECK(domain().nx() > 0);
      read(in, "epsilon", epsilon_);
      read(in, "maxIter", maxItr_);
      read(in, "lambdaPlus", lambdaPlus_);
      read(in, "lambdaMinus", lambdaMinus_);

      allocate();
   }

   // Allocate required memory (called by readParameters)
   void BinaryRelaxIterator::allocate()
   {
      if (isAllocated_) return;
      int nm = mixture().nMonomer();  // number of monomer types
      UTIL_CHECK(nm == 2);
      int nx = domain().nx();         // number of grid points
      UTIL_CHECK(nx > 0);
      cArray_.allocate(nm);
      wArray_.allocate(nm);
      dW_.allocate(nm);
      dWNew_.allocate(nm);
      wFieldsNew_.allocate(nm);
      cFieldsNew_.allocate(nm);
      for (int i = 0; i < nm; ++i) {
         wFieldsNew_[i].allocate(nx);
         cFieldsNew_[i].allocate(nx);
         dW_[i].allocate(nx);
         dWNew_[i].allocate(nx);
      }
      isAllocated_ = true;
   }

   void BinaryRelaxIterator::computeDW(Array<FieldT> const & wOld,
                              Array<FieldT> const & cFields,
                              Array<FieldT> & dW,
                              double & dWNorm)

   {
      double dWm, dWp, c0, c1, w0, w1, wm;
      double dWmNorm = 0.0;
      double dWpNorm = 0.0;

      // The SCF equations require w1 - w0 = (chi01 - chi00)*c0
      // + (chi11 - chi01)*c1. The exchange residual is scaled by the
      // effective chi = chi01 - (chi00 + chi11)/2. With zero diagonal
      // elements, r0 = r1 = 1 and chiEff = chi01.
      Interaction const & interaction = system().interaction();
      const double chi01 = interaction.chi(0,1);
      const double chiEff = chi01
                  - 0.5*(interaction.chi(0,0) + interaction.chi(1,1));
      if (chiEff == 0.0) {
         UTIL_THROW("BinaryRelaxIterator requires a nonzero effective chi");
      }
      const double r0 = (chi01 - interaction.chi(0,0))/chiEff;
      const double r1 = (chi01 - interaction.chi(1,1))/chiEff;
      int nx = domain().nx();        // number of grid points
      //For AB diblok
      for (int i = 0; i < nx; ++i) {
         c0 = cFields[0][i];
         c1 = cFields[1][i];
         w0 = wOld[0][i];
         w1 = wOld[1][i];
         wm = w1 - w0;
         dWp = c1 + c0 - 1.0;
         dWm = 0.5*( r0*c0 - r1*c1 - wm/chiEff);
         dWpNorm += dWp*dWp;
         dWmNorm += dWm*dWm;
         dWp *= lambdaPlusWork_;
         dWm *= lambdaMinusWork_;
         dW[0][i] = dWp -  dWm;
         dW[1][i] = dWp + dWm;
      }
      dWNorm = (dWpNorm + dWmNorm)/double(nx);
      dWNorm = sqrt(dWNorm);
   }

   void BinaryRelaxIterator::updateWFields(Array<FieldT> const & wOld,
                                  Array<FieldT> const & dW,
                                  Array<FieldT> & wNew)
   {
      int nm = mixture().nMonomer();  // number of monomer types
      UTIL_CHECK(nm == 2);
      int nx = domain().nx();        // number of grid points
      int i;                         // monomer index
      int j;                         // grid point index
      double w0, w1;

      // Add dW
      for (j = 0; j < nx; j++){
         w0 = wOld[0][j];
         w1 = wOld[1][j];
         wNew[0][j] = w0 + dW[0][j];
         wNew[1][j] = w1 + dW[1][j];
       }

      // If canonical, shift such that last element is exactly zero
      if (isCanonical_) {
         double shift = wNew[nm-1][nx-1];
         for (i = 0; i < nm; ++i) {
            for (j = 0; j < nx; ++j) {
               wNew[i][j] -= shift;
            }
         }
      }

   }

   /*
   * Solve the MDE for given fields and compute the increments.
   */
   bool BinaryRelaxIterator::computeTrial(DArray<FieldT> const & wFields,
                                          DArray<FieldT> & cFields,
                                          Array<FieldT> & dW,
                                          double & dWNorm)
   {
      try {
         mixture().compute(wFields, cFields);
      } catch (NanException const &) {
         // The exception has written its message to the log file
         dWNorm = std::numeric_limits<double>::infinity();
         return false;
      }
      computeDW(wFields, cFields, dW, dWNorm);
      return true;
   }

   int BinaryRelaxIterator::solve(bool isContinuation)
   {
      // Start overall timer
      Timer timerTotal;
      timerTotal.start();

      int nm = mixture().nMonomer();  // number of monomer types
      UTIL_CHECK(nm == 2);

      // Start each solve with the mixing parameters from the parameter
      // file; only the working copies are reduced below, so a reduction
      // does not carry over to later solves (or into WRITE_PARAM)
      lambdaPlusWork_ = lambdaPlus_;
      lambdaMinusWork_ = lambdaMinus_;
      int np = mixture().nPolymer();  // number of polymer species
      int ns = mixture().nSolvent();  // number of solvent species
      int nx = domain().nx();         // number of grid points

      // Determine if isCanonical (iff all species ensembles are closed)
      isCanonical_ = true;
      Ensemble ensemble;
      for (int i = 0; i < np; ++i) {
         ensemble = mixture().polymer(i).ensemble();
         if (ensemble == Ensemble::Unknown) {
            UTIL_THROW("Unknown species ensemble");
         }
         if (ensemble == Ensemble::Open) {
            isCanonical_ = false;
         }
      }
      if (ns > 0) {
         for (int i = 0; i < ns; ++i) {
            ensemble = system().mixture().solvent(i).ensemble();
            if (ensemble == Ensemble::Unknown) {
               UTIL_THROW("Unknown species ensemble");
            }
            if (ensemble == Ensemble::Open) {
               isCanonical_ = false;
            }
         }
      }

      // If isCanonical, shift so that last element is zero.
      // Note: This is one of the residuals in this case.
      if (isCanonical_) {
         double shift = wFields()[nm-1][nx-1];
         int i, j;
         for (i = 0; i < nm; ++i) {
            for (j = 0; j < nx; ++j) {
               wFields()[i][j] -= shift;
            }
         }
      }

      // Compute initial dWNorm.
      if (!computeTrial(system().wFields(), system().cFields(), dW_,
                        dWNorm_)) {
         Log::file() << "Non-finite initial residual" << std::endl;
         return 1;
      }

      // Iterative loop
      int i, j, k;
      for (i = 0; i < maxItr_; ++i) {
         Log::file() << "iteration " << i
         << " , error = " << dWNorm_
         << std::endl;

         if (dWNorm_ < epsilon_) {
            // Stop timers
            timerTotal.stop();

            Log::file() << "The epsilon is " << epsilon_<< std::endl;
            Log::file() << "Converged" << std::endl;
            system().computeFreeEnergy();
            // Success
            Log::file() << "\n\n";
            // Output timing resultsl;
            Log::file() << "Total time:                             "
                        << timerTotal.time()   << " s  "  << std::endl;
            Log::file() << "Average time cost of each iteration:    "
                        << timerTotal.time()/i  << " s  " << std::endl;
            Log::file() << "\n\n";

            // Success (iterated to convergence)
            return 0;
         }

         // Try full Fd relaxation update
         updateWFields(system().wFields(), dW_, wFieldsNew_);
         computeTrial(wFieldsNew_, cFieldsNew_, dWNew_, dWNormNew_);

         // Decrease increment if necessary
         j = 0;
         while (dWNormNew_ > dWNorm_ && j < 3) {
            Log::file() << "      error = "
                        << dWNormNew_ << ", decreasing increment" << std::endl;
            lambdaPlusWork_ *= 0.5;
            lambdaMinusWork_ *= 0.5;
            //Print the reduced mixing parameters
            Log::file() << "      lambdaPlus = "
                        << lambdaPlusWork_ << std::endl;
            Log::file() << "      lambdaMinus = "
                        << lambdaMinusWork_<< std::endl;
            computeDW(system().wFields(), system().cFields(),
                      dWNew_, dWNormNew_);
            updateWFields(system().wFields(), dWNew_, wFieldsNew_);
            computeTrial(wFieldsNew_, cFieldsNew_, dWNew_, dWNormNew_);
            ++j;
         }

         // Accept update or declare failure and break iteration loop
         if (dWNormNew_ < dWNorm_) {
            // Update system fields
            for (j = 0; j < nm; ++j) {
               for (k = 0; k < nx; ++k) {
                  system().wField(j)[k] = wFieldsNew_[j][k];
                  system().cField(j)[k] = cFieldsNew_[j][k];
                  dW_[j][k] = dWNew_[j][k];
               }
            }
            dWNorm_ = dWNormNew_;
         } else {
            Log::file() << "Iteration failed, norm = "
                        << dWNormNew_ << std::endl;
            break;
         }

      }

      // Failure (no convergence after maxIter iterations)
      return 1;
   }


} // namespace R1
} // namespace Pscf
