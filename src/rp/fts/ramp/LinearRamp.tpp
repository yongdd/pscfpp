#ifndef RP_LINEAR_RAMP_TPP
#define RP_LINEAR_RAMP_TPP

#include "LinearRamp.h"
#include <rp/fts/simulator/Simulator.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/MixtureModifier.h>
#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   LinearRamp<D,T>::LinearRamp(Simulator<D,T>& simulator)
    : RampT(simulator)
   {  ParamComposite::setClassName("LinearRamp"); }

   template <int D, class T>
   void LinearRamp<D,T>::readParameters(std::istream& in)
   {
      // Read the number of ramp parameters, allocate parameters_ array
      ParamComposite::read(in, "nParameter", nParameter_);
      UTIL_CHECK(nParameter_ > 0);
      parameters_.allocate(nParameter_);

      // Read array of RampParameters, using << operator for each
      ParamComposite::readDArray(in, "parameters", parameters_, 
                                 nParameter_);

      // Verify net zero change in volume fractions, if these are swept
      // and every species is closed (open species take up the rest of
      // the volume, so closed volume fractions may then change freely)
      if (Ramp<D,T>::simulator().system().mixture().isCanonical()) {
         double sum = 0.0;
         for (int i = 0; i < nParameter_; ++i) {
            if (parameters_[i].type() == "phi_polymer" ||
                parameters_[i].type() == "phi_solvent")
            {
               sum += parameters_[i].change();
            }
         }
         UTIL_CHECK(sum > -0.0000001);
         UTIL_CHECK(sum < 0.0000001);
      }
   }

   template <int D, class T>
   void LinearRamp<D,T>::setup(int nStep)
   {
      Ramp<D,T>::nStep_ = nStep;

      // Set simulator pointer and initial value for each parameter object
      for (int i = 0; i < nParameter_; ++i) {
         parameters_[i].setSimulator(Ramp<D,T>::simulator());
         parameters_[i].getInitial();
      }
   }

   template <int D, class T>
   void LinearRamp<D,T>::setParameters(int iStep)
   {
      // Compute a ramp parameter in the range [0,1]
      double s = double(iStep)/double(Ramp<D,T>::nStep_);

      // Update the system parameter values
      double newVal;
      bool hasPhi = false;
      bool hasChi = false;
      for (int i = 0; i < nParameter_; ++i) {
         newVal = parameters_[i].initial() + s*parameters_[i].change();
         parameters_[i].update(newVal);

         // A simultaneous chi ramp is valid as a complete matrix, not
         // necessarily after each individual entry has been changed.
         if (parameters_[i].type() == "chi") {
            hasChi = true;
         }

         if (parameters_[i].type() == "phi_polymer" ||
             parameters_[i].type() == "phi_solvent")
         {  hasPhi = true; }

      }

      if (hasChi) {
         Ramp<D,T>::simulator().analyzeChi();
      }

      // If every species is closed, the changes in volume fractions
      // cancel only to within the tolerance of readParameters. Rescale
      // the volume fractions to add up to 1, as on input: the FTS
      // compressors could not otherwise remove the k = 0 residual.
      System<D,T>& system = Ramp<D,T>::simulator().system();
      if (hasPhi && system.mixture().isCanonical()) {
         system.mixtureModifier().normalizePhi();
      }
   }

   template <int D, class T>
   bool LinearRamp<D,T>::changesUnitCell() const
   {
      for (int i = 0; i < nParameter_; ++i) {
         if (parameters_[i].type() == "cell_param" &&
             parameters_[i].change() != 0.0)
         {  return true; }
      }
      return false;
   }

   template <int D, class T>
   void LinearRamp<D,T>::output()
   {
      for (int i = 0; i < nParameter_; ++i) {
         Log::file() << "Parameter type: "
                     << parameters_[i].type() << std::endl;
         Log::file() << "Initial: "
                     << parameters_[i].initial() << " ";
         Log::file()<< "Change: "
                    << parameters_[i].change() << std::endl;
      }
   }

}
}
#endif
