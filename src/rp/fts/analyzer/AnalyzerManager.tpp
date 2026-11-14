#ifndef RP_ANALYZER_MANAGER_TPP
#define RP_ANALYZER_MANAGER_TPP

#include <rp/fts/analyzer/AnalyzerManager.h>
#include <rp/fts/analyzer/AnalyzerFactory.h>

#include <util/param/Factory.h>
#include <util/param/ParamComposite.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   AnalyzerManager<D,T>::AnalyzerManager(Simulator<D,T>& simulator,
                                         System<D,T>& system)
    : Base(),
      simulatorPtr_(&simulator),
      systemPtr_(&system)
   {  ParamComposite::setClassName("AnalyzerManager"); }

   /*
   * Return a pointer to a new AnalyzerFactory object.
   */
   template <int D, class T>
   Factory<Analyzer<D,T> >* AnalyzerManager<D,T>::newDefaultFactory() 
   const
   {  return new AnalyzerFactory<D,T>(*simulatorPtr_, *systemPtr_); }

   /*
   * Read body of parameter file block.
   */
   template <int D, class T>
   void AnalyzerManager<D,T>::readParameters(std::istream &in)
   {
      Analyzer<D,T>::baseInterval = 1;
      ParamComposite::readOptional(in, "baseInterval",
                                   Analyzer<D,T>::baseInterval);
      Base::readParameters(in);
   }

   /*
   * Call setup method of each analyzer.
   */
   template <int D, class T>
   void AnalyzerManager<D,T>::setup()
   {
      for (int i = 0; i < Base::size(); ++i) {
         (*this)[i].setup();
      }
   }

   /*
   * Call sample method of each analyzer.
   */
   template <int D, class T>
   void AnalyzerManager<D,T>::sample(long iStep)
   {
      int baseInterval = Analyzer<D,T>::baseInterval;
      UTIL_CHECK(baseInterval > 0);
      UTIL_CHECK(iStep % baseInterval == 0);
      for (int i = 0; i < Base::size(); ++i) {
         (*this)[i].sample(iStep);
      }
   }

   /*
   * Call output method of each analyzer.
   */
   template <int D, class T>
   void AnalyzerManager<D,T>::output()
   {
      for (int i = 0; i < Base::size(); ++i) {
         (*this)[i].output();
      }
   }

   /*
   * Check that no analyzer requires a fixed unit cell.
   */
   template <int D, class T>
   void AnalyzerManager<D,T>::checkVariableUnitCell() const
   {
      for (int i = 0; i < Base::size(); ++i) {
         if ((*this)[i].requiresFixedUnitCell()) {
            std::string message = "Analyzer ";
            message += (*this)[i].className();
            message += " requires a fixed unit cell and cannot be used"
                       " with a ramp that changes the unit cell"
                       " (cell_param)";
            UTIL_THROW(message.c_str());
         }
      }
   }

}
}
#endif
