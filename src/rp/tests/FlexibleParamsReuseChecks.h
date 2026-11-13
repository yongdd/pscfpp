#ifndef RP_FLEXIBLE_PARAMS_REUSE_CHECKS_H
#define RP_FLEXIBLE_PARAMS_REUSE_CHECKS_H

#include <util/tests/LogFileUnitTest.h>
#include <rp/system/System.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <rp/field/Domain.h>
#include <rp/scft/iterator/Iterator.h>
#include <pscf/chem/PolymerModel.h>
#include <pscf/field/RField.h>
#include <pscf/backend/HostSendArray.h>
#include <pscf/backend/Reduce.h>
#include <util/param/BracketPolicy.h>
#include <sstream>
#include <cmath>

namespace Pscf {
namespace Rp {

template <class T>
class FlexibleParamsReuseChecks : public LogFileUnitTest
{
public:
   void setUp()
   {
      bracketPolicy_ = Util::BracketPolicy::get();
      Util::BracketPolicy::set(Util::BracketPolicy::Optional);
   }

   void tearDown()
   {
      LogFileUnitTest::tearDown();
      PolymerModel::setModel(PolymerModel::Thread);
      Util::BracketPolicy::set(bracketPolicy_);
   }

   void check(bool basis)
   {
      this->openLogFile(basis ? "out/FlexibleParamsReuseBasis.log"
                              : "out/FlexibleParamsReuseGrid.log");
      PolymerModel::setModel(PolymerModel::Thread);
      for (bool open : {false, true}) {
      for (bool initial : {false, true}) {
         System<1,T> system;
         std::ostringstream text;
         text << "System{ Mixture{ nMonomer 2 monomers[ 1 1 ] nPolymer 1 "
              << "Polymer{ type linear nBlock 2 blocks[ 0 .5 1 .5 ] "
              << (open ? "mu 5 } " : "phi 1 } ")
              << "ds .05 } Interaction{ chi( 0 1 10 ) } "
              << "Domain{ mesh 8 lattice lamellar groupName P_-1 } "
              << (basis ? "AmIteratorBasis" : "AmIteratorGrid")
              << "{ epsilon 1e-10 maxItr 100 maxHist 4 errorType "
              << (open ? "relNorm " : "maxResid ")
              << "isFlexible " << int(initial) << " } }";
         std::istringstream in(text.str());
         system.readParam(in);
         Util::FSArray<double,6> cell; cell.append(1.0);
         system.setUnitCell(cell);
         Util::DArray<RField<1,T> > fields; fields.allocate(2);
         IntVec<1> dims(8);
         for (int m = 0; m < 2; ++m) fields[m].allocate(dims);

         // Repeatedly grow and shrink the vectors, including continuation
         // calls and unchanged flags. Perturbed starts exercise populated
         // history/basis slots, not just the first vector of an exact solve.
         bool flags[] = {!initial, initial, !initial, !initial, initial};
         for (int step = 0; step < 5; ++step) {
            for (int m = 0; m < 2; ++m) {
               HostSendArray<typename T::Real,T> host;
               host.associate(fields[m]);
               for (int i = 0; i < 8; ++i) {
                  host[i] = 5.0 + (m ? -.02 : .02)
                                  *std::cos(2.0*std::acos(-1.0)*i/8);
               }
               fields[m] = host;
               host.dissociate();
            }
            system.w().setRGrid(fields, basis);
            Util::FSArray<bool,6> flex; flex.append(flags[step]);
            system.iterator().setFlexibleParams(flex);
            TEST_ASSERT(system.iterate(step % 2 == 1) == 0);
            for (int m = 0; m < 2; ++m) {
               const auto& c = system.c().rgrid(m);
               TEST_ASSERT(std::abs(Reduce::min(c) - .5) < 1e-8);
               TEST_ASSERT(std::abs(Reduce::max(c) - .5) < 1e-8);
            }
         }
      }
      }
   }

private:
   Util::BracketPolicy::Type bracketPolicy_;
};

}
}
#endif
