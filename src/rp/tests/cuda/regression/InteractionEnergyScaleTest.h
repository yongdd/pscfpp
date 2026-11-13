#ifndef RPG_INTERACTION_ENERGY_SCALE_TEST_H
#define RPG_INTERACTION_ENERGY_SCALE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <rp/system/System.h>
#include <rp/scft/ScftThermo.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <pscf/crystal/Basis.h>
#include <pscf/field/RField.h>
#include <pscf/backend/VecOp.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/floryHuggins/FhMixture.tpp>
#include <util/tests/LogFileUnitTest.h>
#include <util/misc/Exception.h>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace Util;
using namespace Pscf;

class InteractionEnergyScaleTest : public LogFileUnitTest
{
public:
   void setUp() { setVerbose(0); }

   void testInteractionScaleThermo()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/InteractionEnergyScaleInteractionScale.log");
      const double chis[] = {0.0, 10.0, 1e-9, -1e-9, 1.01e-9};
      for (double chi : chis) {
         const double length = std::abs(chi) < 1e-8 && chi != 0.0 ? 1e10 : 1.0;
         for (int symmetric = 0; symmetric < 2; ++symmetric) {
            std::ostringstream text;
            text << std::setprecision(17)
                 << "System{ Mixture{ nMonomer 2 monomers[ "
                 << 1.0/std::sqrt(length) << " " << 1.0/std::sqrt(length)
                 << " ] nPolymer 1 Polymer{ type linear nBlock 2 blocks[ 0 "
                 << length/2.0 << " 1 " << length/2.0 << " ] phi 1 } ds "
                 << length/20.0 << " } Interaction{ chi( 0 1 " << chi
                 << " ) } Domain{ mesh 8 lattice lamellar groupName P_-1 } }";
            Rp::System<1,CUT> system;
            std::istringstream input(text.str());
            system.readParam(input);
            FSArray<double,6> cell;
            cell.append(1.0);
            system.setUnitCell(cell);
            if (symmetric) {
               DArray< DArray<double> > fields;
               fields.allocate(2);
               const int nb = system.domain().basis().nBasis();
               for (int i = 0; i < 2; ++i) {
                  fields[i].allocate(nb);
                  for (int k = 0; k < nb; ++k) fields[i][k] = 0.0;
                  fields[i][0] = chi/2.0;
               }
               system.w().setBasis(fields);
            } else {
               DArray< RField<1,CUT> > fields;
               fields.allocate(2);
               for (int i = 0; i < 2; ++i) {
                  fields[i].allocate(system.domain().mesh().dimensions());
                  VecOp::eqS(fields[i], chi/2.0);
               }
               system.w().setRGrid(fields);
            }
            system.compute();
            system.scft().compute();
            const double interaction = length*chi/4.0;
            TEST_ASSERT(std::abs(length*system.scft().fInter()-interaction) < 2e-11);
            TEST_ASSERT(std::abs(length*system.scft().fHelmholtz()+1.0-interaction) < 2e-11);
            TEST_ASSERT(std::abs(length*system.scft().pressure()-1.0-interaction) < 2e-11);
         }
      }
   }

   // Orthogonal composition vectors v=(1,-1,0)/sqrt(2), u=(1,1,-2)/sqrt(6).
   // This matrix has pressure eigenvalue3*scale and composition eigenvalues
   // -scale and small*scale, independently of the eigensolver under test.

};

TEST_BEGIN(InteractionEnergyScaleTest)
TEST_ADD(InteractionEnergyScaleTest, testInteractionScaleThermo)
TEST_END(InteractionEnergyScaleTest)

#endif
