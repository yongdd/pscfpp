#ifndef RPG_FILM_DILUTE_POLYMER_TEST_H
#define RPG_FILM_DILUTE_POLYMER_TEST_H
#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <rp/system/System.h>
#include <rp/scft/ScftThermo.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/environment/FilmFieldGenMaskBase.h>
#include <pscf/backend/VecOp.h>
#include <pscf/chem/PolymerModel.h>
#include <util/tests/LogFileUnitTest.h>
#include <sstream>
#include <iomanip>
#include <cmath>
using namespace Util;
using namespace Pscf;

class FilmDilutePolymerTest : public LogFileUnitTest
{
   // Expose the thermodynamic conversion used by the film stress.
   struct Potential : Rp::FilmFieldGenMaskBase<1> {
      using Rp::FilmFieldGenMaskBase<1>::semiGrandPotential;
   };
public:
   void setUp() { setVerbose(0); }
   void tearDown() { PolymerModel::setModel(PolymerModel::Thread); }
   void testOpenPolymerPotential()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmDilutePolymer.log");
      for (int bead=0; bead<2; ++bead) {
         PolymerModel::setModel(bead ? PolymerModel::Bead : PolymerModel::Thread);
         for (double mu : {-16.0, -20.0, -1000.0}) {
            const double phi = std::exp(mu);
            std::ostringstream text;
            text << std::setprecision(17)
                 << "System{ polymerModel " << (bead ? "bead " : "thread ")
                 << "Mixture{ nMonomer 2 monomers[ 1 1 ] nPolymer 1 nSolvent 1 "
                 << "Polymer{ type linear nBlock 1 blocks[ 1 4 ] mu " << mu << " } "
                 << "Solvent{ monomerId 0 size 1 phi " << 1.0-phi << " } "
                 << (bead ? "" : "ds .05 ")
                 << "} Interaction{ chi( 0 1 2 ) } Domain{ mesh 8 lattice lamellar } }";
            Rp::System<1,CUT> system;
            std::istringstream input(text.str());
            system.readParam(input);
            FSArray<double,6> cell;
            cell.append(1.0);
            system.setUnitCell(cell);
            DArray< RField<1,CUT> > fields;
            fields.allocate(2);
            for (int i=0; i<2; ++i) {
               fields[i].allocate(system.domain().mesh().dimensions());
               VecOp::eqS(fields[i], 0.0);
            }
            system.w().setRGrid(fields);
            system.compute();
            system.scft().compute();
            // Homogeneous ideal entropy plus chi*cA*cB, with the open
            // polymer's chemical work removed. Includes underflow to zero.
            const double exact = (1.0-phi)*(std::log1p(-phi)-1.0)
                               - phi/4.0 + 2.0*phi*(1.0-phi);
            const double actual = Potential::semiGrandPotential(
               system.mixture(), system.scft().fHelmholtz());
            TEST_ASSERT(std::isfinite(actual));
            TEST_ASSERT(std::abs(actual-exact) < 1e-12);
         }
      }
   }
};
TEST_BEGIN(FilmDilutePolymerTest)
TEST_ADD(FilmDilutePolymerTest, testOpenPolymerPotential)
TEST_END(FilmDilutePolymerTest)
#endif
