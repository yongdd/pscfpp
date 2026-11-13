#ifndef RPG_DILUTE_THERMODYNAMICS_TEST_H
#define RPG_DILUTE_THERMODYNAMICS_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <rp/system/System.h>
#include <rp/scft/ScftThermo.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/ramp/LinearRamp.h>
#include <rp/fts/ramp/RampParameter.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <pscf/field/RField.h>
#include <pscf/backend/VecOp.h>
#include <pscf/interaction/Interaction.h>
#include <util/tests/LogFileUnitTest.h>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace Util;
using namespace Pscf;

// Exact homogeneous thermodynamics and simultaneous parameter updates.
class DiluteThermodynamicsTest : public LogFileUnitTest
{
public:
   void setUp() { setVerbose(0); }

   void setFields(Rp::System<1,CUT>& system, double value)
   {
      DArray< RField<1,CUT> > fields;
      fields.allocate(2);
      for (int i = 0; i < 2; ++i) {
         fields[i].allocate(system.domain().mesh().dimensions());
         VecOp::eqS(fields[i], value);
      }
      system.w().setRGrid(fields);
   }

   void checkThermo(bool solvent, bool openMinor)
   {
      const double fractions[] = {0.0, 5e-9, 1e-8, 1.01e-8, 0.1};
      for (double phi : fractions) {
         if (openMinor && phi == 0.0) continue;
         std::ostringstream text;
         text << std::setprecision(17)
              << "System{ Mixture{ nMonomer 2 monomers[ 1 1 ] "
              << "nPolymer " << (solvent ? 1 : 2)
              << " nSolvent " << (solvent ? 1 : 0)
              << " Polymer{ type linear nBlock 1 blocks[ 0 1 ] phi "
              << 1.0-phi << " } ";
         if (solvent) {
            text << "Solvent{ monomerId 1 size 1 ";
         } else {
            text << "Polymer{ type linear nBlock 1 blocks[ 1 1 ] ";
         }
         if (openMinor) text << "mu " << std::log(phi);
         else text << "phi " << phi;
         text << " } ds .05 } Interaction{ chi( 0 1 10 ) } "
              << "Domain{ mesh 8 lattice lamellar } }";
         Rp::System<1,CUT> system;
         std::istringstream input(text.str());
         system.readParam(input);
         FSArray<double,6> cell;
         cell.append(1.0);
         system.setUnitCell(cell);
         Rp::BdSimulator<1,CUT> simulator(system);
         std::istringstream simInput("BdSimulator{ seed 1234 }");
         simulator.readParam(simInput);
         simulator.analyzeChi();
         const double entropyMajor = (1.0-phi)*std::log1p(-phi);
         const double entropyMinor = phi > 0.0 ? phi*std::log(phi) : 0.0;
         const double exactF = entropyMajor + entropyMinor - 1.0
                             + 10.0*phi*(1.0-phi);
         const double exactH = 1.5 + entropyMajor
                             + (openMinor ? 0.0 : entropyMinor);
         const int nGauge = openMinor ? 1 : 2;
         for (int k = 0; k < nGauge; ++k) {
            double gauge = k ? 20.0 : 0.0;
            setFields(system, gauge);
            system.compute();
            system.scft().compute();
            TEST_ASSERT(std::abs(system.scft().fHelmholtz()-exactF) < 2e-11);
            // The off-shell thermodynamic pressure shifts with w.
            double exactP = gauge + 1.0 - 10.0*phi*(1.0-phi);
            TEST_ASSERT(std::abs(system.scft().pressure()-exactP) < 2e-11);
            simulator.computeWc();
            simulator.computeHamiltonian();
            TEST_ASSERT(std::abs(simulator.hamiltonian()-exactH) < 2e-11);
         }
      }
   }

   void testDilutePolymers()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/DiluteThermodynamicsNumericalPolymers.log");
      checkThermo(false, false);
      checkThermo(false, true);
   }

   void testDiluteSolvent()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/DiluteThermodynamicsNumericalSolvent.log");
      checkThermo(true, false);
      checkThermo(true, true);
   }

};

TEST_BEGIN(DiluteThermodynamicsTest)
TEST_ADD(DiluteThermodynamicsTest, testDilutePolymers)
TEST_ADD(DiluteThermodynamicsTest, testDiluteSolvent)
TEST_END(DiluteThermodynamicsTest)

#endif
