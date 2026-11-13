#ifndef RPC_FTS_SPECTRUM_SCALE_TEST_H
#define RPC_FTS_SPECTRUM_SCALE_TEST_H

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

class FtsSpectrumScaleTest : public LogFileUnitTest
{
public:
   void setUp() { setVerbose(0); }

   void readTernary(Rp::System<1,CPT>& system)
   {
      std::istringstream input(
         "System{ Mixture{ nMonomer 3 monomers[ 1 1 1 ] nPolymer 1 "
         "Polymer{ type linear nBlock 3 blocks[ 0 1 1 1 2 1 ] phi 1 } "
         "ds .05 } Interaction{ chi( 0 1 1 0 2 1 1 2 1 ) } "
         "Domain{ mesh 8 lattice lamellar } }");
      system.readParam(input);
   }

   void testChiSpectralScaling()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsSpectrumScaleChiSpectralScale.log");
      Rp::System<1,CPT> system;
      readTernary(system);
      Rp::BdSimulator<1,CPT> simulator(system);
      std::istringstream simInput("BdSimulator{ seed 1234 }");
      simulator.readParam(simInput);
      const double scales[] = {1e-12, 1e-9, 1.0, 1e9, 1e12};
      for (double scale : scales) {
         for (int i = 0; i < 3; ++i) {
            for (int j = i+1; j < 3; ++j) system.interaction().setChi(i,j,scale);
         }
         simulator.analyzeChi();
         for (int a = 0; a < 2; ++a) {
            TEST_ASSERT(std::abs(simulator.chiEval(a)/scale+1.0) < 1e-12);
         }
         TEST_ASSERT(simulator.chiEval(2) == 0.0);
         for (int i = 0; i < 3; ++i) {
            TEST_ASSERT(std::abs(simulator.chiEvecs(2,i)-1.0) < 1e-12);
         }
      }
   }

   // Orthogonal composition vectors v=(1,-1,0)/sqrt(2), u=(1,1,-2)/sqrt(6).
   // This matrix has pressure eigenvalue3*scale and composition eigenvalues
   // -scale and small*scale, independently of the eigensolver under test.
   void setSpectrum(Rp::System<1,CPT>& system, double scale, double small)
   {
      const double v[] = {1/std::sqrt(2.0), -1/std::sqrt(2.0), 0.0};
      const double u[] = {1/std::sqrt(6.0), 1/std::sqrt(6.0), -2/std::sqrt(6.0)};
      for (int i = 0; i < 3; ++i) {
         for (int j = i; j < 3; ++j) {
            system.interaction().setChi(i,j,scale*(1.0-v[i]*v[j]+small*u[i]*u[j]));
         }
      }
   }

   void testChiNearNullComposition()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsSpectrumScaleChiNearNull.log");
      Rp::System<1,CPT> system;
      readTernary(system);
      Rp::BdSimulator<1,CPT> simulator(system);
      std::istringstream simInput("BdSimulator{ seed 1234 }");
      simulator.readParam(simInput);
      const double scales[] = {1e-12, 1.0, 1e12};
      for (double scale : scales) {
         setSpectrum(system, scale, -1e-6);
         simulator.analyzeChi();
         double a = simulator.chiEval(0)/scale;
         double b = simulator.chiEval(1)/scale;
         TEST_ASSERT(std::abs(std::min(a,b)+1.0) < 1e-12);
         TEST_ASSERT(std::abs(std::max(a,b)/1e-6+1.0) < 1e-8);
         TEST_ASSERT(simulator.chiEval(2) == 0.0);
      }
   }

   void testChiRejectsInvalidComposition()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsSpectrumScaleChiInvalid.log");
      Rp::System<1,CPT> system;
      readTernary(system);
      Rp::BdSimulator<1,CPT> simulator(system);
      std::istringstream simInput("BdSimulator{ seed 1234 }");
      simulator.readParam(simInput);
      const double scales[] = {1e-12, 1.0, 1e12};
      const double invalid[] = {0.0, 0.1};
      for (double scale : scales) {
         for (double small : invalid) {
            setSpectrum(system, scale, small);
            bool threw = false;
            try { simulator.analyzeChi(); }
            catch (Exception const &) { threw = true; }
            TEST_ASSERT(threw);
            setSpectrum(system, scale, -0.1);
            simulator.analyzeChi();
            TEST_ASSERT(simulator.chiEval(0) < 0.0);
            TEST_ASSERT(simulator.chiEval(1) < 0.0);
         }
      }
   }
};

TEST_BEGIN(FtsSpectrumScaleTest)
TEST_ADD(FtsSpectrumScaleTest, testChiSpectralScaling)
TEST_ADD(FtsSpectrumScaleTest, testChiNearNullComposition)
TEST_ADD(FtsSpectrumScaleTest, testChiRejectsInvalidComposition)
TEST_END(FtsSpectrumScaleTest)

#endif
