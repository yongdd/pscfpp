#ifndef RPC_MASK_STRESS_SCALE_TEST_H
#define RPC_MASK_STRESS_SCALE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <rp/system/System.h>
#include <rp/scft/ScftThermo.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <rp/field/Mask.h>
#include <rp/environment/Environment.h>
#include <pscf/backend/HostSendArray.h>
#include <util/misc/Exception.h>
#include <util/tests/LogFileUnitTest.h>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace Util;
using namespace Pscf;

// Shared small-film fixtures for numerical and public-API regressions.
class MaskStressScaleTestBase : public LogFileUnitTest
{
public:
   void setUp() { setVerbose(0); }

   void initialize(Rp::System<1,CPT>& system, double n, double wallScale,
                   double chiN = 20.0, bool symmetric = false,
                   bool useGroup = false)
   {
      std::ostringstream text;
      text << std::setprecision(17)
           << "System{ Mixture{ nMonomer 2 monomers[ "
           << 1.0/std::sqrt(n) << ' ' << 1.0/std::sqrt(n)
           << " ] nPolymer 1 Polymer{ type linear nBlock 2 blocks[ 0 "
           << n/2.0 << " 1 " << n/2.0 << " ] phi 1 } vMonomer "
           << 1.0/n << " ds " << n/100.0 << " } Interaction{ chi( 0 1 "
           << chiN/n << " ) } Domain{ mesh 64 lattice lamellar ";
      if (useGroup) text << "groupName P_-1 ";
      text << "} FilmEnvironment{ normalVecId 0 interfaceThickness .2 "
           << "excludedThickness .4 fBulk 0 chiBottom[ "
           << wallScale/n << ' ' << 2.0*wallScale/n << " ] chiTop[ "
           << (symmetric ? 1.0 : 5.0)*wallScale/n << ' '
           << (symmetric ? 2.0 : 3.0)*wallScale/n << " ] } "
           << "AmIteratorGrid{ epsilon 1e-9 maxItr 100 maxHist 10 "
           << "isFlexible 1 } }";
      std::istringstream input(text.str());
      system.readParam(input);
      FSArray<double,6> cell;
      cell.append(1.0);
      system.setUnitCell(cell);
   }

   void setFields(Rp::System<1,CPT>& system, double amplitude)
   {
      DArray<RField<1,CPT> > fields;
      fields.allocate(2);
      for (int m = 0; m < 2; ++m) {
         fields[m].allocate(system.domain().mesh().dimensions());
         HostSendArray<CPT::Real,CPT> host;
         host.associate(fields[m]);
         for (int i = 0; i < 64; ++i) {
            host[i] = (m == 0 ? 1.0 : -0.6)*amplitude
                    * std::cos(2.0*std::acos(-1.0)*i/64.0);
         }
         fields[m] = host;
         host.dissociate();
      }
      system.w().setRGrid(fields);
   }

   // Manual's sigmoid at fixed fractional coordinate, L=1, T=.4, t=.2.
   double mask(int i)
   {
      double arg = 4.0*(-0.3 + std::abs(i/64.0-0.5))/0.2;
      return 0.5*(1.0-std::tanh(arg));
   }
};

class MaskStressScaleTest : public MaskStressScaleTestBase
{
public:

   void testMaskStressReferenceScaling()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/MaskStressScaleEnvironmentMaskScaling.log");
      double phi = 0.0;
      double dPhi = 0.0;
      for (int i = 0; i < 64; ++i) {
         double u = i/64.0;
         double arg = 4.0*(-0.3+std::abs(u-0.5))/0.2;
         phi += mask(i)/64.0;
         // Differentiate the manual's tanh formula at fixed u, T and t.
         dPhi -= 2.0*(std::abs(u-0.5)-0.5)
                /(0.2*std::cosh(arg)*std::cosh(arg)*64.0);
      }
      for (double chiN : {20.0, -20.0, 0.0}) {
         // With uniform w=0, c_A=c_B=phi/2 and mixture stress is zero.
         double exactF = std::log(phi)-1.0+chiN*phi/4.0;
         double exactP = 1.0-chiN*phi/4.0;
         double exactStress = 0.6*(exactP+chiN*phi/2.0)*dPhi/phi+exactF;
         for (double n : {1.0, 1e7, 1e8, 1e9}) {
            Rp::System<1,CPT> system;
            initialize(system, n, 0.0, chiN);
            setFields(system, 0.0);
            system.compute(true);
            TEST_ASSERT(!system.h().hasData());
            TEST_ASSERT(std::abs(system.mask().phiTot()-phi) < 1e-13);
            TEST_ASSERT(std::abs(n*system.scft().fHelmholtz()-exactF) < 1e-10);
            TEST_ASSERT(std::abs(n*system.scft().pressure()-exactP) < 1e-10);
            TEST_ASSERT(std::abs(n*system.mixture().stress(0)) < 1e-12);
            TEST_ASSERT(std::abs(n*system.environment().stress(0)-exactStress)
                        < 1e-10);
         }
      }
   }
};

TEST_BEGIN(MaskStressScaleTest)
TEST_ADD(MaskStressScaleTest, testMaskStressReferenceScaling)
TEST_END(MaskStressScaleTest)

#endif
