#ifndef RPC_WALL_FIELD_SCALE_TEST_H
#define RPC_WALL_FIELD_SCALE_TEST_H

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
class WallFieldScaleTestBase : public LogFileUnitTest
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

   std::vector<double> values(RField<1,CPT> const& field)
   {
      return std::vector<double>(field.cArray(),
                                 field.cArray() + field.capacity());
   }

   // Manual's sigmoid at fixed fractional coordinate, L=1, T=.4, t=.2.
   double mask(int i)
   {
      double arg = 4.0*(-0.3 + std::abs(i/64.0-0.5))/0.2;
      return 0.5*(1.0-std::tanh(arg));
   }
};

class WallFieldScaleTest : public WallFieldScaleTestBase
{
public:
   void testWallReferenceScaling()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/WallFieldScaleEnvironmentWallScaling.log");
      for (double sign : {1.0, -1.0, 0.0}) {
         for (double n : {1.0, 1e7, 1e8, 1e9}) {
            Rp::System<1,CPT> system;
            initialize(system, n, sign);
            setFields(system, 0.0);
            system.compute();
            system.scft().compute();
            TEST_ASSERT(system.h().hasData() == (sign != 0.0));
            if (sign == 0.0) {
               TEST_ASSERT(system.scft().fExt() == 0.0);
               continue;
            }
            double expectedEnergy = 0.0;
            double phi = system.mask().phiTot();
            for (int m = 0; m < 2; ++m) {
               std::vector<double> h = values(system.h().rgrid(m));
               std::vector<double> c = values(system.c().rgrid(m));
               for (int i = 0; i < 64; ++i) {
                  double chi = i < 32 ? m+1.0 : (m == 0 ? 5.0 : 3.0);
                  double expected = sign*chi*(1.0-mask(i));
                  TEST_ASSERT(std::abs(n*h[i]-expected) < 1e-11);
                  TEST_ASSERT(std::abs(c[i]-phi/2.0) < 1e-12);
                  expectedEnergy += expected/128.0;
               }
            }
            TEST_ASSERT(std::abs(n*system.scft().fExt()-expectedEnergy)
                        < 1e-10);
         }
      }
   }

   void testSmallWallSymmetry()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/WallFieldScaleEnvironmentWallSymmetry.log");
      // Reflection is compatible with equal walls, including tiny values.
      Rp::System<1,CPT> symmetric;
      initialize(symmetric, 1e8, 1.0, 20.0, true, true);
      setFields(symmetric, 0.0);
      symmetric.compute();
      TEST_ASSERT(symmetric.h().hasData());
      TEST_ASSERT(symmetric.h().isSymmetric());

      // Tiny asymmetric fields must not be labeled reflection-symmetric.
      bool threw = false;
      try {
         Rp::System<1,CPT> asymmetric;
         initialize(asymmetric, 1e8, 1.0, 20.0, false, true);
         setFields(asymmetric, 0.0);
         asymmetric.compute();
      } catch (Exception&) {
         threw = true;
      }
      TEST_ASSERT(threw);
   }

};

TEST_BEGIN(WallFieldScaleTest)
TEST_ADD(WallFieldScaleTest, testWallReferenceScaling)
TEST_ADD(WallFieldScaleTest, testSmallWallSymmetry)
TEST_END(WallFieldScaleTest)

#endif
