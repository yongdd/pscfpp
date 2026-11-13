#ifndef RPC_HOMOGENEOUS_CLUMP_SCALE_TEST_H
#define RPC_HOMOGENEOUS_CLUMP_SCALE_TEST_H

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

class HomogeneousClumpScaleTest : public LogFileUnitTest
{
public:
   void setUp() { setVerbose(0); }

   void testShortClumps()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/HomogeneousClumpScaleShortClumps.log");
      const double scales[] = {1.0, 1e-7, 1e-8, 5e-9, 1e-12};
      for (double scale : scales) {
         std::ostringstream text;
         text << std::setprecision(17)
              << "System{ Mixture{ nMonomer 3 monomers[ 1 1 1 ] nPolymer 1 "
              << "Polymer{ type linear nBlock 2 blocks[ 0 " << scale
              << " 1 " << 2*scale << " ] phi 1 } ds " << scale/10.0
              << " } Interaction{ chi( 0 1 1 ) } "
              << "Domain{ mesh 8 lattice lamellar } }";
         Rp::System<1,CPT> system;
         std::istringstream input(text.str());
         system.readParam(input);
         FhMixture reference;
         reference.initialize(system.mixture());
         TEST_ASSERT(reference.molecule(0).nClump() == 2);
         TEST_ASSERT(reference.molecule(0).clump(0).monomerId() == 0);
         TEST_ASSERT(reference.molecule(0).clump(1).monomerId() == 1);
         TEST_ASSERT(std::abs(reference.molecule(0).size()/(3*scale)-1) < 1e-14);
         DArray<double> phi;
         phi.allocate(1);
         phi[0] = 1.0;
         reference.setComposition(phi);
         TEST_ASSERT(std::abs(reference.c(0)-1.0/3.0) < 1e-14);
         TEST_ASSERT(std::abs(reference.c(1)-2.0/3.0) < 1e-14);
         TEST_ASSERT(reference.c(2) == 0.0);
      }
   }

   // Orthogonal composition vectors v=(1,-1,0)/sqrt(2), u=(1,1,-2)/sqrt(6).
   // This matrix has pressure eigenvalue3*scale and composition eigenvalues
   // -scale and small*scale, independently of the eigensolver under test.

};

TEST_BEGIN(HomogeneousClumpScaleTest)
TEST_ADD(HomogeneousClumpScaleTest, testShortClumps)
TEST_END(HomogeneousClumpScaleTest)

#endif
