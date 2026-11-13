#ifndef PSCF_BEAD_COUNT_VALIDATION_TEST_H
#define PSCF_BEAD_COUNT_VALIDATION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/chem/Edge.h>
#include <pscf/chem/PolymerModel.h>
#include <pscf/chem/PolymerType.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/field/cpp/FFT.h>
#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldDft.h>
#include <pscf/field/cpp/WaveList.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/containers/FSArray.h>
#include <util/misc/Exception.h>
#include <cmath>
#include <sstream>

using namespace Util;
using namespace Pscf;

class BeadCountValidationTest : public UnitTest
{
public:
   void setUp() {}
   void tearDown() { PolymerModel::setModel(PolymerModel::Thread); }

   void testPositiveBeadCounts()
   {
      printMethod(TEST_FUNC);
      PolymerModel::setModel(PolymerModel::Bead);
      for (int count = -1; count <= 1; ++count) {
         Edge edge;
         edge.setPolymerType(PolymerType::Linear);
         bool threw = false;
         try { edge.setNBead(count); } catch (Exception&) { threw = true; }
         TEST_ASSERT(threw == (count <= 0));
         std::stringstream input;
         input << "0 " << count;
         threw = false;
         try { input >> edge; } catch (Exception&) { threw = true; }
         TEST_ASSERT(threw == (count <= 0));
         if (count > 0) TEST_ASSERT(edge.nBead() == count);
      }
   }

   // Rounded right angles are physical perturbations, not exact ties.
   // Every inverse pair must still have the same minimum squared norm.

};

TEST_BEGIN(BeadCountValidationTest)
TEST_ADD(BeadCountValidationTest, testPositiveBeadCounts)
TEST_END(BeadCountValidationTest)

#endif
