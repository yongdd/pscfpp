#ifndef PSCF_CORE_EDGE_TEST_H
#define PSCF_CORE_EDGE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/chem/Edge.h>
#include <pscf/chem/PolymerModel.h>
#include <pscf/chem/PolymerType.h>

#include <cmath>
#include <sstream>

using namespace Pscf;
using namespace Util;

class CoreEdgeTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {  PolymerModel::setModel(PolymerModel::Thread); }

   /*
   * A block length written by operator << and read back by operator >>
   * must be the same double. Lengths were written with 13 significant
   * digits, so a written parameter file did not reproduce the input.
   */
   void testWriteReadLength()
   {
      printMethod(TEST_FUNC);

      PolymerModel::setModel(PolymerModel::Thread);

      double lengths[4];
      lengths[0] = 1.0/3.0;
      lengths[1] = 0.1 + 0.2;
      lengths[2] = std::nextafter(0.56, 1.0);
      lengths[3] = 2.0*std::acos(-1.0);

      Edge v, u;
      v.setPolymerType(PolymerType::Branched);
      u.setPolymerType(PolymerType::Branched);
      v.setMonomerId(1);
      v.setVertexIds(2, 5);
      for (int i = 0; i < 4; ++i) {
         v.setLength(lengths[i]);
         std::stringstream stream;
         stream << v;
         stream >> u;
         TEST_ASSERT(u.monomerId() == 1);
         TEST_ASSERT(u.vertexId(0) == 2);
         TEST_ASSERT(u.vertexId(1) == 5);
         TEST_ASSERT(u.length() == lengths[i]);
      }
   }

};

TEST_BEGIN(CoreEdgeTest)
TEST_ADD(CoreEdgeTest, testWriteReadLength)
TEST_END(CoreEdgeTest)

#endif
