#ifndef PSCF_EDGE_OUTPUT_MODE_TEST_H
#define PSCF_EDGE_OUTPUT_MODE_TEST_H
#include <pscf/chem/Edge.h>
#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <iomanip>
#include <sstream>

class EdgeOutputModeTest : public UnitTest
{
public:
   void setUp() { model_ = Pscf::PolymerModel::model(); }
   void tearDown() { Pscf::PolymerModel::setModel(model_); }
   void testThreadRoundTrip()
   {
      printMethod(TEST_FUNC);
      check(false);
   }
   void testBeadRoundTrip()
   {
      printMethod(TEST_FUNC);
      check(true);
   }
private:
   Pscf::PolymerModel::Type model_;
   void check(bool bead)
   {
      using namespace Pscf;
      PolymerModel::setModel(bead ? PolymerModel::Bead : PolymerModel::Thread);
      for (auto type : {PolymerType::Linear, PolymerType::Branched}) {
         for (int mode = 0; mode < 4; ++mode) {
            Edge source, target;
            source.setId(3); source.setMonomerId(1);
            source.setPolymerType(type); source.setVertexIds(3,4);
            if (bead) source.setNBead(7); else source.setLength(1.2345678901234567);
            std::ostringstream out;
            if (mode == 1) out << std::fixed;
            if (mode == 2) out << std::scientific;
            if (mode == 3) out << std::hexfloat;
            out << source << " " << 314159;
            if (!bead) TEST_ASSERT(out.str().find("0x") == std::string::npos);
            target.setId(3); target.setPolymerType(type);
            if (type == PolymerType::Linear) target.setVertexIds(3,4);
            std::istringstream in(out.str());
            int marker = 0;
            in >> target >> marker;
            TEST_ASSERT(!in.fail() && marker == 314159);
            TEST_ASSERT(target.id() == source.id());
            TEST_ASSERT(target.monomerId() == source.monomerId());
            TEST_ASSERT(target.polymerType() == source.polymerType());
            TEST_ASSERT(target.vertexId(0) == 3 && target.vertexId(1) == 4);
            if (bead) {
               TEST_ASSERT(target.nBead() == source.nBead());
            } else {
               TEST_ASSERT(target.length() == source.length());
            }
         }
      }
   }
};
TEST_BEGIN(EdgeOutputModeTest)
TEST_ADD(EdgeOutputModeTest, testThreadRoundTrip)
TEST_ADD(EdgeOutputModeTest, testBeadRoundTrip)
TEST_END(EdgeOutputModeTest)
#endif
