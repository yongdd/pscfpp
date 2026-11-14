#ifndef PSCF_DECIMAL_OUTPUT_MODE_TEST_H
#define PSCF_DECIMAL_OUTPUT_MODE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/math/RealVec.h>
#include <pscf/chem/Monomer.h>
#include <pscf/floryHuggins/FhClump.h>
#include <iomanip>
#include <sstream>
#include <string>

class PscfDecimalOutputModeTest : public UnitTest
{
public:
   void testRealVectors()
   {
      printMethod(TEST_FUNC);
      checkVector<1, double>();
      checkVector<2, double>();
      checkVector<3, double>();
      checkVector<3, float>();
   }

   void testMonomers()
   {
      printMethod(TEST_FUNC);
      for (int mode = 0; mode < 4; ++mode) {
         Pscf::Monomer first, second, a, b;
         first.setId(3); first.setKuhn(1.25);
         second.setId(4); second.setKuhn(2.5);
         a.setId(7); b.setId(8);
         std::ostringstream out;
         setMode(out, mode);
         out << first << ' ' << second << ' ' << 314159;
         TEST_ASSERT(out.str().find("0x") == std::string::npos);
         std::istringstream in(out.str());
         int marker = 0;
         in >> a >> b >> marker;
         TEST_ASSERT(!in.fail() && marker == 314159);
         TEST_ASSERT(a.kuhn() == first.kuhn());
         TEST_ASSERT(b.kuhn() == second.kuhn());
         // Text records contain only the length, preserving assigned ids.
         TEST_ASSERT(a.id() == 7 && b.id() == 8);
      }
   }

   void testClumps()
   {
      printMethod(TEST_FUNC);
      for (int mode = 0; mode < 4; ++mode) {
         Pscf::FhClump first, second, a, b;
         first.setMonomerId(1); first.setSize(1.25);
         second.setMonomerId(2); second.setSize(2.5);
         std::ostringstream out;
         setMode(out, mode);
         out << first << ' ' << second << ' ' << 314159;
         TEST_ASSERT(out.str().find("0x") == std::string::npos);
         std::istringstream in(out.str());
         int marker = 0;
         in >> a >> b >> marker;
         TEST_ASSERT(!in.fail() && marker == 314159);
         TEST_ASSERT(a.monomerId() == first.monomerId());
         TEST_ASSERT(b.monomerId() == second.monomerId());
         TEST_ASSERT(a.size() == first.size());
         TEST_ASSERT(b.size() == second.size());
      }
   }

private:
   void setMode(std::ostream& out, int mode)
   {
      if (mode == 1) out << std::scientific;
      if (mode == 2) out << std::fixed;
      if (mode == 3) out << std::hexfloat;
   }

   template <int D, typename T>
   void checkVector()
   {
      for (int mode = 0; mode < 4; ++mode) {
         Pscf::RealVec<D, T> source(T(0)), target(T(-99));
         for (int i = 0; i < D; ++i) {
            source[i] = T((i % 2 ? -1.0 : 1.0)*1.25*(i+1));
         }
         std::ostringstream out;
         setMode(out, mode);
         out << source << ' ' << 314159;
         TEST_ASSERT(out.str().find("0x") == std::string::npos);
         std::istringstream in(out.str());
         int marker = 0;
         in >> target >> marker;
         TEST_ASSERT(!in.fail() && marker == 314159);
         for (int i = 0; i < D; ++i) {
            TEST_ASSERT(target[i] == source[i]);
         }
      }
   }
};

TEST_BEGIN(PscfDecimalOutputModeTest)
TEST_ADD(PscfDecimalOutputModeTest, testRealVectors)
TEST_ADD(PscfDecimalOutputModeTest, testMonomers)
TEST_ADD(PscfDecimalOutputModeTest, testClumps)
TEST_END(PscfDecimalOutputModeTest)

#endif
