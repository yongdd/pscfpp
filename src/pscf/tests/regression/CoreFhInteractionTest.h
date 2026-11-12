#ifndef PSCF_CORE_FH_INTERACTION_TEST_H
#define PSCF_CORE_FH_INTERACTION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/floryHuggins/FhInteraction.h>
#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>
#include <util/param/BracketPolicy.h>

#include <fstream>

using namespace Pscf;
using namespace Util;

class CoreFhInteractionTest : public UnitTest 
{

public:

   void setUp()
   {  BracketPolicy::set(BracketPolicy::Optional); }

   void tearDown()
   {}

   /*
   * A singular chi matrix is valid input, e.g. one zero chi parameter
   * in a 3-monomer system (here, monomers 0 and 1 do not interact with
   * each other and interact identically with monomer 2). It used to
   * abort the process through the GSL error handler when it was read.
   * Functions that need only chi work; those that need its inverse
   * throw an Exception.
   */
   void testSingularChi3()
   {
      printMethod(TEST_FUNC);

      FhInteraction v;
      v.setNMonomer(3);
      std::ifstream in;
      openInputFile("in/Core/FhInteractionSingular", in);
      v.readParam(in);
      in.close();

      TEST_ASSERT(eq(v.chi(0,1), 0.0));
      TEST_ASSERT(eq(v.chi(0,2), 10.0));
      TEST_ASSERT(eq(v.chi(1,2), 10.0));
      TEST_ASSERT(!v.hasChiInverse());

      // Functions of chi alone
      DArray<double> c;
      DArray<double> w;
      c.allocate(3);
      w.allocate(3);
      c[0] = 0.2;
      c[1] = 0.3;
      c[2] = 0.5;
      v.computeW(c, w);
      TEST_ASSERT(eq(w[0], 5.0));
      TEST_ASSERT(eq(w[1], 5.0));
      TEST_ASSERT(eq(w[2], 5.0));
      TEST_ASSERT(eq(v.fHelmholtz(c), 2.5));

      // Functions of the inverse of chi
      bool thrown = false;
      try {
         v.chiInverse(0, 0);
      } catch (Exception&) {
         thrown = true;
      }
      TEST_ASSERT(thrown);
      thrown = false;
      double xi;
      try {
         v.computeC(w, c, xi);
      } catch (Exception&) {
         thrown = true;
      }
      TEST_ASSERT(thrown);

      // An invertible chi defines the inverse again
      DMatrix<double> chi;
      chi.allocate(3, 3);
      for (int i = 0; i < 3; ++i) {
         for (int j = 0; j < 3; ++j) {
            chi(i, j) = v.chi(i, j);
         }
      }
      chi(0, 1) = 5.0;
      chi(1, 0) = 5.0;
      v.setChi(chi);
      TEST_ASSERT(v.hasChiInverse());
      double sum;
      for (int i = 0; i < 3; ++i) {
         for (int j = 0; j < 3; ++j) {
            sum = 0.0;
            for (int k = 0; k < 3; ++k) {
               sum += v.chi(i, k)*v.chiInverse(k, j);
            }
            TEST_ASSERT(eq(sum, (i == j) ? 1.0 : 0.0));
         }
      }

      // Making chi singular again must not leave a stale inverse
      v.setChi(0, 1, 0.0);
      TEST_ASSERT(!v.hasChiInverse());
      thrown = false;
      try {
         v.computeC(w, c, xi);
      } catch (Exception&) {
         thrown = true;
      }
      TEST_ASSERT(thrown);
   }

   /*
   * For 2 monomers, chi = 0 (an ideal mixture) is valid input, but it
   * has no inverse: chiInverse() used to return values computed by a
   * division by a zero determinant. When chi becomes singular after
   * an invertible chi, computeC must not use the stale inverse.
   */
   void testSingularChi2()
   {
      printMethod(TEST_FUNC);

      FhInteraction v;
      v.setNMonomer(2);
      v.setChi(0, 1, 2.0);
      TEST_ASSERT(v.hasChiInverse());
      TEST_ASSERT(eq(v.chiInverse(0, 0), 0.0));
      TEST_ASSERT(eq(v.chiInverse(0, 1), 0.5));
      TEST_ASSERT(eq(v.chiInverse(1, 1), 0.0));

      v.setChi(0, 1, 0.0);
      TEST_ASSERT(!v.hasChiInverse());
      bool thrown = false;
      try {
         v.chiInverse(0, 1);
      } catch (Exception&) {
         thrown = true;
      }
      TEST_ASSERT(thrown);

      DArray<double> c;
      DArray<double> w;
      c.allocate(2);
      w.allocate(2);
      w[0] = 0.4;
      w[1] = 0.6;
      double xi;
      thrown = false;
      try {
         v.computeC(w, c, xi);
      } catch (Exception&) {
         thrown = true;
      }
      TEST_ASSERT(thrown);
   }

};

TEST_BEGIN(CoreFhInteractionTest)
TEST_ADD(CoreFhInteractionTest, testSingularChi3)
TEST_ADD(CoreFhInteractionTest, testSingularChi2)
TEST_END(CoreFhInteractionTest)

#endif
