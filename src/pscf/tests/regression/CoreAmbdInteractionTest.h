#ifndef PSCF_CORE_AMBD_INTERACTION_TEST_H
#define PSCF_CORE_AMBD_INTERACTION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/interaction/Interaction.h>
#include <pscf/iterator/AmbdInteraction.h>
#include <util/param/BracketPolicy.h>

#include <fstream>
#include <string>

using namespace Pscf;

class CoreAmbdInteractionTest : public UnitTest 
{

public:

   void setUp()
   { Util::BracketPolicy::set(Util::BracketPolicy::Optional); }

   void tearDown()
   {}

   /*
   * The AM residual is defined with the inverse of chi. A singular chi
   * (one zero chi parameter in a 3-monomer system, or chi = 0 for 2
   * monomers) must give an Exception: it aborted the process through
   * the GSL error handler (3 monomers) or gave values divided by a zero
   * determinant (2 monomers).
   */
   void testSingularChi()
   {
      printMethod(TEST_FUNC);
      TEST_ASSERT(updateThrows("in/Core/InteractionSingular", 3));
      TEST_ASSERT(updateThrows("in/Core/InteractionZero2", 2));
   }

   /*
   * An invertible chi whose inverse has a zero sum of elements (zero
   * diagonal, chi_01 = chi_12 = 10, chi_02 = 40): the AM residual
   * divides by this sum, so update() must throw instead of dividing
   * by zero.
   */
   void testZeroSumInverse()
   {
      printMethod(TEST_FUNC);
      TEST_ASSERT(updateThrows("in/Core/InteractionZeroSum", 3));
   }

private:

   /*
   * Read an Interaction and call AmbdInteraction::update with it.
   * Return true if update throws an Exception.
   */
   bool updateThrows(std::string fileName, int nMonomer)
   {
      Interaction v;
      v.setNMonomer(nMonomer);
      std::ifstream in;
      openInputFile(fileName, in);
      v.readParam(in);
      in.close();

      AmbdInteraction u;
      u.setNMonomer(nMonomer);
      bool thrown = false;
      try {
         u.update(v);
      } catch (Util::Exception&) {
         thrown = true;
      }
      return thrown;
   }

};

TEST_BEGIN(CoreAmbdInteractionTest)
TEST_ADD(CoreAmbdInteractionTest, testSingularChi)
TEST_ADD(CoreAmbdInteractionTest, testZeroSumInverse)
TEST_END(CoreAmbdInteractionTest)

#endif
