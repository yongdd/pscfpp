#ifndef PSCF_CORE_MIX_AND_MATCH_ENV_TEST_H
#define PSCF_CORE_MIX_AND_MATCH_ENV_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/environment/MixAndMatchEnvTmpl.h>
#include <pscf/environment/EnvironmentBase.h>
#include <pscf/environment/FieldGeneratorBase.h>
#include <util/misc/Exception.h>

#include <sstream>
#include <string>

using namespace Util;
using namespace Pscf;

/*
* Field generator of a given type that reads one parameter.
*/
class CoreTestFieldGenerator : public FieldGeneratorBase
{

public:

   CoreTestFieldGenerator(Type type, std::string const & label)
    : label_(label),
      value_(0.0)
   {
      type_ = type;
      setClassName("CoreTestFieldGenerator");
   }

   void readParameters(std::istream& in)
   {  read(in, label_.c_str(), value_); }

   bool needsUpdate() const
   {  return false; }

   void checkCompatibility()
   {}

   double value() const
   {  return value_; }

protected:

   void compute()
   {}

private:

   std::string label_;
   double value_;

};

/*
* MixAndMatchEnvTmpl with one or two generators of given types.
*/
class CoreTestEnvironment
 : public MixAndMatchEnvTmpl<EnvironmentBase, FieldGeneratorBase>
{

public:

   typedef FieldGeneratorBase::Type Type;

   CoreTestEnvironment(Type type1, Type type2, bool hasSecond)
    : type1_(type1),
      type2_(type2),
      hasSecond_(hasSecond)
   {  setClassName("CoreTestEnvironment"); }

protected:

   void createGenerators()
   {
      fieldGenPtr1_ = new CoreTestFieldGenerator(type1_, "value1");
      if (hasSecond_) {
         fieldGenPtr2_ = new CoreTestFieldGenerator(type2_, "value2");
      }
   }

private:

   Type type1_;
   Type type2_;
   bool hasSecond_;

};

/*
* MixAndMatchEnvTmpl::readParameters must accept a first generator of
* type Mask or External, with or without a second generator of the
* other type. A first generator of type External was rejected ("must
* have type Mask or External"), as the type Mask was tested twice.
*/
class CoreMixAndMatchEnvTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   /*
   * Read the parameters of an environment. Returns true if an
   * Exception was thrown.
   */
   bool readThrows(CoreTestEnvironment& env, std::string const & text)
   {
      std::istringstream in(text);
      try {
         env.readParameters(in);
      } catch (Exception& e) {
         return true;
      }
      return false;
   }

   void testExternalThenMask()
   {
      printMethod(TEST_FUNC);
      CoreTestEnvironment env(FieldGeneratorBase::External,
                              FieldGeneratorBase::Mask, true);
      TEST_ASSERT(!readThrows(env, "value1  1.5\nvalue2  2.5\n"));
      TEST_ASSERT(env.generatesExternalFields());
      TEST_ASSERT(env.generatesMask());
      TEST_ASSERT(env.hasFieldGenerator2());
      TEST_ASSERT(env.fieldGenerator1().type()
                  == FieldGeneratorBase::External);
      CoreTestFieldGenerator const & gen1 =
         dynamic_cast<CoreTestFieldGenerator const &>(env.fieldGenerator1());
      CoreTestFieldGenerator const & gen2 =
         dynamic_cast<CoreTestFieldGenerator const &>(env.fieldGenerator2());
      TEST_ASSERT(gen1.value() == 1.5);
      TEST_ASSERT(gen2.value() == 2.5);
   }

   void testExternalOnly()
   {
      printMethod(TEST_FUNC);
      CoreTestEnvironment env(FieldGeneratorBase::External,
                              FieldGeneratorBase::None, false);
      TEST_ASSERT(!readThrows(env, "value1  1.5\n"));
      TEST_ASSERT(env.generatesExternalFields());
      TEST_ASSERT(!env.generatesMask());
      TEST_ASSERT(!env.hasFieldGenerator2());
   }

   /*
   * Controls: the order used by FilmEnvironment, and invalid pairs.
   */
   void testMaskFirst()
   {
      printMethod(TEST_FUNC);
      CoreTestEnvironment env1(FieldGeneratorBase::Mask,
                               FieldGeneratorBase::External, true);
      TEST_ASSERT(!readThrows(env1, "value1  1.5\nvalue2  2.5\n"));
      TEST_ASSERT(env1.generatesExternalFields());
      TEST_ASSERT(env1.generatesMask());

      CoreTestEnvironment env2(FieldGeneratorBase::Mask,
                               FieldGeneratorBase::None, false);
      TEST_ASSERT(!readThrows(env2, "value1  1.5\n"));
      TEST_ASSERT(!env2.generatesExternalFields());
      TEST_ASSERT(env2.generatesMask());

      CoreTestEnvironment env3(FieldGeneratorBase::Mask,
                               FieldGeneratorBase::Mask, true);
      TEST_ASSERT(readThrows(env3, "value1  1.5\nvalue2  2.5\n"));

      CoreTestEnvironment env4(FieldGeneratorBase::External,
                               FieldGeneratorBase::External, true);
      TEST_ASSERT(readThrows(env4, "value1  1.5\nvalue2  2.5\n"));

      CoreTestEnvironment env5(FieldGeneratorBase::None,
                               FieldGeneratorBase::None, false);
      TEST_ASSERT(readThrows(env5, "value1  1.5\n"));
   }

};

TEST_BEGIN(CoreMixAndMatchEnvTest)
TEST_ADD(CoreMixAndMatchEnvTest, testExternalThenMask)
TEST_ADD(CoreMixAndMatchEnvTest, testExternalOnly)
TEST_ADD(CoreMixAndMatchEnvTest, testMaskFirst)
TEST_END(CoreMixAndMatchEnvTest)

#endif
