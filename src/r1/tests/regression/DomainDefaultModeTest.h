#ifndef R1_DOMAIN_DEFAULT_MODE_TEST_H
#define R1_DOMAIN_DEFAULT_MODE_TEST_H
#include <r1/field/Domain.h>
#include <util/math/Constants.h>
#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <sstream>

class DomainDefaultModeTest : public UnitTest
{
   bool echo_;
public:
   void setUp()
   {
      echo_ = Util::ParamComponent::echo();
      Util::ParamComponent::setEcho(false);
      Util::Label::clear();
   }
   void tearDown()
   {
      Util::Label::clear();
      Util::ParamComponent::setEcho(echo_);
   }
   void check(char const* body, Pscf::R1::GeometryMode mode,
              double lower, bool shell, double volume)
   {
      Pscf::R1::Domain domain;
      std::istringstream in(std::string("Domain{ ")+body+" } next 314159\n");
      domain.readParam(in);
      TEST_ASSERT(domain.mode() == mode && domain.nx() == 5);
      TEST_ASSERT(domain.xMin() == lower && domain.xMax() == 2.0);
      TEST_ASSERT(domain.isShell() == shell);
      TEST_ASSERT(eq(domain.dx(), (2.0-lower)/4.0));
      TEST_ASSERT(eq(domain.volume(), volume));
      std::string label; int marker = 0;
      in >> label >> marker;
      TEST_ASSERT(label == "next" && marker == 314159);
      std::ostringstream out;
      domain.writeParam(out);
      Pscf::R1::Domain restored;
      std::istringstream saved(out.str());
      restored.readParam(saved);
      TEST_ASSERT(restored.mode() == mode && restored.isShell() == shell);
      TEST_ASSERT(restored.nx() == 5 && restored.xMin() == lower);
      TEST_ASSERT(eq(restored.volume(), volume));
   }
   void testDefaultPlanar()
   {
      printMethod(TEST_FUNC);
      check("xMax 2 nx 5", Pscf::R1::Planar, 0.0, false, 2.0);
      check("xMin -0.5 xMax 2 nx 5", Pscf::R1::Planar, -0.5, false, 2.5);
   }
   void testExplicitModes()
   {
      printMethod(TEST_FUNC);
      check("mode planar xMin 0 xMax 2 nx 5", Pscf::R1::Planar, 0, false, 2);
      check("mode cylindrical xMax 2 nx 5", Pscf::R1::Cylindrical, 0, false,
            4*Util::Constants::Pi);
      check("mode spherical xMax 2 nx 5", Pscf::R1::Spherical, 0, false,
            32*Util::Constants::Pi/3);
      check("mode cylindrical xMin 0.5 xMax 2 nx 5", Pscf::R1::Cylindrical,
            0.5, true, 3.75*Util::Constants::Pi);
   }
};
TEST_BEGIN(DomainDefaultModeTest)
TEST_ADD(DomainDefaultModeTest, testDefaultPlanar)
TEST_ADD(DomainDefaultModeTest, testExplicitModes)
TEST_END(DomainDefaultModeTest)
#endif
