#ifndef RPC_SCFT_THERMO_REF_TEST_H
#define RPC_SCFT_THERMO_REF_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/system/SystemConstRef.h>
#include <rp/scft/ScftThermo.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <rp/field/Mask.h>

#include <pscf/interaction/Interaction.h>

#include <util/misc/FileMaster.h>

using namespace Util;
using namespace Pscf;

/*
* Access to the System components that a ScftThermo refers to. These
* accessors, inherited from SystemConstRef, are protected in ScftThermo.
*/
template <int D>
class ScftThermoRefProbe : public Rp::ScftThermo<D,CPT>
{

public:

   typedef Rp::SystemConstRef<D,CPT> RefT;

   /*
   * Does every reference held by thermo point to the component of
   * system that the System accessor returns?
   */
   static bool refersTo(Rp::ScftThermo<D,CPT> const & thermo,
                        Rp::System<D,CPT> const & system)
   {
      Rp::System<D,CPT> const & (RefT::*getSystem)() const
                                = &ScftThermoRefProbe::system;
      Rp::Mixture<D,CPT> const & (RefT::*getMixture)() const
                                = &ScftThermoRefProbe::mixture;
      Interaction const & (RefT::*getInteraction)() const
                                = &ScftThermoRefProbe::interaction;
      Rp::Domain<D,CPT> const & (RefT::*getDomain)() const
                                = &ScftThermoRefProbe::domain;
      Rp::CFields<D,CPT> const & (RefT::*getC)() const
                                = &ScftThermoRefProbe::c;
      Rp::WFields<D,CPT> const & (RefT::*getW)() const
                                = &ScftThermoRefProbe::w;
      Rp::WFields<D,CPT> const & (RefT::*getH)() const
                                = &ScftThermoRefProbe::h;
      Rp::Mask<D,CPT> const & (RefT::*getMask)() const
                                = &ScftThermoRefProbe::mask;
      FileMaster const & (RefT::*getFileMaster)() const
                                = &ScftThermoRefProbe::fileMaster;

      bool ok = true;
      ok = ok && (&(thermo.*getSystem)() == &system);
      ok = ok && (&(thermo.*getMixture)() == &system.mixture());
      ok = ok && (&(thermo.*getInteraction)() == &system.interaction());
      ok = ok && (&(thermo.*getDomain)() == &system.domain());
      ok = ok && (&(thermo.*getC)() == &system.c());
      ok = ok && (&(thermo.*getW)() == &system.w());
      ok = ok && (&(thermo.*getH)() == &system.h());
      ok = ok && (&(thermo.*getMask)() == &system.mask());
      ok = ok && (&(thermo.*getFileMaster)() == &system.fileMaster());
      return ok;
   }

};

/*
* Regression tests of the associations made in the System constructor.
*/
class ScftThermoRefTest : public UnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * The ScftThermo object of a new System refers to the components of
   * that System, including its FileMaster (the ScftThermo was created
   * before the FileMaster and bound a reference through a null
   * pointer; with UTIL_DEBUG the System constructor aborted).
   */
   void testConstructorReferences()
   {
      printMethod(TEST_FUNC);

      Rp::System<1,CPT> system1;
      TEST_ASSERT(ScftThermoRefProbe<1>::refersTo(system1.scft(),
                                                  system1));

      Rp::System<2,CPT> system2;
      TEST_ASSERT(ScftThermoRefProbe<2>::refersTo(system2.scft(),
                                                  system2));

      Rp::System<3,CPT> system3;
      TEST_ASSERT(ScftThermoRefProbe<3>::refersTo(system3.scft(),
                                                  system3));
   }

};

TEST_BEGIN(ScftThermoRefTest)
TEST_ADD(ScftThermoRefTest, testConstructorReferences)
TEST_END(ScftThermoRefTest)

#endif
