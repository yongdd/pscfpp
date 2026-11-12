#ifndef PSCF_CORE_GROUP_INPUT_TEST_H
#define PSCF_CORE_GROUP_INPUT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/crystal/SpaceGroup.h>
#include <pscf/crystal/Basis.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/interaction/Interaction.h>
#include <util/containers/FSArray.h>
#include <util/param/BracketPolicy.h>
#include <util/misc/Exception.h>

#include <fstream>
#include <string>

using namespace Util;
using namespace Pscf;

/*
* Parameter and group file input that is read by the util library.
*/
class CoreGroupInputTest : public UnitTest 
{

public:

   void setUp()
   {  BracketPolicy::set(BracketPolicy::Optional); }

   void tearDown()
   {}

   /*
   * A group file equal to c_m, with the translations written as 0/2
   * and 2/4 (in/Core/c_m_unreduced), must be accepted and give the
   * same basis as c_m. Fractions read unreduced made the group check
   * fail ('Product of two element is not in group').
   */
   void testUnreducedTranslations()
   {
      printMethod(TEST_FUNC);

      SpaceGroup<2> group;
      std::ifstream in;
      openInputFile("in/Core/c_m_unreduced", in);
      in >> group;
      in.close();
      TEST_ASSERT(group.isValid());
      TEST_ASSERT(group.size() == 4);

      SpaceGroup<2> reference;
      readGroup("c_m", reference);
      TEST_ASSERT(reference.size() == 4);

      UnitCell<2> unitCell;
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      parameters.append(3.0);
      unitCell.set(UnitCell<2>::Rectangular, parameters);
      IntVec<2> d;
      d[0] = 24;
      d[1] = 24;
      Mesh<2> mesh(d);

      Basis<2> basis;
      basis.makeBasis(mesh, unitCell, group);
      Basis<2> basisRef;
      basisRef.makeBasis(mesh, unitCell, reference);
      TEST_ASSERT(basis.isValid());
      TEST_ASSERT(basis.nStar() == basisRef.nStar());
      TEST_ASSERT(basis.nBasis() == basisRef.nBasis());
   }

   /*
   * A chi parameter with 1-based indices, 'chi( 2 1 20.0 )' for two
   * monomers, must be rejected: it was accepted, left chi at zero and
   * wrote past the arrays.
   */
   void testChiIndexOutOfRange()
   {
      printMethod(TEST_FUNC);

      Interaction v;
      v.setNMonomer(2);
      std::ifstream in;
      openInputFile("in/Core/InteractionOneBased", in);
      bool thrown = false;
      std::string message;
      try {
         v.readParam(in);
      } catch (Exception& e) {
         thrown = true;
         message = e.message();
      }
      in.close();
      TEST_ASSERT(thrown);
      TEST_ASSERT(message.find("out of range") != std::string::npos);
   }

};

TEST_BEGIN(CoreGroupInputTest)
TEST_ADD(CoreGroupInputTest, testUnreducedTranslations)
TEST_ADD(CoreGroupInputTest, testChiIndexOutOfRange)
TEST_END(CoreGroupInputTest)

#endif
