#ifndef PSCF_CORE_GROUP_FILE_TEST_H
#define PSCF_CORE_GROUP_FILE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/crystal/SpaceGroup.h>
#include <util/misc/Exception.h>

#include <cstdio>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

using namespace Util;
using namespace Pscf;

/*
* Lookup of group files by readGroup: a file in the current directory
* is read before the standard group file of the same name.
*/
class CoreGroupFileTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   /*
   * A directory in the current directory with the name of a standard
   * group (e.g. an output directory) must not be taken for a group
   * file: the standard group is read.
   */
   template <int D>
   void checkDirectory(std::string const & name)
   {
      SpaceGroup<D> reference;
      readGroup(name, reference);
      TEST_ASSERT(reference.isValid());

      rmdir(name.c_str());
      TEST_ASSERT(mkdir(name.c_str(), 0755) == 0);
      SpaceGroup<D> group;
      bool thrown = false;
      try {
         readGroup(name, group);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(rmdir(name.c_str()) == 0);
      TEST_ASSERT(!thrown);
      TEST_ASSERT(group.isValid());
      TEST_ASSERT(group.size() == reference.size());
   }

   void testDirectory2D()
   {
      printMethod(TEST_FUNC);
      checkDirectory<2>("p_6_m_m");
   }

   void testDirectory3D()
   {
      printMethod(TEST_FUNC);
      checkDirectory<3>("I_m_-3_m");
   }

   /*
   * A regular file in the current directory is still read first, also
   * if it has the name of a standard group.
   */
   void testLocalFile()
   {
      printMethod(TEST_FUNC);

      SpaceGroup<2> p2;
      readGroup("p_2", p2);
      TEST_ASSERT(p2.size() == 2);
      SpaceGroup<2> p6mm;
      readGroup("p_6_m_m", p6mm);
      TEST_ASSERT(p6mm.size() == 12);

      std::string name = "p_6_m_m";
      writeGroup(name, p2);
      SpaceGroup<2> group;
      bool thrown = false;
      try {
         readGroup(name, group);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(std::remove(name.c_str()) == 0);
      TEST_ASSERT(!thrown);
      TEST_ASSERT(group.size() == 2);
   }

};

TEST_BEGIN(CoreGroupFileTest)
TEST_ADD(CoreGroupFileTest, testDirectory2D)
TEST_ADD(CoreGroupFileTest, testDirectory3D)
TEST_ADD(CoreGroupFileTest, testLocalFile)
TEST_END(CoreGroupFileTest)

#endif
