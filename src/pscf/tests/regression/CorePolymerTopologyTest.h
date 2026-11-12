#ifndef PSCF_CORE_POLYMER_TOPOLOGY_TEST_H
#define PSCF_CORE_POLYMER_TOPOLOGY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include "../solvers/PolymerStub.h"
#include <util/param/BracketPolicy.h>

#include <fstream>
#include <string>

using namespace Pscf;

class CorePolymerTopologyTest : public UnitTest 
{

public:

   void setUp()
   {  Util::BracketPolicy::set(Util::BracketPolicy::Optional); }

   void tearDown()
   {}

   /*
   * Block graphs that are not trees must be rejected when read. For a
   * graph with a cycle, makePlan() never terminated, so on the code
   * without the check these tests hang (run them under timeout).
   */
   void testReadParamRing()
   {
      printMethod(TEST_FUNC);

      // Two-block ring: both blocks join vertices 0 and 1
      PolymerStub p;
      TEST_ASSERT(readParamThrows(p, "in/Core/PolymerRing"));
   }

   void testReadParamTriangle()
   {
      printMethod(TEST_FUNC);

      // Blocks 0, 1, 2 form the loop 0-1-2-0 (vertex 3 is unused)
      PolymerStub p;
      TEST_ASSERT(readParamThrows(p, "in/Core/PolymerTriangle"));
   }

   /*
   * A vertex id >= nVertex = nBlock + 1 must be rejected before it is
   * used as an index. The unchecked code then called Vertex::addEdge
   * on an element beyond the end of the vertex array, which may also
   * throw, by chance, so the message is checked too.
   */
   void testReadParamVertexRange()
   {
      printMethod(TEST_FUNC);

      PolymerStub p;
      std::string message;
      TEST_ASSERT(readParamThrows(p, "in/Core/PolymerVertexRange", 
                                  &message));
      TEST_ASSERT(message.find("out of range") != std::string::npos);
   }

private:

   /*
   * Read a polymer from a file. Return true if this throws an Exception,
   * and copy its message to *message if message is not null.
   */
   bool readParamThrows(PolymerStub& p, std::string fileName,
                        std::string* message = 0)
   {
      std::ifstream in;
      openInputFile(fileName, in);
      bool thrown = false;
      try {
         p.readParam(in);
      } catch (Util::Exception& e) {
         thrown = true;
         if (message) *message = e.message();
      }
      in.close();
      return thrown;
   }

};

TEST_BEGIN(CorePolymerTopologyTest)
TEST_ADD(CorePolymerTopologyTest, testReadParamRing)
TEST_ADD(CorePolymerTopologyTest, testReadParamTriangle)
TEST_ADD(CorePolymerTopologyTest, testReadParamVertexRange)
TEST_END(CorePolymerTopologyTest)

#endif
