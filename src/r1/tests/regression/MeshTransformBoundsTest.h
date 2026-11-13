#ifndef R1_MESH_TRANSFORM_BOUNDS_TEST_H
#define R1_MESH_TRANSFORM_BOUNDS_TEST_H

#include "R1TestBase.h"
#include <test/UnitTestRunner.h>
#include <r1/field/FieldIo.h>
#include <r1/field/Domain.h>
#include <limits>

class MeshTransformBoundsTest : public R1TestBase
{
public:

   void testMeshTransformBounds()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/MeshTransformBoundsMeshBounds.log");
      System sys;
      readSystem(sys, "in/R1/planar_nr2.prm");
      FieldIo io;
      io.associate(sys.domain(), sys.fileMaster());
      runCommands(sys, "READ_W in/R1/planar.w\nFINISH\n");
      for (int nx = -1; nx <= 1; ++nx) {
         std::ostringstream out; out << "unchanged";
         bool thrown = false;
         try { io.remesh(sys.wFields(), nx, out); }
         catch (Exception&) { thrown = true; }
         TEST_ASSERT(thrown);
         TEST_ASSERT(out.str() == "unchanged");
      }
      for (int m : {-1, std::numeric_limits<int>::max()}) {
         std::ostringstream out; out << "unchanged";
         bool thrown = false;
         try { io.extend(sys.wFields(), m, out); }
         catch (Exception&) { thrown = true; }
         TEST_ASSERT(thrown);
         TEST_ASSERT(out.str() == "unchanged");
      }
      for (int nx : {2, 3}) {
         std::ostringstream out;
         io.remesh(sys.wFields(), nx, out);
         checkRows(out.str(), nx);
      }
      for (int m : {0, 2}) {
         std::ostringstream out;
         io.extend(sys.wFields(), m, out);
         checkRows(out.str(), sys.domain().nx() + m);
      }
      // Filename overloads must reject before truncating existing output.
      runCommands(sys, "WRITE_W out/meshtransformbounds.keep\nFINISH\n");
      const std::string old = fileContents("out/meshtransformbounds.keep");
      for (const std::string command : {"REMESH_W 1 out/meshtransformbounds.keep\n",
                                       "EXTEND_W -1 out/meshtransformbounds.keep\n"}) {
         bool thrown = false;
         try { runCommands(sys, command); }
         catch (Exception&) { thrown = true; }
         TEST_ASSERT(thrown);
         TEST_ASSERT(fileContents("out/meshtransformbounds.keep") == old);
      }
   }

private:
   void checkRows(std::string const& text, int expected)
   {
      std::istringstream in(text);
      std::string label;
      int nx, nm, index;
      in >> label >> nx;
      TEST_ASSERT(label == "nx" && nx == expected);
      in >> label >> nm;
      TEST_ASSERT(label == "nm" && nm == 2);
      for (int i = 0; i < expected; ++i) {
         double a, b;
         in >> index >> a >> b;
         TEST_ASSERT(!in.fail());
         TEST_ASSERT(index == i);
      }
      in >> label;
      TEST_ASSERT(in.eof());
   }
};

TEST_BEGIN(MeshTransformBoundsTest)
TEST_ADD(MeshTransformBoundsTest, testMeshTransformBounds)
TEST_END(MeshTransformBoundsTest)

#endif
