#ifndef PSCF_EDGE_ARCHIVE_STATE_TEST_H
#define PSCF_EDGE_ARCHIVE_STATE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/chem/Edge.h>
#include <pscf/chem/PolymerModel.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/TextFileIArchive.h>
#include <util/archives/TextFileOArchive.h>
#include <util/tests/LogFileUnitTest.h>
#include <util/misc/Exception.h>
#include <sstream>

using namespace Util;
using namespace Pscf;

class EdgeArchiveStateTest : public LogFileUnitTest
{
public:
   void tearDown()
   {
      LogFileUnitTest::tearDown();
      PolymerModel::setModel(PolymerModel::Thread);
   }

   void testCompleteDescriptorRoundTrip()
   {
      printMethod(TEST_FUNC);
      for (bool bead : {false, true}) {
         for (bool linear : {false, true}) {
            for (bool initialized : {false, true}) {
               check(bead, linear, initialized, false);
               check(bead, linear, initialized, true);
            }
         }
      }
   }

   void testLegacyRecordIsRejected()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/EdgeArchiveLegacy.log");
      PolymerModel::setModel(PolymerModel::Bead);
      Edge target;
      target.setId(9); target.setNBead(2);
      target.setPolymerType(PolymerType::Linear);

      // Exact old layout, followed by another record. Missing descriptor
      // members cannot be reconstructed from these old bytes.
      int id = 3, monomer = 2, marker = 31415;
      Pair<int> vertices; vertices[0] = 3; vertices[1] = 4;
      double length = -1.0;
      MemoryOArchive out; out.allocate(1024);
      out & id & monomer & vertices & length & marker;
      MemoryIArchive in; in = out;
      bool threw = false;
      try { in & target; }
      catch (Exception&) { threw = true; }
      in.release();
      TEST_ASSERT(threw);
      TEST_ASSERT(target.id() == 9);
      TEST_ASSERT(target.nBead() == 2);
      TEST_ASSERT(target.polymerType() == PolymerType::Linear);
   }

private:
   void check(bool bead, bool linear, bool initialized, bool text)
   {
      PolymerModel::setModel(bead ? PolymerModel::Bead : PolymerModel::Thread);
      Edge source, target;
      source.setId(3); source.setMonomerId(2); source.setVertexIds(3, 4);
      source.setPolymerType(linear ? PolymerType::Linear : PolymerType::Branched);
      if (bead) source.setNBead(7); else source.setLength(1.25);
      if (initialized) {
         target.setId(9); target.setMonomerId(1); target.setVertexIds(8, 9);
         target.setPolymerType(linear ? PolymerType::Branched : PolymerType::Linear);
         if (bead) target.setNBead(2); else target.setLength(0.5);
      }
      int marker = 31415, restoredMarker = 0;
      if (text) {
         TextFileOArchive out;
         openOutputFile("out/EdgeArchiveState.txt", out.file());
         out & source & marker;
         out.file().close();
         TextFileIArchive in;
         openInputFile("out/EdgeArchiveState.txt", in.file());
         in & target & restoredMarker;
      } else {
         MemoryOArchive out; out.allocate(1024);
         out & source & marker;
         MemoryIArchive in; in = out;
         in & target & restoredMarker;
         in.release();
      }
      TEST_ASSERT(target.id() == source.id());
      TEST_ASSERT(target.monomerId() == source.monomerId());
      TEST_ASSERT(target.vertexId(0) == source.vertexId(0));
      TEST_ASSERT(target.vertexId(1) == source.vertexId(1));
      TEST_ASSERT(target.polymerType() == source.polymerType());
      if (bead) {
         TEST_ASSERT(target.nBead() == 7);
      } else {
         TEST_ASSERT(target.length() == 1.25);
      }
      TEST_ASSERT(restoredMarker == marker);
      std::ostringstream before, after;
      before << source; after << target;
      TEST_ASSERT(before.str() == after.str());
   }
};

TEST_BEGIN(EdgeArchiveStateTest)
TEST_ADD(EdgeArchiveStateTest, testCompleteDescriptorRoundTrip)
TEST_ADD(EdgeArchiveStateTest, testLegacyRecordIsRejected)
TEST_END(EdgeArchiveStateTest)

#endif
