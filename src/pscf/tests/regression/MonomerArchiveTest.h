#ifndef PSCF_MONOMER_ARCHIVE_TEST_H
#define PSCF_MONOMER_ARCHIVE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/chem/Monomer.h>
#include <util/archives/MemoryCounter.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/TextFileIArchive.h>
#include <util/archives/TextFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>
#include <util/archives/BinaryFileOArchive.h>
#include <fstream>

class MonomerArchiveTest : public UnitTest
{
public:
   void testMemoryOutput()
   {
      printMethod(TEST_FUNC);
      Pscf::Monomer source;
      source.setId(3);
      source.setKuhn(1.25);
      int marker = 314159;
      Util::MemoryOArchive out;
      out.allocate(Util::memorySize(source) + sizeof(marker));
      out & source & marker;
      TEST_ASSERT(out.cursor() == out.begin() + 2*sizeof(int) + sizeof(double));
      // Inspect the saved record independently of Monomer's reader.
      Util::MemoryIArchive in;
      in = out;
      int id = -1, restoredMarker = -1;
      double kuhn = -1.0;
      in & id & kuhn & restoredMarker;
      in.release();
      TEST_ASSERT(id == 3);
      TEST_ASSERT(kuhn == 1.25);
      TEST_ASSERT(restoredMarker == marker);
      TEST_ASSERT(source.id() == 3 && source.kuhn() == 1.25);
   }

   void testMemoryInput()
   {
      printMethod(TEST_FUNC);
      // A primitive fixture isolates cursor advancement from the writer.
      int id = 3, marker = 314159, restoredMarker = -1;
      double kuhn = 1.25;
      Util::MemoryOArchive out;
      out.allocate(128);
      out & id & kuhn & marker;
      Util::MemoryIArchive in;
      in = out;
      Pscf::Monomer restored;
      restored.setId(9);
      restored.setKuhn(2.5);
      in & restored;
      TEST_ASSERT(in.cursor() == in.begin() + sizeof(int) + sizeof(double));
      in & restoredMarker;
      in.release();
      TEST_ASSERT(restored.id() == id && restored.kuhn() == kuhn);
      TEST_ASSERT(restoredMarker == marker);
   }

   void testTextRoundTrip()
   {
      printMethod(TEST_FUNC);
      checkFile<Util::TextFileOArchive, Util::TextFileIArchive>(
         "out/MonomerArchiveText", false);
      checkFile<Util::TextFileOArchive, Util::TextFileIArchive>(
         "out/MonomerArchiveBorrowedText", true);
   }

   void testBinaryRoundTrip()
   {
      printMethod(TEST_FUNC);
      checkFile<Util::BinaryFileOArchive, Util::BinaryFileIArchive>(
         "out/MonomerArchiveBinary", false);
      checkFile<Util::BinaryFileOArchive, Util::BinaryFileIArchive>(
         "out/MonomerArchiveBorrowedBinary", true);
   }

private:
   template <class Output, class Input>
   void checkFile(const char* name, bool borrowed)
   {
      std::string path = filePrefix() + name;
      Pscf::Monomer first, second, a, b;
      first.setId(3); first.setKuhn(1.25);
      second.setId(7); second.setKuhn(2.5);
      int marker = 314159, restoredMarker = -1;
      if (borrowed) {
         {
            std::ofstream file(path.c_str(), std::ios::binary);
            Output out(file);
            out & first & second & marker;
         }
         {
            std::ifstream file(path.c_str(), std::ios::binary);
            Input in(file);
            in & a & b & restoredMarker;
         }
      } else {
         {
            Output out(path);
            out & first & second & marker;
         }
         {
            Input in(path);
            in & a & b & restoredMarker;
         }
      }
      TEST_ASSERT(a.id() == 3 && a.kuhn() == 1.25);
      TEST_ASSERT(b.id() == 7 && b.kuhn() == 2.5);
      TEST_ASSERT(restoredMarker == marker);
      TEST_ASSERT(first.id() == 3 && first.kuhn() == 1.25);
      TEST_ASSERT(second.id() == 7 && second.kuhn() == 2.5);
   }
};

TEST_BEGIN(MonomerArchiveTest)
TEST_ADD(MonomerArchiveTest, testMemoryOutput)
TEST_ADD(MonomerArchiveTest, testMemoryInput)
TEST_ADD(MonomerArchiveTest, testTextRoundTrip)
TEST_ADD(MonomerArchiveTest, testBinaryRoundTrip)
TEST_END(MonomerArchiveTest)

#endif
