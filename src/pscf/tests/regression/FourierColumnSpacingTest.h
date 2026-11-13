#ifndef PSCF_FOURIER_COLUMN_SPACING_TEST_H
#define PSCF_FOURIER_COLUMN_SPACING_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/fieldIo/rFieldIo.tpp>
#include <pscf/fieldIo/cFieldIo.tpp>
#include <pscf/field/cpp/RFieldDft.h>
#include <util/tests/LogFileUnitTest.h>
#include <util/misc/Exception.h>
#include <sstream>
#include <complex>
#include <vector>

using namespace Util;
using namespace Pscf;

class FourierColumnSpacingTest : public LogFileUnitTest
{
public:
   void setUp() { setVerbose(0); }

   void checkColumns(std::string const& text, int count)
   {
      std::istringstream rows(text);
      std::string line, word;
      int nrow = 0;
      while (std::getline(rows, line)) {
         std::istringstream row(line);
         int n = 0;
         while (row >> word) ++n;
         TEST_ASSERT(n == count);
         ++nrow;
      }
      TEST_ASSERT(nrow == 2);
   }

   void testFourierColumnSeparators()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FourierColumnSpacingIoSeparators.log");
      IntVec<1> dims;
      dims[0] = 2;
      DArray< RFieldDft<1,CPT> > fields;
      fields.allocate(2);
      for (int j = 0; j < 2; ++j) {
         fields[j].allocate(dims);
         for (int i = 0; i < 2; ++i) {
            assign(fields[j][i], 1e-100, -5e-101);
         }
      }
      std::ostringstream single, many;
      writeKGridData(single, fields[0], dims);
      writeKGridData(many, fields, 2, dims);
      checkColumns(single.str(), 3);
      checkColumns(many.str(), 5);
   }
};

TEST_BEGIN(FourierColumnSpacingTest)
TEST_ADD(FourierColumnSpacingTest, testFourierColumnSeparators)
TEST_END(FourierColumnSpacingTest)

#endif
