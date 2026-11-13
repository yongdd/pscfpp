#ifndef PSCF_FIELD_READ_END_OF_FILE_TEST_H
#define PSCF_FIELD_READ_END_OF_FILE_TEST_H

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

class FieldReadEndOfFileTest : public LogFileUnitTest
{
public:
   void setUp() { setVerbose(0); }

   void testRealGridEndOfFile()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FieldReadEndOfFileIoReal.log");
      IntVec<1> dims;
      dims[0] = 2;
      DArray<double> single;
      single.allocate(2);
      std::istringstream one("1 2");
      readRGridData(one, single, dims);
      TEST_ASSERT(single[0] == 1 && single[1] == 2);
      TEST_ASSERT(one.eof() && !one.fail());
      DArray< DArray<double> > fields;
      fields.allocate(2);
      fields[0].allocate(2);
      fields[1].allocate(2);
      std::istringstream many("1 2 3 4");
      readRGridData(many, fields, 2, dims);
      TEST_ASSERT(fields[0][1] == 3 && fields[1][1] == 4);
      for (int variant = 0; variant < 2; ++variant) {
         std::istringstream bad(variant ? "1 2 3 invalid" : "1 2 3");
         bool threw = false;
         try { readRGridData(bad, fields, 2, dims); }
         catch (Exception&) { threw = true; }
         TEST_ASSERT(threw);
      }
   }

   void testFourierGridEndOfFile()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FieldReadEndOfFileIoFourier.log");
      IntVec<1> dims;
      dims[0] = 2;
      RFieldDft<1,CPT> single;
      single.allocate(dims);
      std::istringstream one("0 1 2\n1 3 4");
      readKGridData(one, single, dims);
      TEST_ASSERT(real(single[1]) == 3 && imag(single[1]) == 4);
      TEST_ASSERT(one.eof() && !one.fail());
      DArray< RFieldDft<1,CPT> > fields;
      fields.allocate(2);
      fields[0].allocate(dims);
      fields[1].allocate(dims);
      std::istringstream many("0 1 2 3 4\n1 5 6 7 8");
      readKGridData(many, fields, 2, dims);
      TEST_ASSERT(real(fields[1][1]) == 7 && imag(fields[1][1]) == 8);
      for (int variant = 0; variant < 2; ++variant) {
         std::istringstream bad(variant ? "0 1 2\n1 3 invalid" : "0 1 2\n1 3");
         bool threw = false;
         try { readKGridData(bad, single, dims); }
         catch (Exception&) { threw = true; }
         TEST_ASSERT(threw);
      }
   }

   void testComplexGridEndOfFile()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FieldReadEndOfFileIoComplex.log");
      IntVec<1> dims;
      dims[0] = 2;
      DArray< std::complex<double> > single;
      single.allocate(2);
      std::istringstream one("1 2 3 4");
      readCFieldData(one, single, dims);
      TEST_ASSERT(real(single[1]) == 3 && imag(single[1]) == 4);
      TEST_ASSERT(one.eof() && !one.fail());
      DArray< DArray< std::complex<double> > > fields;
      fields.allocate(2);
      fields[0].allocate(2);
      fields[1].allocate(2);
      std::istringstream many("1 2 3 4 5 6 7 8");
      readCFieldsData(many, fields, dims);
      TEST_ASSERT(real(fields[1][1]) == 7 && imag(fields[1][1]) == 8);
      for (int variant = 0; variant < 2; ++variant) {
         std::istringstream bad(variant ? "1 2 3 invalid" : "1 2 3");
         bool threw = false;
         try { readCFieldData(bad, single, dims); }
         catch (Exception&) { threw = true; }
         TEST_ASSERT(threw);
      }
   }

};

TEST_BEGIN(FieldReadEndOfFileTest)
TEST_ADD(FieldReadEndOfFileTest, testRealGridEndOfFile)
TEST_ADD(FieldReadEndOfFileTest, testFourierGridEndOfFile)
TEST_ADD(FieldReadEndOfFileTest, testComplexGridEndOfFile)
TEST_END(FieldReadEndOfFileTest)

#endif
