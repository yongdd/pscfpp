#ifndef R1_FIELD_READ_TRANSACTION_TEST_H
#define R1_FIELD_READ_TRANSACTION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <r1/field/Domain.h>
#include <r1/field/FieldIo.h>
#include <util/tests/LogFileUnitTest.h>
#include <util/misc/FileMaster.h>
#include <util/misc/Exception.h>
#include <sstream>

using namespace Util;
using namespace Pscf;

class FieldReadTransactionTest : public LogFileUnitTest
{
public:
   void setUp() { setVerbose(0); }

   void checkRead(std::string const& last, bool valid)
   {
      R1::Domain domain;
      domain.setPlanarParameters(0.0, 1.0, 3);
      FileMaster files;
      R1::FieldIo io;
      io.associate(domain, files);
      DArray< DArray<double> > fields;
      fields.allocate(2);
      for (int j = 0; j < 2; ++j) {
         fields[j].allocate(3);
         for (int i = 0; i < 3; ++i) fields[j][i] = 10*j + i;
      }
      std::istringstream in("nx 3\nnm 2\n0 1 2\n1 3 4\n2 5 " + last);
      bool threw = false;
      try { io.readFields(fields, in); }
      catch (Exception&) { threw = true; }
      TEST_ASSERT(threw != valid);
      for (int j = 0; j < 2; ++j) {
         for (int i = 0; i < 3; ++i) {
            double expected = valid ? 2*i + j + 1 : 10*j + i;
            TEST_ASSERT(fields[j][i] == expected);
         }
      }
   }

   void testCompleteFinalNumber()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FieldReadTransactionR1IoComplete.log");
      checkRead("6", true);
      checkRead("6\n", true);
   }

   void testMissingFinalNumber()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FieldReadTransactionR1IoMissing.log");
      checkRead("", false);
   }

   void testInvalidFinalNumber()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FieldReadTransactionR1IoInvalid.log");
      checkRead("invalid", false);
   }
};

TEST_BEGIN(FieldReadTransactionTest)
TEST_ADD(FieldReadTransactionTest, testCompleteFinalNumber)
TEST_ADD(FieldReadTransactionTest, testMissingFinalNumber)
TEST_ADD(FieldReadTransactionTest, testInvalidFinalNumber)
TEST_END(FieldReadTransactionTest)

#endif
