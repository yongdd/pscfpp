#ifndef PSCF_FFTW_COMPLEX_ARRAY_COPY_TEST_H
#define PSCF_FFTW_COMPLEX_ARRAY_COPY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/backend/cpp/BackendArray.h>
#include <util/misc/Memory.h>

class FftwComplexArrayCopyTest : public UnitTest
{
   using Complex = Pscf::CPT::Complex;
   using Array = Pscf::BackendArray<Complex, Pscf::CPT>;

public:
   void testCopyConstructor()
   {
      printMethod(TEST_FUNC);
      const size_t memory = Util::Memory::total();
      {
         Array source(3);
         fill(source);
         Array copy(source);
         TEST_ASSERT(copy.isOwner() && !copy.isAssociated());
         TEST_ASSERT(copy.cArray() != source.cArray());
         source[1][0] = 99.0; source[1][1] = -99.0;
         check(copy);
      }
      TEST_ASSERT(Util::Memory::total() == memory);
   }

   void testAssignment()
   {
      printMethod(TEST_FUNC);
      for (bool allocated : {false, true}) {
         Array source(3), copy;
         fill(source);
         if (allocated) copy.allocate(3);
         copy = source;
         TEST_ASSERT(copy.isOwner() && !copy.isAssociated());
         TEST_ASSERT(copy.cArray() != source.cArray());
         source[1][0] = 99.0; source[1][1] = -99.0;
         check(copy);
         copy = copy;
         check(copy);
      }
   }

   void testArrayViewAssignment()
   {
      printMethod(TEST_FUNC);
      for (bool allocated : {false, true}) {
         Array source(3), copy;
         fill(source);
         if (allocated) copy.allocate(3);
         Util::Array<Complex> const& view = source;
         copy = view;
         TEST_ASSERT(copy.isOwner() && !copy.isAssociated());
         TEST_ASSERT(copy.cArray() != source.cArray());
         source[1][0] = 99.0; source[1][1] = -99.0;
         check(copy);
         Util::Array<Complex> const& self = copy;
         copy = self;
         check(copy);
      }
   }

private:
   void fill(Array& source)
   {
      const double values[3][2] = {{1.25,-.5}, {-2.5,4.0}, {3.75,-8.0}};
      for (int i = 0; i < 3; ++i) {
         source[i][0] = values[i][0]; source[i][1] = values[i][1];
      }
   }
   void check(Array const& copy)
   {
      TEST_ASSERT(copy.capacity() == 3);
      TEST_ASSERT(copy[0][0] == 1.25 && copy[0][1] == -.5);
      TEST_ASSERT(copy[1][0] == -2.5 && copy[1][1] == 4.0);
      TEST_ASSERT(copy[2][0] == 3.75 && copy[2][1] == -8.0);
   }
};
TEST_BEGIN(FftwComplexArrayCopyTest)
TEST_ADD(FftwComplexArrayCopyTest, testCopyConstructor)
TEST_ADD(FftwComplexArrayCopyTest, testAssignment)
TEST_ADD(FftwComplexArrayCopyTest, testArrayViewAssignment)
TEST_END(FftwComplexArrayCopyTest)
#endif
