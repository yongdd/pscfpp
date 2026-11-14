#ifndef PSCF_FIELD_ARCHIVE_ROUND_TRIP_TEST_H
#define PSCF_FIELD_ARCHIVE_ROUND_TRIP_TEST_H

#include <test/UnitTestRunner.h>
#include <pscf/field/RField.h>
#include <pscf/field/RFieldDft.h>
#include <pscf/field/CField.h>
#include <pscf/backend/HostSendArray.h>
#include <pscf/backend/HostRecvArray.h>
#include <pscf/backend/complex.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/TextFileIArchive.h>
#include <util/archives/TextFileOArchive.h>
#include <util/tests/LogFileUnitTest.h>
#include <util/misc/Exception.h>

using namespace Util;
using namespace Pscf;

class FieldArchiveRoundTripTest : public LogFileUnitTest
{
public:
   void testCpuRealRoundTrip()
   { printMethod(TEST_FUNC); checkRealBackend<CPT>(); }

   void testCpuFourierRoundTrip()
   { printMethod(TEST_FUNC); checkFourierBackend<CPT>(); }

#ifdef PSCF_CUDA
   void testGpuRealRoundTrip()
   { printMethod(TEST_FUNC); checkRealBackend<CUT>(); }

   void testGpuFourierRoundTrip()
   { printMethod(TEST_FUNC); checkFourierBackend<CUT>(); }

   void testGpuEmptyRoundTrip()
   {
      printMethod(TEST_FUNC);
      for (bool text : {false, true}) {
         RField<2,CUT> realSource, realTarget;
         roundTrip(realSource, realTarget, text);
         TEST_ASSERT(!realTarget.isAllocated());
         TEST_ASSERT(realTarget.meshDimensions() == IntVec<2>(0));
         RFieldDft<2,CUT> dftSource, dftTarget;
         roundTrip(dftSource, dftTarget, text);
         TEST_ASSERT(!dftTarget.isAllocated());
         TEST_ASSERT(dftTarget.meshDimensions() == IntVec<2>(0));
         TEST_ASSERT(dftTarget.dftDimensions() == IntVec<2>(0));
      }
   }

   void testGpuInvalidMetadataPreservesFields()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FieldArchiveInvalid.log");
      for (int failure = 0; failure < 3; ++failure) {
         RField<1,CUT> realTarget; realTarget.allocate(IntVec<1>(4));
         fillReal(realTarget);
         RFieldDft<1,CUT> dftTarget; dftTarget.allocate(IntVec<1>(4));
         fillFourier(dftTarget);
         for (bool fourier : {false, true}) {
            MemoryOArchive out; out.allocate(1024);
            int capacity = fourier ? 3 : 4;
            if (failure == 0) capacity = -1;
            if (failure == 1) ++capacity;
            out & capacity;
            if (failure == 2) {
               cudaReal value = -9.0;
               for (int i = 0; i < capacity; ++i) {
                  out & value;
                  if (fourier) out & value;
               }
               IntVec<1> invalidDimensions(0); out & invalidDimensions;
            }
            MemoryIArchive in; in = out;
            bool threw = false;
            try {
               if (fourier) in & dftTarget; else in & realTarget;
            } catch (Exception&) { threw = true; }
            in.release();
            TEST_ASSERT(threw);
         }
         checkRealValues(realTarget);
         checkFourierValues(dftTarget);
         TEST_ASSERT(realTarget.meshDimensions() == IntVec<1>(4));
         TEST_ASSERT(dftTarget.meshDimensions() == IntVec<1>(4));
         TEST_ASSERT(dftTarget.dftDimensions() == IntVec<1>(3));
      }
   }

   void testGpuComplexRoundTrip()
   {
      printMethod(TEST_FUNC);
      for (bool text : {false, true}) {
         for (bool reused : {false, true}) {
            checkComplex<1>(reused, text);
            checkComplex<2>(reused, text);
            checkComplex<3>(reused, text);
         }
      }
   }

   void testGpuComplexEmptyAndInvalid()
   {
      printMethod(TEST_FUNC);
      for (bool text : {false, true}) {
         CField<2,CUT> source, target;
         roundTrip(source, target, text);
         TEST_ASSERT(!target.isAllocated());
         TEST_ASSERT(target.meshDimensions() == IntVec<2>(0));
      }

      // Negative capacity, capacity of another size, zero dimension and
      // dimensions that do not match the capacity: the load must throw
      // before it changes the values or dimensions of the target.
      openLogFile("out/FieldArchiveComplexInvalid.log");
      for (int failure = 0; failure < 4; ++failure) {
         for (bool allocated : {false, true}) {
            if (failure == 1 && !allocated) continue;
            CField<1,CUT> target;
            if (allocated) {
               target.allocate(IntVec<1>(4));
               fillComplex(target);
            }
            MemoryOArchive out; out.allocate(1024);
            int capacity = 4;
            if (failure == 0) capacity = -1;
            if (failure == 1) capacity = 5;
            out & capacity;
            if (failure >= 2) {
               cudaReal value = -9.0;
               for (int i = 0; i < 2*capacity; ++i) out & value;
               IntVec<1> invalidDimensions(failure == 2 ? 0 : 3);
               out & invalidDimensions;
            }
            MemoryIArchive in; in = out;
            bool threw = false;
            try {
               in & target;
            } catch (Exception&) { threw = true; }
            in.release();
            TEST_ASSERT(threw);
            TEST_ASSERT(target.isAllocated() == allocated);
            if (allocated) {
               checkComplexValues(target);
               TEST_ASSERT(target.meshDimensions() == IntVec<1>(4));
            }
         }
      }
   }
#endif

private:
   template<class T> void checkRealBackend()
   {
      for (bool text : {false, true}) {
         for (bool reused : {false, true}) {
            checkReal<1,T>(reused, text);
            checkReal<2,T>(reused, text);
            checkReal<3,T>(reused, text);
         }
      }
   }

   template<class T> void checkFourierBackend()
   {
      for (bool text : {false, true}) {
         for (bool reused : {false, true}) {
            checkFourier<1,T>(reused, text);
            checkFourier<2,T>(reused, text);
            checkFourier<3,T>(reused, text);
         }
      }
   }

   template<int D, class T> void fillReal(RField<D,T>& field)
   {
      HostSendArray<typename T::Real,T> host; host.associate(field);
      for (int i = 0; i < host.capacity(); ++i) host[i] = 2.0 + .125*i;
      field = host; host.dissociate();
   }

   template<int D, class T> void fillFourier(RFieldDft<D,T>& field)
   {
      HostSendArray<typename T::Complex,T> host; host.associate(field);
      for (int i = 0; i < host.capacity(); ++i) {
         assign(host[i], typename T::Real(2.0 + .125*i),
                         typename T::Real(-.5 - .25*i));
      }
      field = host; host.dissociate();
   }

   template<int D, class T> void checkRealValues(RField<D,T> const& field)
   {
      HostRecvArray<typename T::Real,T> host(field);
      for (int i = 0; i < host.capacity(); ++i) {
         TEST_ASSERT(host[i] == 2.0 + .125*i);
      }
   }

   template<int D, class T> void checkFourierValues(RFieldDft<D,T> const& field)
   {
      HostRecvArray<typename T::Complex,T> host(field);
      for (int i = 0; i < host.capacity(); ++i) {
         TEST_ASSERT(real(host[i]) == 2.0 + .125*i);
         TEST_ASSERT(imag(host[i]) == -.5 - .25*i);
      }
   }

#ifdef PSCF_CUDA
   template<int D> void fillComplex(CField<D,CUT>& field)
   {
      HostSendArray<cudaComplex,CUT> host; host.associate(field);
      for (int i = 0; i < host.capacity(); ++i) {
         assign(host[i], cudaReal(1.5 - .25*i), cudaReal(3.0 + .125*i));
      }
      field = host; host.dissociate();
   }

   template<int D> void checkComplexValues(CField<D,CUT> const& field)
   {
      HostRecvArray<cudaComplex,CUT> host(field);
      for (int i = 0; i < host.capacity(); ++i) {
         TEST_ASSERT(real(host[i]) == 1.5 - .25*i);
         TEST_ASSERT(imag(host[i]) == 3.0 + .125*i);
      }
   }

   template<int D> void checkComplex(bool reused, bool text)
   {
      // Same number of values, but different dimensions for D > 1
      IntVec<D> dimensions(3), oldDimensions(3);
      dimensions[D-1] = 6; oldDimensions[0] = 6;
      CField<D,CUT> source, target; source.allocate(dimensions);
      fillComplex(source);
      if (reused) {
         target.allocate(oldDimensions);
         HostSendArray<cudaComplex,CUT> host; host.associate(target);
         for (int i = 0; i < host.capacity(); ++i) {
            assign(host[i], cudaReal(-99.0), cudaReal(99.0));
         }
         target = host; host.dissociate();
      }
      roundTrip(source, target, text);
      TEST_ASSERT(target.meshDimensions() == dimensions);
      checkComplexValues(target); checkComplexValues(source);
   }
#endif

   template<class Field> void roundTrip(Field& source, Field& target, bool text)
   {
      int marker = 31415, restoredMarker = 0;
      if (text) {
         TextFileOArchive out;
         openOutputFile("out/FieldArchiveRoundTrip.txt", out.file());
         out & source & marker; out.file().close();
         TextFileIArchive in;
         openInputFile("out/FieldArchiveRoundTrip.txt", in.file());
         in & target & restoredMarker;
      } else {
         MemoryOArchive out; out.allocate(16384);
         out & source & marker;
         MemoryIArchive in; in = out;
         in & target & restoredMarker; in.release();
      }
      TEST_ASSERT(restoredMarker == marker);
      TEST_ASSERT(target.capacity() == source.capacity());
   }

   template<int D, class T> void checkReal(bool reused, bool text)
   {
      IntVec<D> dimensions(3), oldDimensions(3);
      dimensions[D-1] = 6; oldDimensions[0] = 6;
      RField<D,T> source, target; source.allocate(dimensions);
      fillReal(source);
      if (reused) {
         target.allocate(oldDimensions);
         HostSendArray<typename T::Real,T> host; host.associate(target);
         for (int i = 0; i < host.capacity(); ++i) host[i] = -99.0;
         target = host; host.dissociate();
      }
      roundTrip(source, target, text);
      TEST_ASSERT(target.meshDimensions() == dimensions);
      checkRealValues(target); checkRealValues(source);
   }

   template<int D, class T> void checkFourier(bool reused, bool text)
   {
      IntVec<D> dimensions(3), oldDimensions(3), expectedDft(3);
      dimensions[D-1] = 7; expectedDft[D-1] = 4;
      oldDimensions[D-1] = 6;
      if (D > 1) {
         // Same number of complex values, but a different Fourier mesh.
         oldDimensions[0] = 4; oldDimensions[D-1] = 5;
      }
      RFieldDft<D,T> source, target; source.allocate(dimensions);
      fillFourier(source);
      if (reused) {
         target.allocate(oldDimensions);
         HostSendArray<typename T::Complex,T> host; host.associate(target);
         for (int i = 0; i < host.capacity(); ++i) {
            assign(host[i], typename T::Real(-99.0), typename T::Real(99.0));
         }
         target = host; host.dissociate();
      }
      roundTrip(source, target, text);
      TEST_ASSERT(target.meshDimensions() == dimensions);
      TEST_ASSERT(target.dftDimensions() == expectedDft);
      checkFourierValues(target); checkFourierValues(source);
   }
};

TEST_BEGIN(FieldArchiveRoundTripTest)
TEST_ADD(FieldArchiveRoundTripTest, testCpuRealRoundTrip)
TEST_ADD(FieldArchiveRoundTripTest, testCpuFourierRoundTrip)
#ifdef PSCF_CUDA
TEST_ADD(FieldArchiveRoundTripTest, testGpuRealRoundTrip)
TEST_ADD(FieldArchiveRoundTripTest, testGpuFourierRoundTrip)
TEST_ADD(FieldArchiveRoundTripTest, testGpuEmptyRoundTrip)
TEST_ADD(FieldArchiveRoundTripTest, testGpuInvalidMetadataPreservesFields)
TEST_ADD(FieldArchiveRoundTripTest, testGpuComplexRoundTrip)
TEST_ADD(FieldArchiveRoundTripTest, testGpuComplexEmptyAndInvalid)
#endif
TEST_END(FieldArchiveRoundTripTest)

#endif
