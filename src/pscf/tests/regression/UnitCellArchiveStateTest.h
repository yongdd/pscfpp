#ifndef PSCF_UNIT_CELL_ARCHIVE_STATE_TEST_H
#define PSCF_UNIT_CELL_ARCHIVE_STATE_TEST_H

#include <test/UnitTestRunner.h>
#include <pscf/crystal/UnitCell.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/TextFileIArchive.h>
#include <util/archives/TextFileOArchive.h>
#include <util/tests/LogFileUnitTest.h>
#include <util/signal/Signal.h>
#include <util/math/Constants.h>
#include <util/misc/Exception.h>
#include <cmath>

using namespace Util;
using namespace Pscf;

class UnitCellArchiveStateTest : public LogFileUnitTest
{
public:
   void testGeometryRoundTrip()
   {
      printMethod(TEST_FUNC);
      for (bool text : {false, true}) {
         for (bool reused : {false, true}) {
            check<1>(UnitCell<1>::Lamellar, reused, text);
            check<2>(UnitCell<2>::Rectangular, reused, text);
            check<3>(UnitCell<3>::Orthorhombic, reused, text);
         }
      }
   }

   void testInvalidRecordPreservesCell()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/UnitCellArchiveInvalid.log");
      UnitCell<1> target;
      FSArray<double,6> lengths; lengths.append(3.0);
      target.set(UnitCell<1>::Lamellar, lengths);
      Counter count; Signal<void> signal;
      signal.addObserver(count, &Counter::tick); target.setSignal(signal);
      for (int failure = 0; failure < 2; ++failure) {
         MemoryOArchive out; out.allocate(1024);
         UnitCell<1>::LatticeSystem lattice = UnitCell<1>::Lamellar;
         int nParameter = failure == 0 ? 7 : 1;
         double invalidLength = 0.0;
         serializeEnum(out, lattice, 0);
         out & nParameter;
         if (failure == 1) out & invalidLength;
         MemoryIArchive in; in = out;
         bool threw = false;
         try { in & target; } catch (Exception&) { threw = true; }
         in.release();
         TEST_ASSERT(threw);
         TEST_ASSERT(target.isInitialized());
         TEST_ASSERT(target.parameter(0) == 3.0);
         TEST_ASSERT(target.volume() == 3.0);
         TEST_ASSERT(count.n == 0);
      }
   }

   void testEmptyCellRoundTrip()
   {
      printMethod(TEST_FUNC);
      UnitCell<3> source, target;
      MemoryOArchive out; out.allocate(1024);
      int marker = 31415, restoredMarker = 0;
      out & source & marker;
      MemoryIArchive in; in = out;
      in & target & restoredMarker; in.release();
      TEST_ASSERT(restoredMarker == marker);
      TEST_ASSERT(target.lattice() == UnitCell<3>::Null);
      TEST_ASSERT(!target.isInitialized());
      TEST_ASSERT(target.nParameter() == 0);

      // A null archive cannot reset the immutable lattice of a used cell.
      UnitCell<3> initialized;
      FSArray<double,6> lengths;
      lengths.append(1.0); lengths.append(2.0); lengths.append(3.0);
      initialized.set(UnitCell<3>::Orthorhombic, lengths);
      Counter count; Signal<void> signal;
      signal.addObserver(count, &Counter::tick); initialized.setSignal(signal);
      openLogFile("out/UnitCellArchiveNull.log");
      MemoryIArchive second; second = out;
      bool threw = false;
      try { second & initialized; } catch (Exception&) { threw = true; }
      second.release();
      TEST_ASSERT(threw);
      TEST_ASSERT(initialized.isInitialized());
      TEST_ASSERT(initialized.lattice() == UnitCell<3>::Orthorhombic);
      TEST_ASSERT(initialized.volume() == 6.0);
      TEST_ASSERT(initialized.parameter(0) == 1.0);
      TEST_ASSERT(initialized.parameter(1) == 2.0);
      TEST_ASSERT(initialized.parameter(2) == 3.0);
      TEST_ASSERT(count.n == 0);
   }

private:
   struct Counter { int n = 0; void tick() { ++n; } };

   template<int D>
   void check(typename UnitCell<D>::LatticeSystem lattice,
              bool reused, bool text)
   {
      UnitCell<D> source, target;
      FSArray<double,6> lengths, oldLengths;
      double volume = 1.0;
      for (int i = 0; i < D; ++i) {
         lengths.append(i + 2.0); oldLengths.append(1.0);
         volume *= i + 2.0;
      }
      source.set(lattice, lengths);
      if (reused) target.set(lattice, oldLengths);
      Counter count; Signal<void> signal;
      signal.addObserver(count, &Counter::tick); target.setSignal(signal);
      int marker = 31415, restoredMarker = 0;
      if (text) {
         TextFileOArchive out;
         openOutputFile("out/UnitCellArchiveState.txt", out.file());
         out & source & marker; out.file().close();
         TextFileIArchive in;
         openInputFile("out/UnitCellArchiveState.txt", in.file());
         in & target & restoredMarker;
      } else {
         MemoryOArchive out; out.allocate(4096);
         out & source & marker;
         MemoryIArchive in; in = out;
         in & target & restoredMarker; in.release();
      }
      TEST_ASSERT(restoredMarker == marker);
      TEST_ASSERT(target.isInitialized());
      TEST_ASSERT(target.lattice() == lattice);
      TEST_ASSERT(target.nParameter() == D);
      TEST_ASSERT(count.n == 1);
      TEST_ASSERT(target.hasSignal());
      TEST_ASSERT(std::abs(target.volume() - volume) < 1e-12);
      const double pi = Constants::Pi;
      for (int i = 0; i < D; ++i) {
         const double length = i + 2.0;
         IntVec<D> k(0); k[i] = 1;
         TEST_ASSERT(target.parameter(i) == length);
         TEST_ASSERT(source.parameter(i) == length);
         TEST_ASSERT(std::abs(target.ksq(k) - 4*pi*pi/(length*length)) < 1e-12);
         TEST_ASSERT(std::abs(target.dksq(k, i) + 8*pi*pi/(length*length*length)) < 1e-12);
         for (int j = 0; j < D; ++j) {
            TEST_ASSERT(std::abs(target.rBasis(i)[j] - (i == j ? length : 0.0)) < 1e-12);
            TEST_ASSERT(std::abs(target.kBasis(i)[j] - (i == j ? 2*pi/length : 0.0)) < 1e-12);
         }
      }
   }
};

TEST_BEGIN(UnitCellArchiveStateTest)
TEST_ADD(UnitCellArchiveStateTest, testGeometryRoundTrip)
TEST_ADD(UnitCellArchiveStateTest, testInvalidRecordPreservesCell)
TEST_ADD(UnitCellArchiveStateTest, testEmptyCellRoundTrip)
TEST_END(UnitCellArchiveStateTest)

#endif
