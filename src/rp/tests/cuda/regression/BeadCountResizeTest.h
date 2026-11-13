#ifndef RPG_BEAD_COUNT_RESIZE_TEST_H
#define RPG_BEAD_COUNT_RESIZE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <rp/solvers/Block.h>
#include <rp/solvers/Propagator.h>
#include <pscf/field/FFT.h>
#include <pscf/field/WaveList.h>
#include <pscf/field/RField.h>
#include <pscf/backend/HostSendArray.h>
#include <pscf/backend/Reduce.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/chem/PolymerModel.h>
#include <util/tests/LogFileUnitTest.h>
#include <util/misc/Exception.h>
#include <climits>
#include <sstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

class BeadCountResizeTest : public LogFileUnitTest
{
public:
   void tearDown()
   {
      LogFileUnitTest::tearDown();
      PolymerModel::setModel(PolymerModel::Thread);
   }

   void testRepeatedBeadCountChanges()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/BeadCountResize.log");
      for (bool batched : {false, true}) check(batched);
   }

private:
   typedef CUT T;

   void check(bool batched)
   {
      PolymerModel::setModel(PolymerModel::Bead);
      const int nx = 16;
      const double weight = 0.2;
      Mesh<1> mesh;
      IntVec<1> dimensions; dimensions[0] = nx;
      mesh.setDimensions(dimensions);
      FFT<1,T> fft; fft.setup(dimensions);
      UnitCell<1> cell;
      std::istringstream cellInput("lamellar 4.0"); cellInput >> cell;
      WaveList<1,T> waves; waves.allocate(mesh, cell);
      Rp::Block<1,T> block;
      block.setId(0); block.setMonomerId(0); block.setKuhn(1.0);
      Edge& descriptor = block;
      descriptor.setNBead(4);
      descriptor.setNBead(2); // Also exercise updates before allocation.
      block.associate(mesh, fft, cell, waves);
      block.allocate(0.05, batched);
      block.propagator(0).setEndFlags(true, true);
      block.propagator(1).setEndFlags(true, true);

      RField<1,T> w, h0, h1;
      w.allocate(dimensions); h0.allocate(dimensions); h1.allocate(dimensions);
      HostSendArray<typename T::Real,T> hw, ha, hb;
      hw.associate(w); ha.associate(h0); hb.associate(h1);
      for (int i = 0; i < nx; ++i) {
         const double c = std::cos(2.0*std::acos(-1.0)*i/nx);
         hw[i] = weight;
         ha[i] = 1.0 + 0.4*c;
         hb[i] = 1.0 + 0.6*c;
      }
      w = hw; h0 = ha; h1 = hb;
      hw.dissociate(); ha.dissociate(); hb.dissociate();

      // Growth, shrinkage, a single bead, and unchanged counts. The
      // nonuniform supplied heads also give a nonzero stress oracle.
      for (int count : {2, 4, 1, 4, 4, 2}) {
         int oldCount = block.nBead();
         descriptor.setNBead(count); // Verify virtual dispatch via Edge.
         TEST_ASSERT(block.nBead() == count);
         TEST_ASSERT(block.ns() == count + 2);
         for (int direction = 0; direction < 2; ++direction) {
            TEST_ASSERT(block.propagator(direction).ns() == count + 2);
            if (count != oldCount) {
               TEST_ASSERT(!block.propagator(direction).isSolved());
            }
         }
         block.setupSolver(w);
         block.propagator(0).solve(h0);
         block.propagator(1).solve(h1);
         const double k = 2.0*std::acos(-1.0)/4.0;
         const double a = k*k*(count - 1)/6.0;
         const double mode = 0.12*std::exp(-a);
         const double expectedQ = std::exp(-weight*count)*(1.0 + mode);
         double q0, q1;
         block.propagator(0).computeQ(q0);
         block.propagator(1).computeQ(q1);
         TEST_ASSERT(std::abs(q0 - expectedQ) < 1e-11);
         TEST_ASSERT(std::abs(q1 - expectedQ) < 1e-11);
         const double prefactor = 1.0/(count*q0);
         block.computeConcentrationBead(prefactor);
         const Rp::Block<1,T>& constBlock = block;
         TEST_ASSERT(std::abs(Reduce::sum(constBlock.cField())/nx - 1.0) < 1e-11);
         block.computeStressBead(prefactor);
         const double expectedStress = -mode*2.0*a/(4.0*count*(1.0 + mode));
         TEST_ASSERT(std::abs(block.stress(0) - expectedStress) < 1e-11);
      }

      // Invalid and unrepresentable counts must leave the solved chain
      // intact. GPU storage has a further ns*meshSize integer bound.
      for (int invalid : {0, -1, INT_MAX, INT_MAX - 1, INT_MAX/16}) {
         bool threw = false;
         try { descriptor.setNBead(invalid); }
         catch (Exception&) { threw = true; }
         TEST_ASSERT(threw);
         TEST_ASSERT(block.nBead() == 2);
         TEST_ASSERT(block.ns() == 4);
         TEST_ASSERT(block.propagator(0).isSolved());
         TEST_ASSERT(block.propagator(1).isSolved());
      }
   }
};

TEST_BEGIN(BeadCountResizeTest)
TEST_ADD(BeadCountResizeTest, testRepeatedBeadCountChanges)
TEST_END(BeadCountResizeTest)

#endif
