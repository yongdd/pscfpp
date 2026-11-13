#ifndef RPC_BEAD_HEAD_PARTITION_TEST_H
#define RPC_BEAD_HEAD_PARTITION_TEST_H
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
#include <util/misc/Log.h>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace Util;
using namespace Pscf;


class BeadHeadPartitionTest : public UnitTest
{
public:
   void tearDown() { PolymerModel::setModel(PolymerModel::Thread); }

   void testProvidedHeads()
   {
      printMethod(TEST_FUNC);
      for (int mode=0; mode<3; ++mode) {
         check(true, 1, mode, 0.0);
         check(true, 4, mode, 0.0);
         check(false, 0, mode, 0.0);
      }
   }

   void testNonzeroBeadWeight()
   {
      printMethod(TEST_FUNC);
      check(true, 1, 2, 0.7);
      check(true, 4, 1, -0.4);
   }

private:
   typedef CPT T;
   void check(bool bead, int count, int mode, double weight)
   {
      PolymerModel::setModel(bead ? PolymerModel::Bead : PolymerModel::Thread);
      const int nx=16;
      const double length=3.0, kuhn=1.2, cellLength=4.0;
      Mesh<1> mesh; IntVec<1> dimensions; dimensions[0]=nx;
      mesh.setDimensions(dimensions);
      FFT<1,T> fft; fft.setup(dimensions);
      UnitCell<1> cell; std::istringstream cellInput("lamellar 4.0");
      cellInput >> cell;
      WaveList<1,T> waveList; waveList.allocate(mesh,cell);
      Rp::Block<1,T> block;
      block.setId(0); block.setMonomerId(0); block.setKuhn(kuhn);
      if (bead) block.setNBead(count); else block.setLength(length);
      block.associate(mesh,fft,cell,waveList); block.allocate(.05,false);
      block.propagator(0).setEndFlags(true,true);
      block.propagator(1).setEndFlags(true,true);

      RField<1,T> w,h0,h1;
      w.allocate(dimensions); h0.allocate(dimensions); h1.allocate(dimensions);
      {
         HostSendArray<typename T::Real,T> hw,ha,hb;
         hw.associate(w); ha.associate(h0); hb.associate(h1);
         for(int i=0;i<nx;++i) {
            double c=std::cos(2.0*std::acos(-1.0)*i/nx);
            hw[i]=weight;
            ha[i]=(mode==1 ? 2.0 : (mode==2 ? 1.0+.4*c : 1.0));
            hb[i]=(mode==1 ? 3.0 : (mode==2 ? 1.0+.6*c : 1.0));
         }
         w=hw; h0=ha; h1=hb;
         hw.dissociate(); ha.dissociate(); hb.dissociate();
      }
      block.setupSolver(w);
      auto& p0=block.propagator(0); auto& p1=block.propagator(1);
      if (mode==0) { p0.solve(); p1.solve(); }
      else { p0.solve(h0); p1.solve(h1); }

      // At w=0 the only nonconstant mode is cos(2*pi*x/L).  Its heat
      // kernel eigenvalue gives an independent exact partition oracle.
      const double k=2.0*std::acos(-1.0)/cellLength;
      const double transfer=std::exp(-k*k*kuhn*kuhn/6.0*(bead?count-1:length));
      const double oracle=std::exp(-weight*(bead?count:length))
         *(mode==1?6.0:(mode==2?1.0+.12*transfer:1.0));
      double q0,q1;
      p0.computeQ(q0); p1.computeQ(q1);
      double sliceError=0.0;
      const int begin=bead?1:0, end=bead?block.ns()-1:block.ns();
      for(int j=begin;weight==0.0 && j<end;++j) {
         double q=Reduce::innerProduct(p0.q(j),p1.q(block.ns()-1-j))/nx;
         sliceError=std::max(sliceError,std::abs(q-oracle));
      }
      // Verify that the propagated slices themselves remain correct, then
      // quantify the normalization error caused by using public computeQ.
      if (bead) block.computeConcentrationBead(1.0/(count*q0));
      else block.computeConcentrationThread(1.0/(length*q0));
      const double mass=Reduce::sum(block.cField())/nx;

      TEST_ASSERT(std::abs(q0-oracle) < 1e-11);
      TEST_ASSERT(std::abs(q1-oracle) < 1e-11);
      TEST_ASSERT(sliceError < 1e-11);
      TEST_ASSERT(std::abs(mass-1.0) < 1e-11);
   }
};
TEST_BEGIN(BeadHeadPartitionTest)
TEST_ADD(BeadHeadPartitionTest, testProvidedHeads)
TEST_ADD(BeadHeadPartitionTest, testNonzeroBeadWeight)
TEST_END(BeadHeadPartitionTest)
#endif
