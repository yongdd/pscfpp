#ifndef RPG_SCFT_BEAD_DS_TEST_H
#define RPG_SCFT_BEAD_DS_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/solvers/Block.h>
#include <rp/solvers/Propagator.h>

#include <pscf/field/cuda/FFT.h>
#include <pscf/field/cuda/WaveList.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/chem/PolymerModel.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

#include <util/containers/FSArray.h>

#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): the contour step ds() of a bead-model
* block is 1 on the GPU, as on the CPU, whatever ds is passed to
* allocate. The GPU block stored the ds passed to allocate.
*/
class ScftBeadDsTest : public UnitTest
{

public:

   void setUp()
   {  PolymerModel::setModel(PolymerModel::Bead); }

   void tearDown()
   {  PolymerModel::setModel(PolymerModel::Thread); }

   /*
   * Allocate a bead-model block of 20 beads with the given ds.
   */
   void checkBeadBlock(double ds)
   {
      Rp::Block<1,CUT> block;
      block.setId(0);
      block.setNBead(20);
      block.setMonomerId(1);
      block.setKuhn(std::sqrt(6.0));

      Mesh<1> mesh;
      IntVec<1> d;
      d[0] = 32;
      mesh.setDimensions(d);
      FFT<1,CUT> fft;
      fft.setup(mesh.dimensions());

      UnitCell<1> unitCell;
      FSArray<double, 6> p;
      p.append(4.0);
      unitCell.set(UnitCell<1>::Lamellar, p);

      WaveList<1,CUT> wavelist;
      wavelist.allocate(mesh, unitCell);

      block.associate(mesh, fft, unitCell, wavelist);
      block.allocate(ds);
      TEST_ASSERT(block.ns() == 22);
      TEST_ASSERT(block.ds() == 1.0);
   }

   void testBeadDs()
   {
      printMethod(TEST_FUNC);
      checkBeadBlock(1.0);
      checkBeadBlock(0.05);
      checkBeadBlock(2.5);
   }

};

TEST_BEGIN(ScftBeadDsTest)
TEST_ADD(ScftBeadDsTest, testBeadDs)
TEST_END(ScftBeadDsTest)

#endif
