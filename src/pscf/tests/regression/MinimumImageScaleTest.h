#ifndef PSCF_MINIMUM_IMAGE_SCALE_TEST_H
#define PSCF_MINIMUM_IMAGE_SCALE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/chem/Edge.h>
#include <pscf/chem/PolymerModel.h>
#include <pscf/chem/PolymerType.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/field/cpp/FFT.h>
#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldDft.h>
#include <pscf/field/cpp/WaveList.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/containers/FSArray.h>
#include <util/misc/Exception.h>
#include <cmath>
#include <sstream>

using namespace Util;
using namespace Pscf;

class MinimumImageScaleTest : public UnitTest
{
public:
   void setUp() {}
   void tearDown() { PolymerModel::setModel(PolymerModel::Thread); }

   template <int D>
   void checkImages(UnitCell<D> const& unitCell,
                    UnitCell<D> const& scaledCell, double scale)
   {
      IntVec<D> dims(8);
      Mesh<D> mesh(dims);
      WaveList<D,CPT> waves;
      waves.allocate(mesh, scaledCell);
      waves.allowNyquistAverage(false);
      waves.computeMinimumImages();
      waves.computeKSq();
      TEST_ASSERT(waves.kSq()[0] == 0.0);

      MeshIterator<D> iter(dims);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         IntVec<D> v = iter.position();
         IntVec<D> ref = shiftToMinimum(v, dims, unitCell);
         IntVec<D> actual = shiftToMinimum(v, dims, scaledCell);
         TEST_ASSERT(actual == ref);
         double expected = unitCell.ksq(ref);
         double value = scaledCell.ksq(actual)*scale*scale;
         TEST_ASSERT(std::abs(value - expected) < 1e-10*(1.0 + expected));

         // Adding a lattice image must not change the selected alias.
         IntVec<D> alias = v;
         for (int j = 0; j < D; ++j) alias[j] += dims[j];
         TEST_ASSERT(shiftToMinimum(alias, dims, scaledCell) == ref);
      }
   }

   void testMinimumImages2DUnits()
   {
      printMethod(TEST_FUNC);
      FSArray<double,6> p;
      p.append(1.0);
      UnitCell<2> base, scaled;
      base.set(UnitCell<2>::Square, p);
      double scales[] = {1e-6, 1.0, 1e6};
      for (double scale : scales) {
         p[0] = scale;
         scaled.set(UnitCell<2>::Square, p);
         checkImages(base, scaled, scale);
      }
      p[0] = 1.0;
      p.append(1.3);
      p.append(0.7);
      UnitCell<2> oblique, scaledOblique;
      oblique.set(UnitCell<2>::Oblique, p);
      for (double scale : scales) {
         p[0] = scale;
         p[1] = 1.3*scale;
         scaledOblique.set(UnitCell<2>::Oblique, p);
         checkImages(oblique, scaledOblique, scale);
      }
   }

   void testMinimumImages3DUnits()
   {
      printMethod(TEST_FUNC);
      FSArray<double,6> p;
      p.append(1.0);
      UnitCell<3> base, scaled;
      base.set(UnitCell<3>::Cubic, p);
      double scales[] = {1e-6, 1.0, 1e6};
      for (double scale : scales) {
         p[0] = scale;
         scaled.set(UnitCell<3>::Cubic, p);
         checkImages(base, scaled, scale);
      }
      // One Nyquist direction leaves a two-dimensional minimization.
      p[0] = 1.0;
      p.append(1.3);
      p.append(1.7);
      p.append(1.4);
      p.append(0.2);
      p.append(0.1);
      UnitCell<3> triclinic, scaledTriclinic;
      triclinic.set(UnitCell<3>::Triclinic, p);
      for (double scale : scales) {
         p[0] = scale;
         p[1] = 1.3*scale;
         p[2] = 1.7*scale;
         scaledTriclinic.set(UnitCell<3>::Triclinic, p);
         IntVec<3> dims(8);
         for (int j = 0; j < 8; ++j) {
            IntVec<3> v;
            v[0] = 4; v[1] = j; v[2] = 3;
            TEST_ASSERT(shiftToMinimumNyquist(v, dims, scaledTriclinic) ==
                        shiftToMinimumNyquist(v, dims, triclinic));
         }
      }
   }

   // Rounded right angles are physical perturbations, not exact ties.
   // Every inverse pair must still have the same minimum squared norm.
   void testRoundedAngleInversion()
   {
      printMethod(TEST_FUNC);
      double right = 1.5707963268;
      double down = 1.5707963267;
      double scales[] = {1.0, 25.0};
      IntVec<3> dims(4);
      for (double scale : scales) {
         for (int kind = 0; kind < 4; ++kind) {
            FSArray<double,6> p;
            p.append(2.0*scale);
            p.append((kind == 3 ? 2.0 : 2.5)*scale);
            p.append((kind == 3 ? 2.5 : 2.0)*scale);
            UnitCell<3> cell;
            if (kind == 3) {
               p.append(right);
               cell.set(UnitCell<3>::Monoclinic, p);
            } else {
               p.append(kind == 0 ? right : 0.0);
               p.append(0.25);
               p.append(kind == 2 ? down : right);
               cell.set(UnitCell<3>::Triclinic, p);
            }
            MeshIterator<3> iter(dims);
            for (iter.begin(); !iter.atEnd(); ++iter) {
               IntVec<3> v = iter.position();
               IntVec<3> inverse;
               for (int j = 0; j < 3; ++j) {
                  inverse[j] = (dims[j] - v[j]) % dims[j];
               }
               IntVec<3> a = shiftToMinimum(v, dims, cell);
               IntVec<3> b = shiftToMinimum(inverse, dims, cell);
               double norm = cell.ksq(a);
               double inverseNorm = cell.ksq(b);
               TEST_ASSERT(std::abs(norm-inverseNorm) < 1e-8);
               TEST_ASSERT(std::abs(norm-inverseNorm) < 1e-12*(1.0+norm));
            }
         }
      }
   }

};

TEST_BEGIN(MinimumImageScaleTest)
TEST_ADD(MinimumImageScaleTest, testMinimumImages2DUnits)
TEST_ADD(MinimumImageScaleTest, testMinimumImages3DUnits)
TEST_ADD(MinimumImageScaleTest, testRoundedAngleInversion)
TEST_END(MinimumImageScaleTest)

#endif
