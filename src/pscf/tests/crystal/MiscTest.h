#ifndef PSCF_MISC_TEST_H
#define PSCF_MISC_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/field/cpp/FFT.h>
#include <pscf/crystal/sortWaves.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/Sort.h>
#include <util/math/Constants.h>
#include <util/containers/FSArray.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <vector>
#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class MiscTest : public UnitTest 
{

public:

   void setUp()
   {  setVerbose(0); }

   void tearDown()
   {  setVerbose(0); }
 
   void testSortWavesReal() 
   {
      printMethod(TEST_FUNC);

      UnitCell<3> cell;
      FSArray<double, 6> parameters;
      double a = 2.0;
      parameters.append(a);
      cell.set(UnitCell<3>::Cubic, parameters);

      IntVec<3> dimensions;
      dimensions[0] = 8;
      dimensions[1] = 8;
      dimensions[2] = 8;
      IntVec<3> kDimensions;
      int kSize;
      FFT<3,CPT>::computeKMesh(dimensions, kDimensions, kSize);
      Mesh<3> kMesh(kDimensions);

      std::vector< Sort::Item<double> > items;
      GArray< Sort::Bunch > bunches;
      double epsilon = 1.0E-8;
      bool isRealField = true;
      sortWaves(cell, dimensions, items, bunches, epsilon, isRealField);

      const double twoPi = 2.0 * Constants::Pi;
      const double bSqInv = (a*a)/(twoPi * twoPi);

      double oldVal, value;
      int nb = bunches.size();
      int begin, end, size, ib, iw, rank, itot;
      IntVec<3> position;
      itot = 0; 
      for (ib = 0; ib < nb; ++ib) {
         begin = bunches[ib][0];
         end = bunches[ib][1];
         if (ib > 0) {
           TEST_ASSERT(begin == bunches[ib-1][1]);
         }
         oldVal = items[begin].value;
         for (iw = begin; iw < end; ++iw) {
             value = items[iw].value;
             TEST_ASSERT( std::abs( value - oldVal) < epsilon);
             if (iw > 0) {
                TEST_ASSERT(value >= items[iw-1].value);
             } 
             ++itot;
         }
         if (verbose() > 0) {
            size = end - begin;
            std::cout << std::endl << Int(begin) << Int(end) << Int(size);
            for (iw = begin; iw < end; ++iw) {
               value = items[iw].value;
               value *= bSqInv;
               rank = items[iw].id;
               position = kMesh.position(rank);
               position = shiftToMinimum(position, dimensions, cell);
               std::cout << std::endl 
                         << Int(iw)
                         << Int(rank)
                         << Dbl(value);
               for (int k=0; k < 3; ++k) {
                  std::cout << Int(position[k]);
               }
            }
         }
      }
      TEST_ASSERT(itot == (int) items.size());
   }

   void testSortWavesComplex() 
   {
      printMethod(TEST_FUNC);

      UnitCell<3> cell;
      FSArray<double, 6> parameters;
      double a = 2.0;
      parameters.append(a);
      cell.set(UnitCell<3>::Cubic, parameters);

      IntVec<3> dimensions;
      dimensions[0] = 8;
      dimensions[1] = 8;
      dimensions[2] = 8;
      IntVec<3> kDimensions = dimensions;
      Mesh<3> kMesh(kDimensions);

      std::vector< Sort::Item<double> > items;
      GArray< Sort::Bunch > bunches;
      double epsilon = 1.0E-8;
      bool isRealField = false;
      sortWaves(cell, dimensions, items, bunches, epsilon, isRealField);

      const double twoPi = 2.0 * Constants::Pi;
      const double bSqInv = (a*a)/(twoPi * twoPi);

      double oldVal, value;
      int nb = bunches.size();
      int begin, end, size, ib, iw, rank, itot;
      IntVec<3> position;
      itot = 0; 
      for (ib = 0; ib < nb; ++ib) {
         begin = bunches[ib][0];
         end = bunches[ib][1];
         if (ib > 0) {
           TEST_ASSERT(begin == bunches[ib-1][1]);
         }
         oldVal = items[begin].value;
         for (iw = begin; iw < end; ++iw) {
             value = items[iw].value;
             TEST_ASSERT( std::abs( value - oldVal) < epsilon);
             if (iw > 0) {
                TEST_ASSERT(value >= items[iw-1].value);
             } 
             ++itot;
         }
         if (verbose() > 0) {
            size = end - begin;
            std::cout << std::endl << Int(begin) << Int(end) << Int(size);
            for (iw = begin; iw < end; ++iw) {
               value = items[iw].value;
               value *= bSqInv;
               rank = items[iw].id;
               position = kMesh.position(rank);
               position = shiftToMinimum(position, dimensions, cell);
               std::cout << std::endl 
                         << Int(iw)
                         << Int(rank)
                         << Dbl(value);
               for (int k=0; k < 3; ++k) {
                  std::cout << Int(position[k]);
               }
            }
         }
      }
      TEST_ASSERT(itot == (int) items.size());
   }

   /*
   * Skewed cells with strongly anisotropic meshes, where the fixed
   * search window of shiftToMinimum alone misses minimum images.
   */
   void testShiftToMinimumOblique()
   {
      printMethod(TEST_FUNC);

      const double pi = Constants::Pi;
      double gammas[4] = {2.0, 2.4, 2.8, pi - 2.4};
      int meshes[4][2] = {{4, 64}, {8, 64}, {1, 16}, {64, 4}};
      UnitCell<2> cell;
      IntVec<2> d;
      for (int ig = 0; ig < 4; ++ig) {
         FSArray<double, 6> parameters;
         parameters.append(2.0);
         parameters.append(2.0);
         parameters.append(gammas[ig]);
         cell.set(UnitCell<2>::Oblique, parameters);
         for (int im = 0; im < 4; ++im) {
            d[0] = meshes[im][0];
            d[1] = meshes[im][1];
            checkMinimumImages<2>(cell, d, 24);
         }
      }
   }

   void testShiftToMinimumMonoclinic()
   {
      printMethod(TEST_FUNC);

      double betas[2] = {2.4, 2.8};
      int meshes[3][3] = {{4, 2, 32}, {32, 2, 4}, {8, 4, 32}};
      UnitCell<3> cell;
      IntVec<3> d;
      for (int ib = 0; ib < 2; ++ib) {
         FSArray<double, 6> parameters;
         parameters.append(2.0);
         parameters.append(2.0);
         parameters.append(2.0);
         parameters.append(betas[ib]);
         cell.set(UnitCell<3>::Monoclinic, parameters);
         for (int im = 0; im < 3; ++im) {
            for (int i = 0; i < 3; ++i) d[i] = meshes[im][i];
            checkMinimumImages<3>(cell, d, 12);
         }
      }
   }

private:

   /*
   * Check that shiftToMinimum returns, for every mesh point v, an image
   * of v with the smallest |k|^2 among the images v + s*d with all
   * |s[i]| <= q, and that v and -v get images of equal |k|^2.
   */
   template <int D>
   void checkMinimumImages(UnitCell<D> const & cell, IntVec<D> const & d,
                           int q)
   {
      const double epsilon = 1.0E-8;
      IntVec<D> sDims;
      for (int i = 0; i < D; ++i) sDims[i] = 2*q + 1;
      MeshIterator<D> vItr(d);
      MeshIterator<D> sItr(sDims);
      IntVec<D> v, r, u, vNeg, rNeg;
      for (vItr.begin(); !vItr.atEnd(); ++vItr) {
         v = vItr.position();
         r = shiftToMinimum(v, d, cell);
         for (int i = 0; i < D; ++i) {
            TEST_ASSERT((r[i] - v[i]) % d[i] == 0);
         }
         double gMin = cell.ksq(v);
         for (sItr.begin(); !sItr.atEnd(); ++sItr) {
            for (int i = 0; i < D; ++i) {
               u[i] = v[i] + (sItr.position(i) - q)*d[i];
            }
            double g = cell.ksq(u);
            if (g < gMin) gMin = g;
         }
         TEST_ASSERT(cell.ksq(r) < gMin + epsilon);
         vNeg.negate(v);
         rNeg = shiftToMinimum(vNeg, d, cell);
         TEST_ASSERT(std::abs(cell.ksq(rNeg) - cell.ksq(r)) < epsilon);

         // The result depends only on the class of v, not on the
         // representative passed in (Basis passes minimum images and
         // their negatives, as well as DFT indices)
         TEST_ASSERT(shiftToMinimum(r, d, cell) == r);
         u.negate(r);
         TEST_ASSERT(shiftToMinimum(u, d, cell) == rNeg);
      }
   }

};

TEST_BEGIN(MiscTest)
TEST_ADD(MiscTest, testSortWavesReal)
TEST_ADD(MiscTest, testSortWavesComplex)
TEST_ADD(MiscTest, testShiftToMinimumOblique)
TEST_ADD(MiscTest, testShiftToMinimumMonoclinic)
TEST_END(MiscTest)

#endif
