#ifndef PSCF_CORE_NYQUIST_PLANE_CACHE_TEST_H
#define PSCF_CORE_NYQUIST_PLANE_CACHE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/containers/FSArray.h>

#include <cmath>
#include <vector>

using namespace Util;
using namespace Pscf;

/*
* shiftToMinimumNyquist keeps the reduced neighbors of the plane
* lattices of image shifts of the last mesh and cell (per thread), so
* that a WaveList does not repeat the reduction for every wavevector.
* Calls that alternate between cells and meshes must give the images
* of calls made for one cell and mesh at a time, and these images
* minimize ksqNyquist (checked by a search over images).
*/
class CoreNyquistPlaneCacheTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void makeCell(UnitCell<3>& cell, UnitCell<3>::LatticeSystem lattice,
                 double const * p, int nParam)
   {
      FSArray<double, 6> parameters;
      for (int i = 0; i < nParam; ++i) {
         parameters.append(p[i]);
      }
      cell.set(lattice, parameters);
   }

   /*
   * Images of all wavevectors with a Nyquist index, for one cell.
   */
   void images(UnitCell<3> const & cell, IntVec<3> const & d,
               std::vector< IntVec<3> >& waves,
               std::vector< IntVec<3> >& result)
   {
      waves.clear();
      result.clear();
      MeshIterator<3> iter(d);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         if (hasNyquistIndex(iter.position(), d)) {
            waves.push_back(iter.position());
            result.push_back(shiftToMinimumNyquist(iter.position(),
                                                   d, cell));
         }
      }
   }

   /*
   * Minimum of ksqNyquist over the images v + s*d, |s[i]| <= 3.
   */
   double searchMinimum(IntVec<3> const & v, IntVec<3> const & d,
                        UnitCell<3> const & cell)
   {
      double best = -1.0;
      IntVec<3> u;
      for (int s0 = -3; s0 <= 3; ++s0) {
         for (int s1 = -3; s1 <= 3; ++s1) {
            for (int s2 = -3; s2 <= 3; ++s2) {
               u[0] = v[0] + s0*d[0];
               u[1] = v[1] + s1*d[1];
               u[2] = v[2] + s2*d[2];
               double g = ksqNyquist(u, d, cell);
               if (best < 0.0 || g < best) best = g;
            }
         }
      }
      return best;
   }

   void testAlternatingCells()
   {
      printMethod(TEST_FUNC);

      UnitCell<3> cells[3];
      double pa[6] = {1.9, 2.1, 2.3, 1.3, 0.5, 1.2};
      makeCell(cells[0], UnitCell<3>::Triclinic, pa, 6);
      double pb[6] = {2.0, 1.8, 2.4, 0.3, 1.1, 1.9};
      makeCell(cells[1], UnitCell<3>::Triclinic, pb, 6);
      double pc[4] = {1.9, 2.1, 2.3, 2.0};
      makeCell(cells[2], UnitCell<3>::Monoclinic, pc, 4);
      IntVec<3> dims[2];
      dims[0][0] = 6; dims[0][1] = 8; dims[0][2] = 10;
      dims[1] = IntVec<3>(8);

      // Images for one cell and mesh at a time
      std::vector< IntVec<3> > waves[3][2];
      std::vector< IntVec<3> > reference[3][2];
      for (int c = 0; c < 3; ++c) {
         for (int m = 0; m < 2; ++m) {
            images(cells[c], dims[m], waves[c][m], reference[c][m]);
            TEST_ASSERT(waves[c][m].size() > 0);
         }
      }

      // Each image minimizes ksqNyquist
      for (int c = 0; c < 3; ++c) {
         for (int m = 0; m < 2; ++m) {
            IntVec<3> const & d = dims[m];
            for (size_t i = 0; i < waves[c][m].size(); ++i) {
               IntVec<3> const & image = reference[c][m][i];
               double g = ksqNyquist(image, d, cells[c]);
               double best = searchMinimum(waves[c][m][i], d, cells[c]);
               TEST_ASSERT(std::abs(g - best) < 1.0E-10*(1.0 + best));
               for (int j = 0; j < 3; ++j) {
                  TEST_ASSERT((image[j] - waves[c][m][i][j]) % d[j] == 0);
               }
            }
         }
      }

      // Calls alternating between cells and meshes give the same images
      for (int m = 0; m < 2; ++m) {
         for (size_t i = 0; i < waves[0][m].size(); ++i) {
            for (int c = 0; c < 3; ++c) {
               int mm = (c == 1) ? 1 - m : m;
               if (i < waves[c][mm].size()) {
                  IntVec<3> image = shiftToMinimumNyquist(waves[c][mm][i],
                                                    dims[mm], cells[c]);
                  TEST_ASSERT(image == reference[c][mm][i]);
               }
            }
         }
      }
   }

};

TEST_BEGIN(CoreNyquistPlaneCacheTest)
TEST_ADD(CoreNyquistPlaneCacheTest, testAlternatingCells)
TEST_END(CoreNyquistPlaneCacheTest)

#endif
