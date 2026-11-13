#ifndef PSCF_WAVE_SHELL_SCALE_CHECKS_H
#define PSCF_WAVE_SHELL_SCALE_CHECKS_H

#include <test/UnitTest.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/field/WaveList.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/containers/FSArray.h>
#include <limits>
#include <map>
#include <vector>

namespace Pscf {

   template <class T>
   void checkWaveShellScaleCubicShells()
   {
      IntVec<3> dims(8), kdims(8);
      kdims[2] = 5;
      Mesh<3> mesh(dims);
      std::vector<int> norms;
      std::map<int,int> counts;
      MeshIterator<3> iter(kdims);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         int norm = 0;
         for (int j = 0; j < 3; ++j) {
            int k = iter.position()[j];
            if (k > 4) k = 8-k;
            norm += k*k;
         }
         norms.push_back(norm);
         ++counts[norm];
      }
      TEST_ASSERT(counts.size() == 32);
      double scales[] = {1e-6, .01, 1.0, 1e4, 1e5, 1e6};
      for (double scale : scales) {
         FSArray<double,6> p;
         p.append(scale);
         UnitCell<3> cell;
         cell.set(UnitCell<3>::Cubic, p);
         WaveList<3,T> waves;
         waves.allocate(mesh, cell);
         waves.sortWaves();
         TEST_ASSERT(waves.nBunch() == int(counts.size()));
         TEST_ASSERT(waves.sortedBunches()[0][1] == 1);
         TEST_ASSERT(waves.sortedIds()[0] == 0);
         int ib = 0;
         for (auto const & expected : counts) {
            auto const & bunch = waves.sortedBunches()[ib];
            TEST_ASSERT(bunch[1]-bunch[0] == expected.second);
            for (int i = bunch[0]; i < bunch[1]; ++i) {
               int rank = waves.sortedIds()[i];
               TEST_ASSERT(norms[rank] == expected.first);
               TEST_ASSERT(waves.bunchIds()[rank] == ib);
            }
            ++ib;
         }
      }
   }

   template <class T>
   void checkWaveShellScaleDistinctShells()
   {
      // Near-square cells must retain physically distinct shells. The
      // double-precision perturbation matches the earlier W4 regression.
      double delta = std::numeric_limits<typename T::Real>::digits < 50
                   ? 2.0E-4 : 2.0E-9;
      IntVec<2> dims(4);
      Mesh<2> mesh(dims);
      double scales[] = {1e-6, 1.0, 1e6};
      for (double scale : scales) {
         FSArray<double,6> p;
         p.append(3.2*scale);
         p.append(3.2*(1.0+delta)*scale);
         UnitCell<2> cell;
         cell.set(UnitCell<2>::Rectangular, p);
         WaveList<2,T> waves;
         waves.allocate(mesh, cell);
         waves.sortWaves();
         TEST_ASSERT(waves.nBunch() == 9);
         TEST_ASSERT(waves.sortedBunches()[0][1] == 1);
         // The y and x fundamentals have ranks 1 and 3 in the 4x3 k-grid.
         TEST_ASSERT(waves.bunchIds()[1] != waves.bunchIds()[3]);
      }
   }

}
#endif
