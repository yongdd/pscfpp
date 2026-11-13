#ifndef PSCF_ANISOTROPIC_MIRROR_MESH_TEST_H
#define PSCF_ANISOTROPIC_MIRROR_MESH_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/crystal/SpaceGroup.h>
#include <pscf/crystal/Basis.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <util/containers/FSArray.h>
#include <util/misc/Exception.h>
#include <cmath>
#include <limits>
#include <set>

using namespace Util;
using namespace Pscf;

/* Mesh invariance is divisibility, not equality of mixed-axis counts. */
class AnisotropicMirrorMeshTest : public UnitTest
{
public:

   bool accepts(SpaceGroup<2> const & group, int nx, int ny)
   {
      IntVec<2> dimensions;
      dimensions[0] = nx;
      dimensions[1] = ny;
      try {
         group.checkMeshDimensions(dimensions);
      } catch (Exception&) {
         return false;
      }
      return true;
   }

   void checkBasis(bool oblique, int nx, int ny, int expected)
   {
      SpaceGroup<2> group;
      SpaceSymmetry<2> mirror;
      mirror.R(0,1) = -1;
      mirror.R(1,1) = -1;
      group.add(mirror);
      TEST_ASSERT(group.isValid());
      UnitCell<2> cell;
      FSArray<double,6> p;
      p.append(2.0);
      if (oblique) {
         // Cartesian vectors (2,0), (-1,sqrt(8)); horizontal reflection.
         p.append(3.0);
         p.append(std::acos(-1.0/3.0));
         cell.set(UnitCell<2>::Oblique, p);
      } else {
         cell.set(UnitCell<2>::Hexagonal, p);
      }
      // Independently check the physical metric and exact node permutation.
      for (int i = -4; i <= 4; ++i) {
         for (int j = -4; j <= 4; ++j) {
            IntVec<2> k, image;
            k[0] = i; k[1] = j;
            image[0] = i; image[1] = -i-j;
            TEST_ASSERT(std::abs(cell.ksq(k)-cell.ksq(image)) < 1e-10);
         }
      }
      std::set<int> nodes;
      for (int i = 0; i < nx; ++i) {
         for (int j = 0; j < ny; ++j) {
            TEST_ASSERT((j*nx) % ny == 0);
            int x = ((i-j*nx/ny) % nx + nx) % nx;
            int y = (ny-j) % ny;
            nodes.insert(x*ny+y);
         }
      }
      TEST_ASSERT(int(nodes.size()) == nx*ny);
      IntVec<2> dimensions;
      dimensions[0] = nx;
      dimensions[1] = ny;
      Mesh<2> mesh(dimensions);
      Basis<2> basis;
      basis.makeBasis(mesh, cell, group);
      TEST_ASSERT(basis.isValid());
      // Burnside: (node count + number of mirror-fixed nodes)/2.
      TEST_ASSERT(basis.nBasis() == expected);
   }

   void testAnisotropicBasis()
   {
      printMethod(TEST_FUNC);
      checkBasis(true, 16, 8, 72);
      checkBasis(false, 16, 8, 72);
      checkBasis(true, 8, 8, 36);
   }

   void testNonInvariantControls()
   {
      printMethod(TEST_FUNC);
      SpaceGroup<2> mirrorGroup;
      SpaceSymmetry<2> mirror;
      mirror.R(0,1) = -1;
      mirror.R(1,1) = -1;
      mirrorGroup.add(mirror);
      TEST_ASSERT(!accepts(mirrorGroup, 8, 16));
      TEST_ASSERT(!accepts(mirrorGroup, 12, 8));
      TEST_ASSERT(accepts(mirrorGroup, 24, 8));
      // A genuine exchange of axes still requires equal counts.
      SpaceGroup<2> swapGroup;
      SpaceSymmetry<2> swap;
      swap.R(0,0) = 0; swap.R(1,1) = 0;
      swap.R(0,1) = 1; swap.R(1,0) = 1;
      swapGroup.add(swap);
      TEST_ASSERT(swapGroup.isValid());
      TEST_ASSERT(!accepts(swapGroup, 16, 8));
      TEST_ASSERT(accepts(swapGroup, 8, 8));
      // Translation compatibility must be retained.
      SpaceGroup<2> translated;
      SpaceSymmetry<2> halfShift;
      halfShift.t(0) = Rational(1,2);
      translated.add(halfShift);
      TEST_ASSERT(translated.isValid());
      TEST_ASSERT(accepts(translated, 16, 9));
      TEST_ASSERT(!accepts(translated, 15, 9));
   }

   void testLargeIntegerDimensions()
   {
      printMethod(TEST_FUNC);
      SpaceGroup<2> group;
      SpaceSymmetry<2> mirror;
      mirror.R(0,1) = -2;
      mirror.R(1,1) = -1;
      group.add(mirror);
      TEST_ASSERT(group.isValid());
      // Only check arithmetic: these tests do not allocate enormous meshes.
      // Multiplying nx by R(0,1) would overflow a signed int.
      int largest = std::numeric_limits<int>::max();
      TEST_ASSERT(accepts(group, largest, largest));
      TEST_ASSERT(accepts(group, largest, 2));
      TEST_ASSERT(!accepts(group, largest, largest-1));
      TEST_ASSERT(!accepts(group, 0, 8));
      TEST_ASSERT(!accepts(group, 8, -1));
   }
};

TEST_BEGIN(AnisotropicMirrorMeshTest)
TEST_ADD(AnisotropicMirrorMeshTest, testAnisotropicBasis)
TEST_ADD(AnisotropicMirrorMeshTest, testNonInvariantControls)
TEST_ADD(AnisotropicMirrorMeshTest, testLargeIntegerDimensions)
TEST_END(AnisotropicMirrorMeshTest)

#endif
