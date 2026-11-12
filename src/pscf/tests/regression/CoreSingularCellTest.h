#ifndef PSCF_CORE_SINGULAR_CELL_TEST_H
#define PSCF_CORE_SINGULAR_CELL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/crystal/UnitCell.h>
#include <pscf/crystal/unitCellHeader.h>
#include <util/containers/FSArray.h>
#include <util/math/Constants.h>
#include <util/misc/Exception.h>

#include <iostream>
#include <limits>
#include <sstream>
#include <string>

using namespace Util;
using namespace Pscf;

/*
* Unit cells whose Bravais basis vectors are linearly dependent (or
* nearly so) must be rejected wherever the parameters are set. They
* were accepted, and the minimum image and basis algorithms then failed
* or the programs gave results for a cell of volume ~1e-16.
*/
class CoreSingularCellTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   /*
   * Make a parameter array from a string of values.
   */
   FSArray<double, 6> params(std::string const & values)
   {
      FSArray<double, 6> p;
      std::istringstream in(values);
      double x;
      while (in >> x) {
         p.append(x);
      }
      return p;
   }

   /*
   * Does UnitCell::set(lattice, parameters) throw an Exception?
   */
   template <int D>
   bool setThrows(typename UnitCell<D>::LatticeSystem lattice,
                  FSArray<double, 6> const & p)
   {
      UnitCell<D> cell;
      bool thrown = false;
      try {
         cell.set(lattice, p);
      } catch (Exception& e) {
         thrown = true;
      }
      if (thrown) {
         TEST_ASSERT(!cell.isInitialized());
      } else {
         TEST_ASSERT(cell.isInitialized());
      }
      return thrown;
   }

   /*
   * Does setParameters(p) on a valid cell of this lattice throw?
   */
   template <int D>
   bool setParametersThrows(typename UnitCell<D>::LatticeSystem lattice,
                            FSArray<double, 6> const & valid,
                            FSArray<double, 6> const & p)
   {
      UnitCell<D> cell;
      cell.set(lattice, valid);
      TEST_ASSERT(cell.isInitialized());
      bool thrown = false;
      try {
         cell.setParameters(p);
      } catch (Exception& e) {
         thrown = true;
      }
      return thrown;
   }

   /*
   * Does reading a unit cell from a string with operator >> throw?
   */
   template <int D>
   bool extractThrows(std::string const & text)
   {
      UnitCell<D> cell;
      std::istringstream in(text);
      bool thrown = false;
      try {
         in >> cell;
      } catch (Exception& e) {
         thrown = true;
      }
      return thrown;
   }

   void test1D()
   {
      printMethod(TEST_FUNC);
      double nan = std::numeric_limits<double>::quiet_NaN();
      double inf = std::numeric_limits<double>::infinity();
      FSArray<double, 6> p;

      TEST_ASSERT(!setThrows<1>(UnitCell<1>::Lamellar, params("1.5")));
      TEST_ASSERT(!setThrows<1>(UnitCell<1>::Lamellar, params("1.0e-8")));
      TEST_ASSERT(setThrows<1>(UnitCell<1>::Lamellar, params("0.0")));
      p.clear();
      p.append(nan);
      TEST_ASSERT(setThrows<1>(UnitCell<1>::Lamellar, p));
      p.clear();
      p.append(inf);
      TEST_ASSERT(setThrows<1>(UnitCell<1>::Lamellar, p));
      TEST_ASSERT(extractThrows<1>("lamellar 0.0"));
      TEST_ASSERT(!extractThrows<1>("lamellar 2.0"));
   }

   void test2D()
   {
      printMethod(TEST_FUNC);
      std::string pi = "3.141592653589793";

      // Valid cells, including a strongly skewed one
      TEST_ASSERT(!setThrows<2>(UnitCell<2>::Square, params("1.0")));
      TEST_ASSERT(!setThrows<2>(UnitCell<2>::Rectangular,
                                params("1.0 2.0")));
      TEST_ASSERT(!setThrows<2>(UnitCell<2>::Hexagonal, params("1.0")));
      TEST_ASSERT(!setThrows<2>(UnitCell<2>::Rhombic, params("1.0 1.0")));
      TEST_ASSERT(!setThrows<2>(UnitCell<2>::Oblique,
                                params("1.0 2.0 1.0e-3")));

      // Zero lengths
      TEST_ASSERT(setThrows<2>(UnitCell<2>::Square, params("0.0")));
      TEST_ASSERT(setThrows<2>(UnitCell<2>::Rectangular,
                               params("1.0 0.0")));
      TEST_ASSERT(setThrows<2>(UnitCell<2>::Hexagonal, params("0.0")));
      TEST_ASSERT(setThrows<2>(UnitCell<2>::Rhombic, params("0.0 1.0")));
      TEST_ASSERT(setThrows<2>(UnitCell<2>::Oblique,
                               params("0.0 2.0 1.0")));

      // Parallel basis vectors (gamma = 0 or pi), and nearly parallel
      TEST_ASSERT(setThrows<2>(UnitCell<2>::Rhombic, params("1.0 0.0")));
      TEST_ASSERT(setThrows<2>(UnitCell<2>::Rhombic, params("1.0 " + pi)));
      TEST_ASSERT(setThrows<2>(UnitCell<2>::Oblique,
                               params("1.0 2.0 0.0")));
      TEST_ASSERT(setThrows<2>(UnitCell<2>::Oblique,
                               params("1.0 2.0 " + pi)));
      TEST_ASSERT(setThrows<2>(UnitCell<2>::Oblique,
                               params("1.0 2.0 1.0e-9")));

      // Other paths that set the parameters
      TEST_ASSERT(setParametersThrows<2>(UnitCell<2>::Oblique,
                                         params("1.0 2.0 1.0"),
                                         params("1.0 2.0 0.0")));
      TEST_ASSERT(!setParametersThrows<2>(UnitCell<2>::Oblique,
                                          params("1.0 2.0 1.0"),
                                          params("1.0 2.0 2.0")));
      TEST_ASSERT(extractThrows<2>("oblique 1.0 2.0 " + pi));
      TEST_ASSERT(!extractThrows<2>("oblique 1.0 2.0 2.0"));
   }

   void test3D()
   {
      printMethod(TEST_FUNC);
      std::string pi = "3.141592653589793";
      std::string halfPi = "1.5707963267948966";
      std::string twoPi3 = "2.0943951023931957"; // 2 pi/3

      // Valid cells
      TEST_ASSERT(!setThrows<3>(UnitCell<3>::Cubic, params("1.0")));
      TEST_ASSERT(!setThrows<3>(UnitCell<3>::Tetragonal,
                                params("1.0 2.0")));
      TEST_ASSERT(!setThrows<3>(UnitCell<3>::Orthorhombic,
                                params("1.0 2.0 3.0")));
      TEST_ASSERT(!setThrows<3>(UnitCell<3>::Hexagonal,
                                params("1.0 2.0")));
      TEST_ASSERT(!setThrows<3>(UnitCell<3>::Rhombohedral,
                                params("1.0 1.0")));
      TEST_ASSERT(!setThrows<3>(UnitCell<3>::Monoclinic,
                                params("1.0 2.0 3.0 1.2")));
      TEST_ASSERT(!setThrows<3>(UnitCell<3>::Triclinic,
                                params("1.0 1.0 1.0 0.0 0.2 1.45")));
      TEST_ASSERT(!setThrows<3>(UnitCell<3>::Triclinic,
                                params("1.0 1.0 1.0 0.3 1.57 1.0")));

      // Zero lengths
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Cubic, params("0.0")));
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Tetragonal,
                               params("1.0 0.0")));
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Orthorhombic,
                               params("1.0 0.0 3.0")));
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Hexagonal,
                               params("0.0 2.0")));
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Rhombohedral,
                               params("0.0 1.0")));
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Monoclinic,
                               params("1.0 2.0 0.0 1.2")));
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Triclinic,
                               params("1.0 1.0 0.0 0.0 0.2 1.45")));

      // Rhombohedral: beta = 0 (parallel vectors), 2 pi/3 (coplanar)
      // and beyond 2 pi/3 (no such cell: NaN)
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Rhombohedral,
                               params("1.0 0.0")));
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Rhombohedral,
                               params("1.0 " + twoPi3)));
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Rhombohedral,
                               params("1.0 2.5")));

      // Monoclinic: c parallel to a (beta = 0 or pi)
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Monoclinic,
                               params("1.0 2.0 3.0 0.0")));
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Monoclinic,
                               params("1.0 2.0 3.0 " + pi)));

      // Triclinic: c in the plane of a and b (theta = pi/2), a parallel
      // to b (gamma = 0), and nearly coplanar vectors
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Triclinic,
                         params("1.0 1.0 1.0 0.0 " + halfPi + " 1.45")));
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Triclinic,
                         params("1.0 1.0 1.0 0.3 " + halfPi + " 1.0")));
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Triclinic,
                         params("1.0 1.0 1.0 0.0 0.2 0.0")));
      TEST_ASSERT(setThrows<3>(UnitCell<3>::Triclinic,
                         params("1.0 1.0 1.0 0.0 1.5707963258 1.45")));

      // Other paths that set the parameters
      TEST_ASSERT(setParametersThrows<3>(UnitCell<3>::Triclinic,
                         params("1.0 1.0 1.0 0.0 0.2 1.45"),
                         params("1.0 1.0 1.0 0.0 " + halfPi + " 1.45")));
      TEST_ASSERT(!setParametersThrows<3>(UnitCell<3>::Triclinic,
                         params("1.0 1.0 1.0 0.0 0.2 1.45"),
                         params("1.0 1.0 1.0 0.0 0.3 1.45")));
      TEST_ASSERT(extractThrows<3>("triclinic 3.5 3.5 3.5 0.0 "
                                    + halfPi + " 1.45"));
      TEST_ASSERT(!extractThrows<3>("triclinic 3.5 3.5 3.5 0.0 0.2 1.45"));

      // Unit cell section of a field file header
      std::string header = "crystal_system\n"
                           "              triclinic\n"
                           "N_cell_param\n"
                           "              6\n"
                           "cell_param\n"
                           "    3.5 3.5 3.5 0.0 " + halfPi + " 1.45\n";
      UnitCell<3> cell;
      std::istringstream in(header);
      bool thrown = false;
      try {
         readUnitCellHeader(in, cell);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(thrown);
      TEST_ASSERT(!cell.isInitialized());
   }

};

TEST_BEGIN(CoreSingularCellTest)
TEST_ADD(CoreSingularCellTest, test1D)
TEST_ADD(CoreSingularCellTest, test2D)
TEST_ADD(CoreSingularCellTest, test3D)
TEST_END(CoreSingularCellTest)

#endif
