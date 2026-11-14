#ifndef PSCF_CPP_WAVE_LIST_TEST_H
#define PSCF_CPP_WAVE_LIST_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/field/cpp/WaveList.h>
#include <pscf/field/cpp/RField.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/DArray.h>
#include <util/math/Constants.h>

#include <cmath>

using namespace Util;
using namespace Pscf;

class CppWaveListTest : public UnitTest
{

private:

   // Error tolerance for array equality
   constexpr static double tolerance_ = 1E-10;

   Mesh<1> mesh1;
   Mesh<2> mesh2;
   Mesh<3> mesh3;

   UnitCell<1> cell1;
   UnitCell<2> cell2;
   UnitCell<3> cell3;

   IntVec<1> meshDims1;
   IntVec<2> meshDims2;
   IntVec<3> meshDims3;

   IntVec<1> kMeshDims1;
   IntVec<2> kMeshDims2;
   IntVec<3> kMeshDims3;

   int kSize1, kSize2, kSize3;

public:

   void setUp()
   {
      // set up 1D mesh
      // IntVec<1> meshDims1; 
      meshDims1[0] = 32; 
      mesh1.setDimensions(meshDims1);
      kMeshDims1[0] = (meshDims1[0] / 2) + 1;
      kSize1 = kMeshDims1[0];

      // Set up 2D meshes
      //IntVec<2> meshDims2; 
      meshDims2[0] = 32; 
      meshDims2[1] = 48; 
      mesh2.setDimensions(meshDims2);
      kMeshDims2 = meshDims2;
      kMeshDims2[1] = meshDims2[1] / 2 + 1;
      kSize2 = kMeshDims2[0] * kMeshDims2[1];

      // Set up 3D meshes
      //IntVec<3> meshDims3; 
      meshDims3[0] = 24; 
      meshDims3[1] = 64; 
      meshDims3[2] = 21;
      mesh3.setDimensions(meshDims3);
      kMeshDims3 = meshDims3;
      kMeshDims3[2] = meshDims3[2] / 2 + 1;
      kSize3 = kMeshDims3[0] * kMeshDims3[1] * kMeshDims3[2];

      // set up 1D, 2D, and 3D unit cell objects
      std::ifstream in1, in2, in3;
      openInputFile("in/Lamellar", in1);
      in1 >> cell1;

      openInputFile("in/Oblique", in2);
      in2 >> cell2;

      openInputFile("in/Triclinic", in3);
      in3 >> cell3;
   }

   void tearDown()
   {}

   /*
   * Reference for a 2D wavevector v whose kSq the WaveList averages
   * over the sign of its Nyquist index (a lattice system with a variable
   * angle, an index d/2, two nonzero indices; see ksqNyquist).
   *
   * Then kSq = (1 - w)*|k(r)|^2 + w*A, where r is the minimum image, A
   * is the average of |k|^2 of the image u with indices reduced to
   * (-d/2, d/2] and of u with the sign of a Nyquist index reversed, and
   * w is a smooth step of the absolute cosine t of the angle between the
   * reciprocal basis vectors: 1 for t <= 0.1, 0 for t >= 0.25, and
   * 1 - x^3 (10 - 15x + 6x^2) with x = (t - 0.1)/0.15 between. The
   * stored image is u if w >= 1/2, r otherwise. If dksq is not null,
   * dksq[n] is set to the derivative of kSq with respect to cell 
   * parameter n (without the factor 2 of implicit wavevectors).
   *
   * Returns false and sets nothing if w = 0 or if v has no Nyquist
   * index or a zero index (kSq is then |k|^2 of the minimum image).
   */
   bool nyquistReference(IntVec<2> const & v, IntVec<2> const & d,
                         UnitCell<2> const & cell, IntVec<2>& image,
                         double& ksq, double* dksq)
   {
      // Reduced indices u, and u2 = u with a Nyquist index reversed
      IntVec<2> u, u2;
      int nyquist = -1;
      for (int i = 0; i < 2; ++i) {
         u[i] = ((v[i] % d[i]) + d[i]) % d[i];
         if (2*u[i] > d[i]) u[i] -= d[i];
         if (2*u[i] == d[i]) nyquist = i;
      }
      if (nyquist < 0 || u[0] == 0 || u[1] == 0) return false;
      u2 = u;
      u2[nyquist] = -u[nyquist];

      // Weight w, and its derivative dw with respect to t
      double g[2][2];
      for (int i = 0; i < 2; ++i) {
         for (int j = 0; j < 2; ++j) {
            g[i][j] = cell.kBasis(i)[0]*cell.kBasis(j)[0]
                    + cell.kBasis(i)[1]*cell.kBasis(j)[1];
         }
      }
      double norm = std::sqrt(g[0][0]*g[1][1]);
      double c = g[0][1]/norm;
      double t = std::abs(c);
      if (t >= 0.25) return false;
      double w = 1.0;
      double dw = 0.0;
      if (t > 0.1) {
         double x = (t - 0.1)/0.15;
         w = 1.0 - x*x*x*(10.0 - 15.0*x + 6.0*x*x);
         dw = -30.0*x*x*(1.0 - x)*(1.0 - x)/0.15;
      }

      IntVec<2> r = shiftToMinimum(v, d, cell);
      double kMin = cell.ksq(r);
      double kAverage = 0.5*(cell.ksq(u) + cell.ksq(u2));
      ksq = (1.0 - w)*kMin + w*kAverage;
      image = (w >= 0.5) ? u : r;
      if (dksq) {
         for (int n = 0; n < cell.nParameter(); ++n) {
            double dc = cell.dkkBasis(n, 0, 1)/norm
                      - 0.5*c*( cell.dkkBasis(n, 0, 0)/g[0][0]
                              + cell.dkkBasis(n, 1, 1)/g[1][1] );
            double dt = (c < 0.0) ? -dc : dc;
            dksq[n] = (1.0 - w)*cell.dksq(r, n)
                    + 0.5*w*(cell.dksq(u, n) + cell.dksq(u2, n))
                    + (kAverage - kMin)*dw*dt;
         }
      }
      return true;
   }

   void testAllocate()
   {
      printMethod(TEST_FUNC);

      WaveList<1,CPT> wavelist1;
      wavelist1.allocate(mesh1, cell1);
      TEST_ASSERT(wavelist1.isAllocated());
      TEST_ASSERT(!wavelist1.hasMinImages());

      WaveList<2,CPT> wavelist2;
      wavelist2.allocate(mesh2, cell2);
      TEST_ASSERT(wavelist2.isAllocated());
      TEST_ASSERT(!wavelist2.hasMinImages());

      WaveList<3,CPT> wavelist3;
      wavelist3.allocate(mesh3, cell3);
      TEST_ASSERT(wavelist3.isAllocated());
      TEST_ASSERT(!wavelist3.hasMinImages());
   }

   void testComputeMinimumImages1D()
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(mesh1.dimension(0) == 32);
      TEST_ASSERT(kMeshDims1[0] == 17);
      TEST_ASSERT(kSize1 == 17);

      // set up wavelist object
      WaveList<1,CPT> wavelist;
      wavelist.allocate(mesh1, cell1);

      // Compute minimum images (and ksq)
      wavelist.computeMinimumImages(); 
      DArray< IntVec<1> > const & minImages = wavelist.minImages();
      RField<1,CPT> const & ksq = wavelist.kSq();
      DArray<bool> const & implicit = wavelist.implicitInverse();

      // Compute minimum images (and ksq) locally and compare
      IntVec<1> temp, vec;
      MeshIterator<1> iter;
      double val;
      int rank;
      bool flag;
      iter.setDimensions(kMeshDims1);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh1.dimensions(), cell1);
         flag = implicit[rank];
         val = cell1.ksq(vec);
         
         TEST_ASSERT(vec == minImages[rank]);
         TEST_ASSERT(abs(val - ksq[rank]) < tolerance_);
         TEST_ASSERT((flag == FFT<1,CPT>::hasImplicitInverse(temp, meshDims1)));
      }

   }

   void testComputeMinimumImages2D()
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(mesh2.dimension(0) == 32);
      TEST_ASSERT(mesh2.dimension(1) == 48);
      TEST_ASSERT(kMeshDims2[0] == 32);
      TEST_ASSERT(kMeshDims2[1] == 25);
      TEST_ASSERT(kSize2 == 32*25);

      // Set up wavelist object
      WaveList<2,CPT> wavelist;
      wavelist.allocate(mesh2, cell2);

      // compute minimum images (and ksq)
      wavelist.computeMinimumImages(); 
      DArray< IntVec<2> > const & minImages = wavelist.minImages();
      RField<2,CPT> ksq = wavelist.kSq();
      TEST_ASSERT(minImages.capacity() == kSize2);
      TEST_ASSERT(ksq.capacity() == kSize2);

      // compute minimum images (and ksq) locally and compare
      IntVec<2> temp, vec;
      MeshIterator<2> iter;
      double val;
      iter.setDimensions(kMeshDims2);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh2.dimensions(), cell2);
         val = cell2.ksq(vec);
         nyquistReference(temp, mesh2.dimensions(), cell2, vec, val, nullptr);
         
         TEST_ASSERT(vec == minImages[iter.rank()]);
         TEST_ASSERT(abs(val - ksq[iter.rank()]) < tolerance_);
      }
   }

   void testComputeMinimumImages3D()
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(mesh3.dimension(0) == 24);
      TEST_ASSERT(mesh3.dimension(1) == 64);
      TEST_ASSERT(mesh3.dimension(2) == 21);
      TEST_ASSERT(kMeshDims3[0] == 24);
      TEST_ASSERT(kMeshDims3[1] == 64);
      TEST_ASSERT(kMeshDims3[2] == 11);
      TEST_ASSERT(kSize3 == 24*64*11);

      // Set up wavelist object
      WaveList<3,CPT> wavelist;
      wavelist.allocate(mesh3, cell3);

      // Compute minimum images (and ksq)
      wavelist.computeMinimumImages(); 
      DArray< IntVec<3> > const & minImages = wavelist.minImages();
      RField<3,CPT> const & ksq = wavelist.kSq();

      // Compute minimum images (and ksq) locally and compare
      IntVec<3> temp, vec;
      MeshIterator<3> iter;
      double val;
      iter.setDimensions(kMeshDims3);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh3.dimensions(), cell3);
         val = cell3.ksq(vec);

         TEST_ASSERT(vec == minImages[iter.rank()]);
         TEST_ASSERT(abs(val - ksq[iter.rank()]) < tolerance_);
      }
   }

   void testComputeKSq1D()
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(mesh1.dimension(0) == 32);
      TEST_ASSERT(kMeshDims1[0] == 17);
      TEST_ASSERT(kSize1 == 17);

      // set up wavelist object
      WaveList<1,CPT> wavelist;
      wavelist.allocate(mesh1, cell1);

      // Compute kSq two different ways
      wavelist.computeMinimumImages(); // calculates kSq
      RField<1,CPT> const & ksq = wavelist.kSq();
      wavelist.clearUnitCellData(); // resets kSq but not min images
      wavelist.computeKSq(); // recalculates kSq 
      RField<1,CPT> const & ksq2 = wavelist.kSq();

      // Compute kSq locally and compare
      IntVec<1> temp, vec;
      double val;
      MeshIterator<1> iter;
      iter.setDimensions(kMeshDims1);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh1.dimensions(), cell1);
         val = cell1.ksq(vec);

         TEST_ASSERT(abs(val - ksq[iter.rank()]) < tolerance_);
         TEST_ASSERT(abs(val - ksq2[iter.rank()]) < tolerance_);
      }
   }

   void testComputeKSq2D()
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(mesh2.dimension(0) == 32);
      TEST_ASSERT(mesh2.dimension(1) == 48);
      TEST_ASSERT(kMeshDims2[0] == 32);
      TEST_ASSERT(kMeshDims2[1] == 25);
      TEST_ASSERT(kSize2 == 32*25);

      // Set up unit cell with no flexible angles
      // (if there are flexible angles, computeKSq never gets to run,
      //  because computeMinimumImages is always called instead.)
      UnitCell<2> cell;
      std::ifstream in;
      openInputFile("in/Rectangular", in);
      in >> cell;

      // Set up wavelist object
      WaveList<2,CPT> wavelist;
      wavelist.allocate(mesh2, cell);

      // Compute kSq two different ways
      wavelist.computeMinimumImages(); // calculates kSq
      RField<2,CPT> const & ksq = wavelist.kSq();
      wavelist.clearUnitCellData(); // resets kSq but not min images
      wavelist.computeKSq(); // recalculates kSq using a different kernel
      RField<2,CPT> const & ksq2 = wavelist.kSq();

      // Compute kSq in wavelist and test
      IntVec<2> pos, vec;
      double val;
      MeshIterator<2> iter;
      iter.setDimensions(kMeshDims2);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         pos = iter.position();
         vec = shiftToMinimum(pos, mesh2.dimensions(), cell);
         val = cell.ksq(vec);

         TEST_ASSERT(abs(val - ksq[iter.rank()]) < tolerance_);
         TEST_ASSERT(abs(val - ksq2[iter.rank()]) < tolerance_);
      }
   }

   void testComputeKSq3D()
   {
      printMethod(TEST_FUNC);

      // Set up unit cell with no flexible angles
      // (if there are flexible angles, computeKSq never gets to run,
      //  because computeMinimumImages is always called instead.)
      UnitCell<3> cell;
      std::ifstream in;
      openInputFile("in/Hexagonal", in);
      in >> cell;

      // Set up wavelist object
      WaveList<3,CPT> wavelist;
      wavelist.allocate(mesh3, cell);

      // Compute kSq two different ways
      wavelist.computeMinimumImages(); // calculates kSq
      RField<3,CPT> const & ksq = wavelist.kSq();
      wavelist.clearUnitCellData(); // resets kSq but not min images
      wavelist.computeKSq(); // recalculates kSq using a different kernel
      RField<3,CPT> const & ksq2 = wavelist.kSq();

      // Compute kSq locally and compare
      IntVec<3> temp, vec;
      double val;
      MeshIterator<3> iter;
      iter.setDimensions(kMeshDims3);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh3.dimensions(), cell);
         val = cell.ksq(vec);

         TEST_ASSERT(abs(val - ksq[iter.rank()]) < tolerance_);
         TEST_ASSERT(abs(val - ksq2[iter.rank()]) < tolerance_);
      }
   }

   void testComputedKSq1D()
   {
      printMethod(TEST_FUNC);

      // Set up wavelist object
      WaveList<1,CPT> wavelist;
      wavelist.allocate(mesh1, cell1);

      // Compute dKSq
      wavelist.computedKSq();
      RField<1,CPT> dksq = wavelist.dKSq(0);

      // Compute dKSq locally and compare
      IntVec<1> pos, vec;
      double val;
      MeshIterator<1> iter;
      iter.setDimensions(kMeshDims1);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         pos = iter.position();
         vec = shiftToMinimum(pos, mesh1.dimensions(), cell1);
         val = cell1.dksq(vec, 0);
         if (pos[0] != 0) {
            if (mesh1.dimension(0) - pos[0] > mesh1.dimension(0)/2) {
               val *= 2.0;
            }
         }
         TEST_ASSERT(abs(val - dksq[iter.rank()]) < tolerance_);
      }
   }

   void testComputedKSq2D()
   {
      printMethod(TEST_FUNC);

      WaveList<2,CPT> wavelist;
      wavelist.allocate(mesh2, cell2);
      wavelist.computedKSq();

      // compute dKSq locally and compare
      IntVec<2> pos, vec;
      double val, ksqRef, dksqRef[6];
      MeshIterator<2> iter;
      iter.setDimensions(kMeshDims2);
      for (int n = 0; n < cell2.nParameter() ; ++n) {
         RField<2,CPT> const & dksq = wavelist.dKSq(n);
         for (iter.begin(); !iter.atEnd(); ++iter) {
            pos = iter.position();
            vec = shiftToMinimum(pos, mesh2.dimensions(), cell2);
            val = cell2.dksq(vec, n);
            if (nyquistReference(pos, mesh2.dimensions(), cell2, vec,
                                 ksqRef, dksqRef)) {
               val = dksqRef[n];
            }
            if (pos[1] != 0) {
               if (mesh2.dimension(1) - pos[1] > mesh2.dimension(1)/2){
                  val *= 2.0;
               }
            }
            TEST_ASSERT(abs(val - dksq[iter.rank()]) < tolerance_);
         }
      }
   }

   void testComputedKSq3D()
   {
      printMethod(TEST_FUNC);

      WaveList<3,CPT> wavelist;
      wavelist.allocate(mesh3, cell3);
      wavelist.computedKSq();

      // compute dKSq locally and compare
      IntVec<3> pos, vec;
      double val;
      MeshIterator<3> iter;
      iter.setDimensions(kMeshDims3);
      for (int n = 0; n < cell3.nParameter() ; ++n) {
         RField<3,CPT> const & dksq = wavelist.dKSq(n);
         for (iter.begin(); !iter.atEnd(); ++iter) {
            pos = iter.position();
            vec = shiftToMinimum(pos, mesh3.dimensions(), cell3);
            val = cell3.dksq(vec, n);
            if (pos[2] != 0) {
               if (mesh3.dimension(2) - pos[2] > mesh3.dimension(2)/2){
                  val *= 2.0;
               }
            }
            TEST_ASSERT(abs(val - dksq[iter.rank()]) < tolerance_);
         }
      }
   }

   void testSortWaves3D()
   {
      printMethod(TEST_FUNC);

      // Set up unit cell with no flexible angles
      // (if there are flexible angles, computeKSq never gets to run,
      //  because computeMinimumImages is always called instead.)
      UnitCell<3> cell;
      std::ifstream in;
      openInputFile("in/Hexagonal", in);
      in >> cell;

      // Set up wavelist object
      WaveList<3,CPT> wavelist;
      wavelist.allocate(mesh3, cell);

      // Compute kSq 
      wavelist.computeMinimumImages(); 
      wavelist.computeKSq(); 
      RField<3,CPT> const & ksq = wavelist.kSq();
      int kSize = wavelist.kSize();
      TEST_ASSERT(kSize > 0);

      wavelist.sortWaves();
      int nBunch = wavelist.nBunch(); 
     
      // Test sortedIds 
      DArray<int> const & sortedIds = wavelist.sortedIds();
      for (int i = 1; i < kSize; ++i) {
         TEST_ASSERT(ksq[sortedIds[i]] >= ksq[sortedIds[i-1]]);
      }

      // Test sortedBunches
      double epsilon = 1.0E-8;
      double value = 0.0;
      double beginValue = 0.0;
      int begin, end, size, ib, iw;
      GArray< Pair<int> > const & sortedBunches = wavelist.sortedBunches();
      for (ib = 0; ib < nBunch; ++ib) {
         begin = sortedBunches[ib][0];
         end = sortedBunches[ib][1];
         size = end - begin;
         value = ksq[sortedIds[begin]];
         if (ib > 0) {
            TEST_ASSERT(std::abs(value - beginValue) > epsilon);
         }
         beginValue = value;
         if (size > 1) {
            for (iw = begin + 1; iw < end; ++iw) {
               value = ksq[sortedIds[iw]];
               TEST_ASSERT(std::abs(value - beginValue) < 2.0*epsilon);
            }
         }
      }

      // Test bunchIds
      TEST_ASSERT(nBunch > 0);
      DArray<int> const & bunchIds = wavelist.bunchIds();
      TEST_ASSERT(bunchIds.capacity() == kSize);
      DArray<double> oldVal(nBunch);
      DArray<bool> found(nBunch);
      double newVal;
      for (ib = 0; ib < nBunch; ++ib) {
         found[ib] = false;
      }
      for (iw = 0; iw < kSize; ++iw) {
         newVal = ksq[iw];
         ib = bunchIds[iw];
         if (!found[ib]) {
            oldVal[ib] = newVal;
            found[ib] = true;
         } else {
            TEST_ASSERT(std::abs(newVal - oldVal[ib]) < 2.0E-8);
         }
      }
   }

   void testComplex()
   {
      printMethod(TEST_FUNC);

      WaveList<3,CPT> wavelist(false);
      wavelist.allocate(mesh3, cell3);
      wavelist.computedKSq(); // computes min images, ksq, and dksq

      DArray< IntVec<3> > const & minImages = wavelist.minImages();
      RField<3,CPT> const & ksq = wavelist.kSq();
      DArray< RField<3,CPT> > const & dksq = wavelist.dKSq();

      // Check that array sizes are correct
      TEST_ASSERT(minImages.capacity() == mesh3.size());
      TEST_ASSERT(ksq.capacity() == mesh3.size());
      TEST_ASSERT(dksq.capacity() == cell3.nParameter());
      for (int i = 0; i < cell3.nParameter(); i++) {
         TEST_ASSERT(dksq[i].capacity() == mesh3.size());
      }

      // Compute minimum images, ksq, and dksq locally and compare
      IntVec<3> temp, vec;
      MeshIterator<3> iter;
      double val;
      iter.setDimensions(mesh3.dimensions());
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh3.dimensions(), cell3);
         val = cell3.ksq(vec);

         TEST_ASSERT(vec == minImages[iter.rank()]);
         TEST_ASSERT(abs(val - ksq[iter.rank()]) < tolerance_);

         for (int i = 0; i < cell3.nParameter(); i++) {
            val = cell3.dksq(vec, i);
            TEST_ASSERT(abs(val - dksq[i][iter.rank()]) < tolerance_);
         }
      }
   }

};

TEST_BEGIN(CppWaveListTest)
TEST_ADD(CppWaveListTest, testAllocate)
TEST_ADD(CppWaveListTest, testComputeMinimumImages1D)
TEST_ADD(CppWaveListTest, testComputeMinimumImages2D)
TEST_ADD(CppWaveListTest, testComputeMinimumImages3D)
TEST_ADD(CppWaveListTest, testComputeKSq1D)
TEST_ADD(CppWaveListTest, testComputeKSq2D)
TEST_ADD(CppWaveListTest, testComputeKSq3D)
TEST_ADD(CppWaveListTest, testComputedKSq1D)
TEST_ADD(CppWaveListTest, testComputedKSq2D)
TEST_ADD(CppWaveListTest, testComputedKSq3D)
TEST_ADD(CppWaveListTest, testSortWaves3D)
TEST_ADD(CppWaveListTest, testComplex)
TEST_END(CppWaveListTest)

#endif
