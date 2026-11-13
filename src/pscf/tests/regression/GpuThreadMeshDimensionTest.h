#ifndef PSCF_GPU_THREAD_MESH_DIMENSION_TEST_H
#define PSCF_GPU_THREAD_MESH_DIMENSION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/backend/cuda/ThreadMesh.h>
#include <pscf/backend/cuda/HostSendArray.h>
#include <pscf/backend/cuda/cudaErrorCheck.h>

using namespace Util;
using namespace Pscf;

namespace GpuThreadMeshDimensionTestDetail {

   // A non-idempotent, bounded per-node operation exposes duplicate launches.
   template <int D>
   __global__ void countVisits(int* visits, dim3 mesh)
   {
      int x = blockIdx.x*blockDim.x + threadIdx.x;
      int y = blockIdx.y*blockDim.y + threadIdx.y;
      int z = blockIdx.z*blockDim.z + threadIdx.z;
      int id = x;
      if (D > 1) id += mesh.x*y;
      if (D > 2) id += mesh.x*mesh.y*z;
      if (x < mesh.x && y < mesh.y && z < mesh.z) atomicAdd(visits+id, 1);
   }

}

/* Dropping dimensions must replace the cached axes and launch each node once. */
class GpuThreadMeshDimensionTest : public UnitTest
{
public:

   template <int D>
   void checkVisits(IntVec<D> const & dimensions, bool invert)
   {
      int n = 1;
      for (int i = 0; i < D; ++i) n *= dimensions[i];
      BackendArray<int,CUT> counts(n);
      cudaErrorCheck(cudaMemset(counts.cArray(), 0, n*sizeof(int)));
      dim3 mesh = ThreadMesh::meshDims();
      dim3 block = ThreadMesh::blockDims();
      dim3 grid = ThreadMesh::gridDims();
      GpuThreadMeshDimensionTestDetail::countVisits<D><<<grid,block>>>
         (counts.cArray(), mesh);
      cudaErrorCheck(cudaGetLastError());
      HostSendArray<int,CUT> host(n);
      host = counts;
      for (int i = 0; i < n; ++i) TEST_ASSERT(host[i] == 1);
      TEST_ASSERT(mesh.x == unsigned(dimensions[invert ? D-1 : 0]));
      if (D > 1) TEST_ASSERT(mesh.y == unsigned(dimensions[invert ? D-2 : 1]));
      if (D > 2) TEST_ASSERT(mesh.z == unsigned(dimensions[invert ? D-3 : 2]));
      if (D < 3) {
         TEST_ASSERT(mesh.z == 1);
         TEST_ASSERT(block.z == 1);
         TEST_ASSERT(grid.z == 1);
      }
      if (D < 2) {
         TEST_ASSERT(mesh.y == 1);
         TEST_ASSERT(block.y == 1);
         TEST_ASSERT(grid.y == 1);
      }
   }

   void checkTransitions(bool invert)
   {
      IntVec<3> d3;
      d3[0] = invert ? 4 : 16;
      d3[1] = 8;
      d3[2] = invert ? 16 : 4;
      IntVec<2> d2;
      d2[0] = invert ? 8 : 16;
      d2[1] = invert ? 16 : 8;
      IntVec<1> d1(16);
      // The remaining x/y dimensions are identical across these calls.
      ThreadMesh::setConfig(d3, invert, 128);
      checkVisits(d3, invert);
      ThreadMesh::setConfig(d2, invert, 128);
      checkVisits(d2, invert);
      // Default block size permits cache reuse, so this catches a dropped y.
      ThreadMesh::setConfig(d1, invert);
      checkVisits(d1, invert);
      ThreadMesh::setConfig(d3, invert, 128);
      checkVisits(d3, invert);
      ThreadMesh::setConfig(d1, invert);
      checkVisits(d1, invert);
      // Repeat an identical request as a valid cache-hit control.
      dim3 block = ThreadMesh::blockDims();
      dim3 grid = ThreadMesh::gridDims();
      ThreadMesh::setConfig(d1, invert);
      checkVisits(d1, invert);
      TEST_ASSERT(ThreadMesh::blockDims().x == block.x);
      TEST_ASSERT(ThreadMesh::gridDims().x == grid.x);
   }

   void testInvertedDimensions()
   {
      printMethod(TEST_FUNC);
      checkTransitions(true);
   }

   void testNaturalDimensions()
   {
      printMethod(TEST_FUNC);
      checkTransitions(false);
   }
};

TEST_BEGIN(GpuThreadMeshDimensionTest)
TEST_ADD(GpuThreadMeshDimensionTest, testInvertedDimensions)
TEST_ADD(GpuThreadMeshDimensionTest, testNaturalDimensions)
TEST_END(GpuThreadMeshDimensionTest)

#endif
