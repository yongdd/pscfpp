#ifndef RPG_FILM_CPU_GPU_TEST_H
#define RPG_FILM_CPU_GPU_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/Mask.h>
#include <rp/solvers/Mixture.h>
#include <rp/environment/MixAndMatchEnv.h>

#include <pscf/backend/cuda/CUT.h>
#include <pscf/backend/cuda/HostRecvArray.h>
#include <pscf/backend/cpp/CPT.h>
#include <pscf/field/cuda/RField.h>
#include <pscf/field/cpp/RField.h>

#include <util/containers/FSArray.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests that compare thin-film fields of the CPU and GPU code.
*/
class FilmCpuGpuTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Read a parameter file to initialize a System (CPU or GPU).
   */
   template <int D, class T>
   void setupSystem(Rp::System<D,T>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

   /*
   * Generate the mask and the external fields h of a 1D film with
   * asymmetric walls on the CPU and on the GPU, and return the maximum
   * absolute difference between them (L = 1.5).
   */
   double compareCpuGpuFields(std::string fname)
   {
      Rp::System<1,CPT> cpuSystem;
      setupSystem(cpuSystem, fname);

      Rp::System<1,CUT> gpuSystem;
      setupSystem(gpuSystem, fname);

      FSArray<double, 6> parameters;
      parameters.append(1.5);
      cpuSystem.setUnitCell(parameters);
      gpuSystem.setUnitCell(parameters);

      dynamic_cast<Rp::MixAndMatchEnv&>(cpuSystem.environment()).generate();
      dynamic_cast<Rp::MixAndMatchEnv&>(gpuSystem.environment()).generate();
      TEST_ASSERT(cpuSystem.h().hasData());
      TEST_ASSERT(gpuSystem.h().hasData());

      int nx = cpuSystem.domain().mesh().size();
      int nMonomer = cpuSystem.mixture().nMonomer();
      HostRecvArray<cudaReal,CUT> gpuField(nx);
      double diff;
      double maxDiff = 0.0;

      // Mask
      gpuField = gpuSystem.mask().rgrid();
      RField<1,CPT> const & cpuMask = cpuSystem.mask().rgrid();
      for (int j = 0; j < nx; ++j) {
         diff = std::abs(gpuField[j] - cpuMask[j]);
         if (diff > maxDiff) maxDiff = diff;
      }

      // External fields
      for (int i = 0; i < nMonomer; ++i) {
         gpuField = gpuSystem.h().rgrid(i);
         RField<1,CPT> const & cpuField = cpuSystem.h().rgrid(i);
         for (int j = 0; j < nx; ++j) {
            diff = std::abs(gpuField[j] - cpuField[j]);
            if (diff > maxDiff) maxDiff = diff;
         }
      }
      Log::file() << fname << ": max |GPU - CPU| = " << maxDiff
                  << std::endl;
      return maxDiff;
   }

   /*
   * The mask and external fields of the GPU must equal those of the CPU,
   * for an even and an odd number of grid points normal to the walls.
   * (Grid points with index >= N/2 are in the top half of the film.)
   */
   void testGenerateCpuGpu()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmGenerateCpuGpu.log");

      double diff = compareCpuGpuFields("in/Film/mid1D_N32");
      TEST_ASSERT(diff < 1.0E-10);
      diff = compareCpuGpuFields("in/Film/mid1D_N33");
      TEST_ASSERT(diff < 1.0E-10);
   }

};

TEST_BEGIN(FilmCpuGpuTest)
TEST_ADD(FilmCpuGpuTest, testGenerateCpuGpu)
TEST_END(FilmCpuGpuTest)

#endif
