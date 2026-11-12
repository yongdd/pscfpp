#ifndef RPC_FTS_LR_RESPONSE_TEST_H
#define RPC_FTS_LR_RESPONSE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/compressor/IntraCorrelation.h>
#include <rp/fts/compressor/LrCompressor.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/field/cpp/RField.h>
#include <util/containers/DArray.h>
#include <util/containers/FSArray.h>

#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): the intramolecular correlation used
* by the LR compressors is the linear response of the MDE solution.
*
* In a homogeneous state, a change eps*cos(k x) of all w fields changes
* the total concentration by -eps*v*omega(k)*cos(k x) to first order,
* where omega is the intramolecular correlation. The tests compute this
* response from MDE solutions and compare it with IntraCorrelation, up
* to the Nyquist wavenumber of meshes for which k*k*b*b*ds/6 >> 1 there.
*/
class FtsLrResponseTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   void initSystem(Rp::System<1,CPT>& system, std::string filename,
                   double length)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile(filename, in);
      system.readParam(in);
      in.close();

      FSArray<double, 6> parameters;
      parameters.append(length);
      system.setUnitCell(parameters);
   }

   /*
   * Set all w fields to eps*(cos(2 pi n1 j/M) + cos(2 pi n2 j/M)).
   */
   void setW(Rp::System<1,CPT>& system, double eps, int n1, int n2)
   {
      const int nMonomer = system.mixture().nMonomer();
      const int meshSize = system.domain().mesh().size();
      DArray< RField<1,CPT> > w;
      w.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w[i].allocate(system.domain().mesh().dimensions());
         for (int j = 0; j < meshSize; ++j) {
            w[i][j] = 0.0;
            if (n1 >= 0) {
               w[i][j] += eps*std::cos(2.0*M_PI*double(n1*j)/double(meshSize));
            }
            if (n2 >= 0) {
               w[i][j] += eps*std::cos(2.0*M_PI*double(n2*j)/double(meshSize));
            }
         }
      }
      system.w().setRGrid(w);
   }

   /*
   * Maximum relative deviation, over all wavenumbers of the mesh, of
   * the response -dphi/(eps v) of the MDE solution from omega.
   */
   double maxResponseDeviation(Rp::System<1,CPT>& system)
   {
      const double eps = 1.0E-4;
      const int meshSize = system.domain().mesh().size();
      const int nMonomer = system.mixture().nMonomer();
      const double vMonomer = system.mixture().vMonomer();

      IntVec<1> kMeshDimensions;
      kMeshDimensions[0] = meshSize/2 + 1;
      RField<1,CPT> omega;
      omega.allocate(kMeshDimensions);
      Rp::IntraCorrelation<1,CPT> intra(system);
      intra.computeOmegaTotal(omega);

      double maxDev = 0.0;
      for (int n = 1; n <= meshSize/2; ++n) {
         setW(system, eps, n, -1);
         system.compute();
         double dphi = 0.0;
         for (int j = 0; j < meshSize; ++j) {
            double phi = 0.0;
            for (int i = 0; i < nMonomer; ++i) {
               phi += system.c().rgrid(i)[j];
            }
            dphi += phi*std::cos(2.0*M_PI*double(n*j)/double(meshSize));
         }
         dphi *= (2*n == meshSize) ? 1.0/double(meshSize)
                                   : 2.0/double(meshSize);
         double response = -dphi/(eps*vMonomer);
         double dev = std::abs(response/omega[n] - 1.0);
         if (verbose() > 0) {
            Log::file() << "n = " << n << "  response = " << response
                        << "  omega = " << omega[n] << std::endl;
         }
         if (!(dev <= maxDev)) maxDev = dev;
      }
      return maxDev;
   }

   /*
   * Linear diblock, thread model, with unequal segment lengths.
   */
   void testThreadDiblock()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsLrResponseThreadDiblock.log");

      Rp::System<1,CPT> system;
      initSystem(system, "in/Fts/param_system_lrThread_1D", 0.5);
      double dev = maxResponseDeviation(system);
      Log::file() << "Max relative deviation = " << dev << std::endl;
      TEST_ASSERT(dev < 1.0E-5);
   }

   /*
   * Branched polymer (three arms) and a solvent, thread model.
   */
   void testThreadStar()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsLrResponseThreadStar.log");

      Rp::System<1,CPT> system;
      initSystem(system, "in/Fts/param_system_lrStar_1D", 0.5);
      double dev = maxResponseDeviation(system);
      Log::file() << "Max relative deviation = " << dev << std::endl;
      TEST_ASSERT(dev < 1.0E-5);
   }

   /*
   * Linear diblock, bead model.
   */
   void testBeadDiblock()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsLrResponseBeadDiblock.log");

      Rp::System<1,CPT> system;
      initSystem(system, "in/Fts/param_system_lrBead_1D", 2.0);
      double dev = maxResponseDeviation(system);
      Log::file() << "Max relative deviation = " << dev << std::endl;
      TEST_ASSERT(dev < 1.0E-5);
   }

   /*
   * The LrCompressor converges from a perturbed homogeneous state that
   * includes a Nyquist mode with k*k*b*b*ds/6 >> 1.
   */
   void testLrCompressorNyquist()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsLrResponseCompressor.log");

      Rp::System<1,CPT> system;
      initSystem(system, "in/Fts/param_system_lrThread_1D", 0.5);
      Rp::LrCompressor<1,CPT> compressor(system);
      std::ifstream in;
      openInputFile("in/Fts/param_LrCompressor_tight", in);
      compressor.readParam(in);
      in.close();

      const int meshSize = system.domain().mesh().size();
      setW(system, 1.0E-3, 1, meshSize/2);
      int error = compressor.compress();
      TEST_ASSERT(error == 0);
   }

};

TEST_BEGIN(FtsLrResponseTest)
TEST_ADD(FtsLrResponseTest, testThreadDiblock)
TEST_ADD(FtsLrResponseTest, testThreadStar)
TEST_ADD(FtsLrResponseTest, testBeadDiblock)
TEST_ADD(FtsLrResponseTest, testLrCompressorNyquist)
TEST_END(FtsLrResponseTest)

#endif
