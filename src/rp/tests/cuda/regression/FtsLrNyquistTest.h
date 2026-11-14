#ifndef RPG_FTS_LR_NYQUIST_TEST_H
#define RPG_FTS_LR_NYQUIST_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/compressor/IntraCorrelation.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/field/cuda/WaveList.h>
#include <pscf/field/cuda/FFT.h>
#include <pscf/field/cuda/RField.h>
#include <pscf/field/cuda/RFieldDft.h>
#include <pscf/backend/cuda/HostRecvArray.h>
#include <pscf/backend/cuda/HostSendArray.h>
#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/complex.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/DArray.h>
#include <util/containers/FSArray.h>
#include <util/math/Constants.h>
#include <util/format/Dbl.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Intramolecular correlation of the LR compressors at wavevectors with
* Nyquist indices, in a skewed (triclinic) cell on an even mesh.
*
* The LR estimate of the Jacobian, vMonomer*omega(k), must use the kSq
* of the WaveList, which the MDE solver uses; IntraCorrelation takes it
* from the WaveList. For the bead model, the response of the total
* concentration to a small change of all w fields at one wavevector is
* exactly -vMonomer*omega(k) (to linear order).
*/
class FtsLrNyquistTest : public LogFileUnitTest
{

public:

   typedef CUT BT;

   void setUp()
   {  setVerbose(0); }

   void setupSystem(Rp::System<3,BT>& system)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile("in/Fts/nyq_lr_param", in);
      system.readParam(in);
      in.close();
      FSArray<double, 6> p;
      p.append(5.0);
      p.append(4.0);
      p.append(6.0);
      p.append(0.4);
      p.append(0.5);
      p.append(1.2);
      system.setUnitCell(p);
   }

   /*
   * Response -dc(k)/dw(k) of the total concentration to a change of
   * all w fields by eps*cos(k.r), at the wavevector with indices v.
   */
   double response(Rp::System<3,BT>& system, IntVec<3> const & v)
   {
      double eps = 1.0E-4;
      Mesh<3> const & mesh = system.domain().mesh();
      IntVec<3> const & d = mesh.dimensions();
      int nMonomer = system.mixture().nMonomer();

      // Perturbed w fields (host, then device)
      HostSendArray<cudaReal,CUT> wh(mesh.size());
      MeshIterator<3> iter(d);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         double phase = 0.0;
         for (int j = 0; j < 3; ++j) {
            phase += double(v[j]*iter.position(j))/double(d[j]);
         }
         wh[iter.rank()] = eps*std::cos(2.0*Constants::Pi*phase);
      }
      DArray< RField<3,BT> > w;
      w.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w[i].allocate(d);
         w[i] = wh;
      }
      system.w().setRGrid(w);
      system.compute();

      // Fourier coefficients of w and of the total concentration
      RField<3,BT> cTot;
      cTot.allocate(d);
      VecOp::addVV(cTot, system.c().rgrid(0), system.c().rgrid(1));
      RFieldDft<3,BT> wk, ck;
      wk.allocate(d);
      ck.allocate(d);
      system.domain().fft().forwardTransform(w[0], wk);
      system.domain().fft().forwardTransform(cTot, ck);
      HostRecvArray<typename BT::Complex,BT> wkh(wk);
      HostRecvArray<typename BT::Complex,BT> ckh(ck);
      Mesh<3> kMesh(system.waveList().kMeshDimensions());
      int rank = kMesh.rank(v);
      return -real(ckh[rank])/real(wkh[rank]);
   }

   /*
   * vMonomer*omega of IntraCorrelation equals the response of the
   * solver at wavevectors with and without Nyquist indices, including
   * (0,3,5), whose kSq the WaveList averages over the signs of its
   * Nyquist indices (by 3.6% above the minimum-image |k|^2).
   */
   void testNyquistOmega()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsLrNyquist.log");

      Rp::System<3,BT> system;
      setupSystem(system);
      IntVec<3> const & d = system.domain().mesh().dimensions();

      // Omega from IntraCorrelation
      WaveList<3,BT>& waveList = system.waveList();
      IntVec<3> const & kd = waveList.kMeshDimensions();
      Rp::IntraCorrelation<3,BT> intra(system);
      RField<3,BT> omega;
      omega.allocate(kd);
      intra.computeOmegaTotal(omega);
      HostRecvArray<cudaReal,BT> omegaH(omega);

      // Response of the solver
      double vMonomer = system.mixture().vMonomer();
      Mesh<3> kMesh(kd);
      int vs[6][3] = {{1, 1, 1}, {4, 1, 2}, {1, 3, 2}, {2, 1, 5},
                      {4, 3, 1}, {0, 3, 5}};
      double maxResponseError = 0.0;
      for (int n = 0; n < 6; ++n) {
         IntVec<3> v;
         for (int j = 0; j < 3; ++j) v[j] = vs[n][j];
         if (n == 5) {
            TEST_ASSERT(usesNyquistAverage(v, d, system.domain().unitCell()));
         }
         double r = response(system, v);
         int rank = kMesh.rank(v);
         double lr = vMonomer*omegaH[rank];
         double error = std::abs(r - lr)/std::abs(r);
         Log::file() << "wave " << v
                     << "  Nyquist " << hasNyquistIndex(v, d)
                     << "  response " << Dbl(r, 20, 12)
                     << "  vMonomer*omega " << Dbl(lr, 20, 12)
                     << "  rel. error " << error << std::endl;
         maxResponseError = std::max(maxResponseError, error);
      }
      TEST_ASSERT(maxResponseError < 1.0E-6);
   }

};

TEST_BEGIN(FtsLrNyquistTest)
TEST_ADD(FtsLrNyquistTest, testNyquistOmega)
TEST_END(FtsLrNyquistTest)

#endif
