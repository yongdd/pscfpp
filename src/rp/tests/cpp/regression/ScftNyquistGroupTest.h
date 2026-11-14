#ifndef RPC_SCFT_NYQUIST_GROUP_TEST_H
#define RPC_SCFT_NYQUIST_GROUP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <pscf/field/cpp/WaveList.h>
#include <pscf/backend/cpp/HostRecvArray.h>
#include <pscf/crystal/Basis.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/DArray.h>
#include <util/containers/FSArray.h>
#include <util/math/Constants.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Space groups with 3- or 6-fold rotations in a hexagonal setting, used
* with a variable-angle lattice (rhombic, triclinic) at 120 degrees, on
* an even mesh.
*
* kSq of wavevectors with Nyquist indices was averaged over the signs of
* these indices for every variable-angle lattice. These rotations map
* Nyquist wavevectors onto others, so that waves of one star had
* different kSq (8x8 mesh: (4,1) and (-1,-3)), and the concentration
* computed from a symmetric field was not symmetric. The average is now
* used only near right angles (see ksqNyquist), so that kSq is the
* minimum-image |k|^2 (within rounding) at the metrics of these groups.
* A group whose operations are signed permutations (R_-3_m:R on a
* rhombohedral lattice) is also checked.
*/
class ScftNyquistGroupTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   template <int D>
   void setupSystem(Rp::System<D,CPT>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

   /*
   * Read a system, set the unit cell, and check kSq and the symmetry
   * of the concentrations computed from a rough symmetric w field.
   *
   * If average is true, the reference kSq of a wavevector with a
   * Nyquist index is ksqNyquist; otherwise the reference kSq of every
   * wavevector is the minimum-image |k|^2.
   */
   template <int D>
   void checkSystem(std::string fname, FSArray<double, 6> const & p,
                    bool average)
   {
      Rp::System<D,CPT> system;
      setupSystem<D>(system, fname);
      system.setUnitCell(p);
      UnitCell<D> const & cell = system.domain().unitCell();
      IntVec<D> const & d = system.domain().mesh().dimensions();

      // kSq of every wavevector: minimum image or average
      WaveList<D,CPT>& waveList = system.waveList();
      waveList.computeKSq();
      HostRecvArray<double,CPT> kSq(waveList.kSq());
      Mesh<D> kMesh(waveList.kMeshDimensions());
      MeshIterator<D> iter(waveList.kMeshDimensions());
      int nNyquist = 0;
      double maxError = 0.0;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         IntVec<D> v = iter.position();
         double ref;
         if (average && hasNyquistIndex(v, d)) {
            ref = ksqNyquist(shiftToMinimumNyquist(v, d, cell), d, cell);
         } else {
            ref = cell.ksq(shiftToMinimum(v, d, cell));
         }
         if (hasNyquistIndex(v, d)) ++nNyquist;
         double error = std::abs(kSq[iter.rank()] - ref)/(1.0 + ref);
         maxError = std::max(maxError, error);
      }
      TEST_ASSERT(nNyquist > 0);
      Log::file() << fname << "  max kSq error " << maxError << std::endl;
      TEST_ASSERT(maxError < 1.0E-10);

      // Equal kSq within each star
      Basis<D> const & basis = system.domain().basis();
      double maxStarDiff = 0.0;
      for (int is = 0; is < basis.nStar(); ++is) {
         typename Basis<D>::Star const & star = basis.star(is);
         double kSq0 = -1.0;
         for (int iw = star.beginId; iw < star.endId; ++iw) {
            typename Basis<D>::Wave const & wave = basis.wave(iw);
            if (wave.implicit) continue;
            double value = kSq[kMesh.rank(wave.indicesStd)];
            if (kSq0 < 0.0) {
               kSq0 = value;
            } else {
               double diff = std::abs(value - kSq0)/(1.0 + kSq0);
               maxStarDiff = std::max(maxStarDiff, diff);
            }
         }
      }
      Log::file() << fname << "  max kSq difference in a star "
                  << maxStarDiff << std::endl;
      TEST_ASSERT(maxStarDiff < 1.0E-10);

      // Rough symmetric w fields: the c fields must be symmetric
      int nMonomer = system.mixture().nMonomer();
      int nBasis = basis.nBasis();
      DArray< DArray<double> > w;
      w.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w[i].allocate(nBasis);
         w[i][0] = 0.0;
         for (int j = 1; j < nBasis; ++j) {
            w[i][j] = 2.0*std::sin(1.3*j + 0.7*i + 0.2);
         }
      }
      system.w().setBasis(w);
      system.compute();
      for (int i = 0; i < nMonomer; ++i) {
         bool isSymmetric = system.domain().fieldIo().hasSymmetry(
                                 system.c().rgrid(i), 1.0E-8, false);
         Log::file() << fname << "  c field " << i << " symmetric "
                     << isSymmetric << std::endl;
         TEST_ASSERT(isSymmetric);
      }
   }

   /*
   * Rhombic lattice at 120 degrees with p_6_m_m, 8x8 mesh.
   */
   void testRhombicHexagonal()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftNyquistGroup2.log");
      FSArray<double, 6> p;
      p.append(3.0);
      p.append(2.0*Constants::Pi/3.0);
      checkSystem<2>("in/Scft/nyq_rhombic_p6mm", p, false);
   }

   /*
   * Triclinic lattice with a = b, gamma = 120 degrees and c normal to
   * a and b, with P_6%m_m_m, 8x8x6 mesh.
   */
   void testTriclinicHexagonal()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftNyquistGroup3.log");
      FSArray<double, 6> p;
      p.append(3.0);
      p.append(3.0);
      p.append(2.5);
      p.append(0.0);
      p.append(0.0);
      p.append(2.0*Constants::Pi/3.0);
      checkSystem<3>("in/Scft/nyq_triclinic_p6mmm", p, false);
   }

   /*
   * Rhombohedral lattice with R_-3_m:R (signed permutations): kSq of
   * Nyquist wavevectors is ksqNyquist, and c is symmetric.
   */
   void testRhombohedral()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftNyquistGroupR.log");
      FSArray<double, 6> p;
      p.append(3.0);
      p.append(1.2);
      checkSystem<3>("in/Scft/nyq_rhombohedral_r3m", p, true);
   }

};

TEST_BEGIN(ScftNyquistGroupTest)
TEST_ADD(ScftNyquistGroupTest, testRhombicHexagonal)
TEST_ADD(ScftNyquistGroupTest, testTriclinicHexagonal)
TEST_ADD(ScftNyquistGroupTest, testRhombohedral)
TEST_END(ScftNyquistGroupTest)

#endif
