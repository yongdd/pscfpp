#ifndef RPG_FTS_COMPRESSOR_TEST_H
#define RPG_FTS_COMPRESSOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/compressor/LrCompressor.h>
#include <rp/fts/compressor/LrAmCompressor.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/field/cuda/RField.h>
#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/Reduce.h>
#include <util/containers/FSArray.h>

#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests (pscfpp fixes): LR compressors: intramolecular correlation of the current state.
*/
class FtsCompressorTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   template <int D>
   void initSystem(Rp::System<D,CUT>& system, std::string filename)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile(filename, in);
      system.readParam(in);
      in.close();
   }

   template <typename CompressorT>
   void initCompressor(CompressorT& compressor, std::string filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      compressor.readParam(in);
      in.close();
   }

   /*
   * Root-mean-square deviation of the total concentration from 1.
   */
   template <int D>
   double incompressibilityError(Rp::System<D,CUT> const & system)
   {
      RField<D,CUT> error;
      error.allocate(system.domain().mesh().dimensions());
      VecOp::eqS(error, -1.0);
      int nMonomer = system.mixture().nMonomer();
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::addEqV(error, system.c().rgrid(i));
      }
      int meshSize = system.domain().mesh().size();
      return std::sqrt(Reduce::sumSq(error)/double(meshSize));
   }

   /*
   * Read a system, a compressor and w fields, then compress.
   *
   * Return the value returned by compress() (0 if converged).
   */
   template <int D, typename CompressorT>
   int readAndCompress(Rp::System<D,CUT>& system,
                       CompressorT& compressor,
                       std::string systemFileName,
                       std::string compressorFileName,
                       std::string wFileName)
   {
      initSystem(system, systemFileName);
      initCompressor(compressor, compressorFileName);
      system.w().readRGrid(wFileName);
      return compressor.compress();
   }

   /*
   * LR compressors in a grand-canonical blend (two open homopolymers),
   * compressed right after reading w, before any MDE solution.
   */
   void testLrCompressorOpen()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsLrCompressorOpen.log");

      Rp::System<1,CUT> system1;
      Rp::LrCompressor<1,CUT> lrCompressor(system1);
      int error = readAndCompress(system1, lrCompressor,
                                  "in/Fts/param_system_gc_1D",
                                  "in/Fts/param_LrCompressor_small",
                                  "in/Fts/w_gc_1D.rf");
      TEST_ASSERT(error == 0);
      TEST_ASSERT(incompressibilityError(system1) < 1.0E-6);

      Rp::System<1,CUT> system2;
      Rp::LrAmCompressor<1,CUT> lrAmCompressor(system2);
      error = readAndCompress(system2, lrAmCompressor,
                              "in/Fts/param_system_gc_1D",
                              "in/Fts/param_LrAmCompressor_small",
                              "in/Fts/w_gc_1D.rf");
      TEST_ASSERT(error == 0);
      TEST_ASSERT(incompressibilityError(system2) < 1.0E-6);
   }

   /*
   * LrCompressor used again after the unit cell has changed.
   */
   void testLrCompressorCellChange()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsLrCompressorCellChange.log");

      Rp::System<3,CUT> system;
      Rp::LrCompressor<3,CUT> compressor(system);
      int error = readAndCompress(system, compressor,
                                  "in/Fts/param_system_small",
                                  "in/Fts/param_LrCompressor_small",
                                  "in/Fts/w_small.rf");
      TEST_ASSERT(error == 0);

      // Triple the cell size (3.0 to 9.0) and compress again
      FSArray<double, 6> parameters;
      parameters.append(9.0);
      system.setUnitCell(parameters);
      error = compressor.compress();
      TEST_ASSERT(error == 0);
      TEST_ASSERT(incompressibilityError(system) < 1.0E-6);
   }

   /*
   * LrCompressor in a mixed ensemble: closed diblock (phi 0.8) and a
   * short open homopolymer. The k = 0 response comes from the open
   * species only.
   */
   void testLrCompressorMixed()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsLrCompressorMixed.log");

      Rp::System<1,CUT> system;
      Rp::LrCompressor<1,CUT> compressor(system);
      int error = readAndCompress(system, compressor,
                                  "in/Fts/param_system_mixed_1D",
                                  "in/Fts/param_LrCompressor_small",
                                  "in/Fts/w_mixed_1D.rf");
      TEST_ASSERT(error == 0);
      TEST_ASSERT(incompressibilityError(system) < 1.0E-6);
   }

   /*
   * LR compressors in a mixed ensemble (bead model: open diblock,
   * closed solvent) from a converged w shifted by +0.03. The volume
   * fraction of the open diblock at the first MDE solution (about 0.38)
   * differs from its final value (0.70), so the correlation must follow
   * it.
   */
   void testLrCompressorOpenPhi()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsLrCompressorOpenPhi.log");

      Rp::System<3,CUT> system1;
      Rp::LrCompressor<3,CUT> lrCompressor(system1);
      int error = readAndCompress(system1, lrCompressor,
                                  "in/Fts/param_system_bead_open",
                                  "in/Fts/param_LrCompressor_bead",
                                  "in/Fts/w_bead_open.rf");
      TEST_ASSERT(error == 0);
      TEST_ASSERT(incompressibilityError(system1) < 1.0E-6);

      Rp::System<3,CUT> system2;
      Rp::LrAmCompressor<3,CUT> lrAmCompressor(system2);
      error = readAndCompress(system2, lrAmCompressor,
                              "in/Fts/param_system_bead_open",
                              "in/Fts/param_LrAmCompressor_bead",
                              "in/Fts/w_bead_open.rf");
      TEST_ASSERT(error == 0);
      TEST_ASSERT(incompressibilityError(system2) < 1.0E-6);
   }

};

TEST_BEGIN(FtsCompressorTest)
TEST_ADD(FtsCompressorTest, testLrCompressorOpen)
TEST_ADD(FtsCompressorTest, testLrCompressorCellChange)
TEST_ADD(FtsCompressorTest, testLrCompressorMixed)
TEST_ADD(FtsCompressorTest, testLrCompressorOpenPhi)
TEST_END(FtsCompressorTest)

#endif
