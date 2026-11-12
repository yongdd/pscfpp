#ifndef RPG_FTS_REL_NORM_TEST_H
#define RPG_FTS_REL_NORM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/compressor/AmCompressor.h>
#include <rp/fts/compressor/LrAmCompressor.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/field/cuda/RField.h>
#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/Reduce.h>

#include <util/misc/FileMaster.h>
#include <util/tests/LogFileUnitTest.h>

#include <gsl/gsl_errno.h>

#include <fstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

namespace {

   // Number of calls of the GSL error handler
   int ftsRelNormGslErrors = 0;

   // GSL error handler that counts errors instead of aborting
   void ftsRelNormGslHandler(const char*, const char*, int, int)
   {  ++ftsRelNormGslErrors; }

}

/*
* Regression tests (pscfpp fixes): errorType relNorm in the AM
* compressors.
*/
class FtsRelNormTest : public LogFileUnitTest
{

public:

   void setUp()
   {  
      setVerbose(0); 
      ftsRelNormGslErrors = 0;
      oldHandler_ = gsl_set_error_handler(&ftsRelNormGslHandler);
   }

   void tearDown()
   {
      gsl_set_error_handler(oldHandler_);
      LogFileUnitTest::tearDown();
   }

   template <typename CompressorT>
   void init(Rp::System<3,CUT>& system, CompressorT& compressor,
             std::string compressorFileName)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile("in/Fts/param_system_small", in);
      system.readParam(in);
      in.close();
      openInputFile(compressorFileName, in);
      compressor.readParam(in);
      in.close();
   }

   /*
   * Set homogeneous w fields, w_i = values[i].
   */
   void setHomogeneous(Rp::System<3,CUT>& system, double const * values)
   {
      int nMonomer = system.mixture().nMonomer();
      DArray< RField<3,CUT> > w;
      w.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w[i].allocate(system.domain().mesh().dimensions());
         VecOp::eqS(w[i], values[i]);
      }
      system.w().setRGrid(w);
   }

   /*
   * Do the w fields have only finite values?
   */
   bool isFinite(Rp::System<3,CUT>& system)
   {
      int nMonomer = system.mixture().nMonomer();
      for (int i = 0; i < nMonomer; ++i) {
         if (!std::isfinite(Reduce::sumSq(system.w().rgrid(i)))) {
            return false;
         }
      }
      return true;
   }

   /*
   * A homogeneous state of the diblock melt is incompressible. The
   * compressor must accept it at the first iteration (the change of
   * the pressure field vanishes there, and the residual vanishes for
   * w = 0: an error |R|/|X| gives NaN or inf).
   */
   template <typename CompressorT>
   void checkHomogeneous(std::string compressorFileName)
   {
      double values[3][2] = {{0.0, 0.0}, {7.2, 4.8}, {1.0, 1.0}};
      for (int k = 0; k < 3; ++k) {
         Rp::System<3,CUT> system;
         CompressorT compressor(system);
         init(system, compressor, compressorFileName);
         system.w().readRGrid("in/Fts/w_small.rf"); // sets unit cell
         setHomogeneous(system, values[k]);
         int error = compressor.compress();
         TEST_ASSERT(error == 0);
         TEST_ASSERT(isFinite(system));
      }
      TEST_ASSERT(ftsRelNormGslErrors == 0);
   }

   void testAmHomogeneous()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsRelNormAmHomogeneous.log");
      checkHomogeneous< Rp::AmCompressor<3,CUT> >(
                                 "in/Fts/param_AmCompressor_relNorm");
   }

   void testLrAmHomogeneous()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsRelNormLrAmHomogeneous.log");
      checkHomogeneous< Rp::LrAmCompressor<3,CUT> >(
                                 "in/Fts/param_LrAmCompressor_relNorm");
   }

   /*
   * A second compress() of a compressed state must accept it at the
   * first iteration, without further iterations.
   */
   template <typename CompressorT>
   void checkCompressTwice(std::string compressorFileName)
   {
      Rp::System<3,CUT> system;
      CompressorT compressor(system);
      init(system, compressor, compressorFileName);
      system.w().readRGrid("in/Fts/w_small.rf");
      TEST_ASSERT(compressor.compress() == 0);
      int mdeCounter = compressor.mdeCounter();
      TEST_ASSERT(compressor.compress() == 0);
      TEST_ASSERT(compressor.mdeCounter() == mdeCounter + 1);
      TEST_ASSERT(isFinite(system));
      TEST_ASSERT(ftsRelNormGslErrors == 0);
   }

   void testAmCompressTwice()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsRelNormAmCompressTwice.log");
      checkCompressTwice< Rp::AmCompressor<3,CUT> >(
                                 "in/Fts/param_AmCompressor_relNorm");
   }

   void testLrAmCompressTwice()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FtsRelNormLrAmCompressTwice.log");
      checkCompressTwice< Rp::LrAmCompressor<3,CUT> >(
                                 "in/Fts/param_LrAmCompressor_relNorm");
   }

private:

   gsl_error_handler_t* oldHandler_;

};

TEST_BEGIN(FtsRelNormTest)
TEST_ADD(FtsRelNormTest, testAmHomogeneous)
TEST_ADD(FtsRelNormTest, testLrAmHomogeneous)
TEST_ADD(FtsRelNormTest, testAmCompressTwice)
TEST_ADD(FtsRelNormTest, testLrAmCompressTwice)
TEST_END(FtsRelNormTest)

#endif
