#ifndef PSCF_CORE_AM_ITERATOR_TEST_H
#define PSCF_CORE_AM_ITERATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/iterator/AmIteratorTmpl.h>
#include <pscf/iterator/NanException.h>
#include <pscf/backend/cpp/VecOp.h>
#include <pscf/backend/cpp/Reduce.h>
#include <util/param/ParamComposite.h>
#include <util/containers/DArray.h>
#include <util/global.h>

#include <pscf/iterator/AmIteratorTmpl.tpp>

#include <gsl/gsl_errno.h>

#include <string>
#include <vector>
#include <cmath>

using namespace Util;
using namespace Pscf;

namespace {

   // Number of calls of the GSL error handler
   int coreAmGslErrors = 0;

   // GSL error handler that counts errors instead of aborting
   void coreAmGslHandler(const char*, const char*, int, int)
   {  ++coreAmGslErrors; }

}

/*
* Base class (template argument Iterator) of the toy AM iterator.
*/
class CoreAmToyBase : public ParamComposite
{};

/*
* Toy problem for AmIteratorTmpl, with state and residual vectors of
* n elements. The first evaluations return prescribed residuals; later
* ones return R = target - X, which simple mixing solves in one step.
*/
class CoreAmToy : public AmIteratorTmpl< CoreAmToyBase, DArray<double> >
{

public:

   CoreAmToy(int n, std::string errorType)
    : n_(n),
      nEval_(0)
   {
      epsilon_ = 1.0E-10;
      maxItr_ = 100;
      maxHist_ = 10;
      verbose_ = 0;
      errorType_ = errorType;
      UTIL_CHECK(isValidErrorType());
      x_.allocate(n);
      r_.allocate(n);
      target_.allocate(n);
      for (int i = 0; i < n; ++i) {
         x_[i] = 0.0;
         r_[i] = 0.0;
         target_[i] = 0.0;
      }
   }

   // State vector
   DArray<double> x_;

   // Target state (solution after the prescribed residuals)
   DArray<double> target_;

   // Prescribed residuals of the first evaluations
   std::vector< std::vector<double> > prescribed_;

   // Number of evaluations
   int nEval() const
   {  return nEval_; }

   // Largest absolute difference between the state and the target
   double maxError() const
   {
      double error = 0.0;
      for (int i = 0; i < n_; ++i) {
         double diff = std::abs(x_[i] - target_[i]);
         if (!(diff <= error)) error = diff; // also NaN
      }
      return error;
   }

private:

   DArray<double> r_;
   int n_;
   int nEval_;

   int nElements()
   {  return n_; }

   bool hasInitialGuess()
   {  return true; }

   void getCurrent(DArray<double>& curr)
   {
      for (int i = 0; i < n_; ++i) curr[i] = x_[i];
   }

   void evaluate()
   {
      if (nEval_ < (int)prescribed_.size()) {
         for (int i = 0; i < n_; ++i) {
            r_[i] = prescribed_[nEval_][i];
         }
      } else {
         for (int i = 0; i < n_; ++i) {
            r_[i] = target_[i] - x_[i];
         }
      }
      ++nEval_;
   }

   void getResidual(DArray<double>& resid)
   {
      for (int i = 0; i < n_; ++i) resid[i] = r_[i];
   }

   void update(DArray<double>& newGuess)
   {
      for (int i = 0; i < n_; ++i) x_[i] = newGuess[i];
   }

   void outputToLog()
   {}

};

/*
* Regression tests (pscfpp fixes): AmIteratorTmpl.
*/
class CoreAmIteratorTest : public UnitTest
{

public:

   void setUp()
   {
      coreAmGslErrors = 0;
      oldHandler_ = gsl_set_error_handler(&coreAmGslHandler);
   }

   void tearDown()
   {  gsl_set_error_handler(oldHandler_); }

   /*
   * errorType relNorm with a vanishing state vector. If the residual
   * also vanishes, the state is a solution and must be accepted at the
   * first iteration (|R|/|X| = 0/0 gave NaN). If not, the error is
   * infinite, and the iteration must proceed to the solution.
   */
   void testRelNormZeroState()
   {
      printMethod(TEST_FUNC);

      // X = 0 is the solution
      CoreAmToy toy1(4, "relNorm");
      int error = toy1.solve();
      TEST_ASSERT(error == 0);
      TEST_ASSERT(toy1.nEval() == 1);
      TEST_ASSERT(toy1.maxError() == 0.0);

      // X = 0 is not the solution
      CoreAmToy toy2(4, "relNorm");
      for (int i = 0; i < 4; ++i) {
         toy2.target_[i] = double(i + 1);
      }
      error = toy2.solve();
      TEST_ASSERT(error == 0);
      TEST_ASSERT(toy2.maxError() < 1.0E-8);

      TEST_ASSERT(coreAmGslErrors == 0);
   }

   /*
   * Exactly singular U matrix with two basis vectors. The residuals of
   * the first three evaluations are a, b, a, so that the two residual
   * basis vectors a - b and b - a cancel exactly, as after an update
   * that increases the residual by many orders of magnitude. The LU
   * solution of U c = v would call the GSL error handler (which aborts
   * by default). The iteration must continue instead, and converge.
   */
   void testSingularBasis()
   {
      printMethod(TEST_FUNC);

      double a[4] = {1.0, 2.0, 3.0, 4.0};
      double b[4] = {2.0, 0.0, 1.0, 5.0};
      CoreAmToy toy(4, "rms");
      toy.prescribed_.resize(3);
      for (int i = 0; i < 4; ++i) {
         toy.prescribed_[0].push_back(a[i]);
         toy.prescribed_[1].push_back(b[i]);
         toy.prescribed_[2].push_back(a[i]);
         toy.target_[i] = double(i + 1);
      }
      int error = toy.solve();
      TEST_ASSERT(coreAmGslErrors == 0);
      TEST_ASSERT(error == 0);
      TEST_ASSERT(toy.maxError() < 1.0E-8);
   }

   /*
   * U(0,0) = 0 with one basis vector: the residual of the second
   * evaluation equals the first one. The coefficient v/U(0,0) is not
   * defined (it gave NaN fields). The iteration must continue instead,
   * and converge.
   */
   void testZeroBasis()
   {
      printMethod(TEST_FUNC);

      double a[4] = {1.0, 2.0, 3.0, 4.0};
      CoreAmToy toy(4, "rms");
      toy.prescribed_.resize(2);
      for (int i = 0; i < 4; ++i) {
         toy.prescribed_[0].push_back(a[i]);
         toy.prescribed_[1].push_back(a[i]);
         toy.target_[i] = double(i + 1);
      }
      int error = toy.solve();
      TEST_ASSERT(coreAmGslErrors == 0);
      TEST_ASSERT(error == 0);
      TEST_ASSERT(toy.maxError() < 1.0E-8);
   }

private:

   gsl_error_handler_t* oldHandler_;

};

TEST_BEGIN(CoreAmIteratorTest)
TEST_ADD(CoreAmIteratorTest, testRelNormZeroState)
TEST_ADD(CoreAmIteratorTest, testSingularBasis)
TEST_ADD(CoreAmIteratorTest, testZeroBasis)
TEST_END(CoreAmIteratorTest)

#endif
