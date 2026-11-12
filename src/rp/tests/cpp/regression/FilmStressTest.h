#ifndef RPC_FILM_STRESS_TEST_H
#define RPC_FILM_STRESS_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/scft/ScftThermo.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/Mask.h>
#include <rp/environment/MixAndMatchEnv.h>
#include <rp/environment/FieldGenerator.h>

#include <util/containers/FSArray.h>
#include <util/format/Dbl.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests for the stress of thin films.
*/
class FilmStressTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Read a parameter file to initialize a System.
   */
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
   * Compare the stress of a 1D film (bead model, rigid cell, athermal
   * walls) with the derivative of fHelmholtz with respect to the film
   * thickness L, computed by central finite differences.
   *
   * The stress is the sum of the Mixture stress and the stress terms
   * of the field generators. The mask term contains d(phiTot)/dL.
   */
   void testStressFiniteDifference1D()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmStressFiniteDifference1D.log");

      // Converged solution for a film of thickness L
      Rp::System<1,CPT> system;
      setupSystem(system, "in/Film/fd1D_param");
      system.w().readBasis("in/Film/fd1D_w.bf");
      double L = system.domain().unitCell().parameter(0);
      int error = system.iterate();
      TEST_ASSERT(!error);

      // Stress, as used by the iterator for a flexible film thickness
      system.computeStress();
      Rp::MixAndMatchEnv const & env
            = dynamic_cast<Rp::MixAndMatchEnv const &>(system.environment());
      double stress = system.mixture().stress(0);
      stress += env.fieldGenerator1().stress(0);
      if (env.hasFieldGenerator2()) {
         stress += env.fieldGenerator2().stress(0);
      }

      // Free energies for thicknesses L + h and L - h
      double h = 1.0E-4;
      FSArray<double, 6> parameters;
      parameters.append(L + h);
      system.setUnitCell(parameters);
      error = system.iterate();
      TEST_ASSERT(!error);
      double fp = system.scft().fHelmholtz();

      system.w().readBasis("in/Film/fd1D_w.bf");
      parameters[0] = L - h;
      system.setUnitCell(parameters);
      error = system.iterate();
      TEST_ASSERT(!error);
      double fm = system.scft().fHelmholtz();

      double stressFD = (fp - fm)/(2.0*h);
      double diff = std::abs(stress - stressFD);
      Log::file() << "\nstress            = " << Dbl(stress, 20, 12)
                  << "\nstress (fin diff) = " << Dbl(stressFD, 20, 12)
                  << "\ndifference        = " << Dbl(diff, 20, 12)
                  << std::endl;
      if (verbose() > 0) {
         std::cout << "\nstress difference = " << diff << "\n";
      }
      TEST_ASSERT(diff < 1.0E-7);
   }

};

TEST_BEGIN(FilmStressTest)
TEST_ADD(FilmStressTest, testStressFiniteDifference1D)
TEST_END(FilmStressTest)

#endif
