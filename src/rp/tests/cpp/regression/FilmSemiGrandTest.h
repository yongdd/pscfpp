#ifndef RPC_FILM_SEMI_GRAND_TEST_H
#define RPC_FILM_SEMI_GRAND_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/scft/ScftThermo.h>
#include <rp/scft/iterator/Iterator.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Solvent.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/environment/MixAndMatchEnv.h>
#include <rp/environment/FilmFieldGenMaskBase.h>

#include <pscf/chem/Ensemble.h>

#include <util/containers/DArray.h>
#include <util/containers/FSArray.h>
#include <util/format/Dbl.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests for the optimisation of the thickness of thin films
* (dynamic mask with fBulk) in mixed ensembles.
*
* The film contains a closed AB diblock (8 + 12 beads, phi = 0.8), an
* A homopolymer (10 beads) and a B solvent (size 2). The input files in
* in/Film/ were made as follows (bead model, mesh 48, chi = 1 per bead,
* walls chi = 0.5 for A, interfaceThickness 0.2, excludedThickness 0.4):
*  - bulk lamellar phase with all species closed (phi 0.8, 0.1, 0.1),
*    flexible cell: fBulk = fHelmholtz and the chemical potentials;
*  - sg_mixed_param: homopolymer and solvent open at these chemical
*    potentials, fBulk = bulk g = fHelmholtz - sum over open species
*    of phi*mu/N;
*  - sg_closed_param: all species closed, fBulk = bulk fHelmholtz;
*  - sg_*_w.bf: rigid solutions at L = 2.2, the starting points.
*/
class FilmSemiGrandTest : public LogFileUnitTest
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
   * Semi-grand potential per monomer reference volume (bead model):
   * g = fHelmholtz - sum over open species of phi*mu/N, with N the
   * number of beads of a polymer or the size of a solvent.
   */
   template <int D>
   double semiGrand(Rp::System<D,CPT> const & system)
   {
      Rp::Mixture<D,CPT> const & mixture = system.mixture();
      double g = system.scft().fHelmholtz();
      for (int i = 0; i < mixture.nPolymer(); ++i) {
         if (mixture.polymer(i).ensemble() == Ensemble::Open) {
            g -= mixture.polymer(i).phi()*mixture.polymer(i).mu()
                 / double(mixture.polymer(i).nBead());
         }
      }
      for (int i = 0; i < mixture.nSolvent(); ++i) {
         if (mixture.solvent(i).ensemble() == Ensemble::Open) {
            g -= mixture.solvent(i).phi()*mixture.solvent(i).mu()
                 / mixture.solvent(i).size();
         }
      }
      return g;
   }

   /*
   * Relax the thickness L of a 1D film, then return the derivative of
   * (L - t)(g - fBulk) with respect to L at the relaxed thickness, from
   * central finite differences of rigid solutions (t is the excluded
   * thickness). On return, L is the relaxed thickness.
   */
   double relaxAndDifferentiate(std::string paramFile, std::string wFile,
                                double& L)
   {
      Rp::System<1,CPT> system;
      setupSystem(system, paramFile);
      system.w().readBasis(wFile);
      int error = system.iterate();
      TEST_ASSERT(!error);
      L = system.domain().unitCell().parameter(0);

      Rp::MixAndMatchEnv const & env
            = dynamic_cast<Rp::MixAndMatchEnv const &>(system.environment());
      Rp::FilmFieldGenMaskBase<1> const & mask
            = dynamic_cast< Rp::FilmFieldGenMaskBase<1> const &>
                                                (env.fieldGenerator1());
      TEST_ASSERT(mask.hasFBulk());
      double t = mask.excludedThickness();
      double gBulk = mask.fBulk();

      // Rigid solutions for thicknesses L + h and L - h, in a second
      // system whose unit cell is made rigid before its first solve
      DArray< DArray<double> > w0;
      w0 = system.w().basis();
      Rp::System<1,CPT> rigidSystem;
      setupSystem(rigidSystem, paramFile);
      FSArray<bool, 6> rigid;
      rigid.append(false);
      rigidSystem.iterator().setFlexibleParams(rigid);
      double h = 1.0E-3;
      FSArray<double, 6> parameters;
      parameters.append(L + h);
      rigidSystem.setUnitCell(parameters);
      rigidSystem.w().setBasis(w0);
      error = rigidSystem.iterate();
      TEST_ASSERT(!error);
      double gp = semiGrand(rigidSystem);

      parameters[0] = L - h;
      rigidSystem.setUnitCell(parameters);
      rigidSystem.w().setBasis(w0);
      error = rigidSystem.iterate();
      TEST_ASSERT(!error);
      double gm = semiGrand(rigidSystem);

      double dG = ((L + h - t)*(gp - gBulk) - (L - h - t)*(gm - gBulk))
                  / (2.0*h);
      Log::file() << "\nrelaxed L = " << Dbl(L, 20, 12)
                  << "\nd[(L-t)(g-gBulk)]/dL = " << Dbl(dG, 20, 12)
                  << std::endl;
      return dG;
   }

   /*
   * Mixed ensemble (closed diblock, open homopolymer and solvent): at
   * the relaxed thickness, (L - t)(g - gBulk) is stationary. For a
   * closed mixture (g = fHelmholtz), the relaxed thickness is unchanged.
   */
   void testRelaxThickness1D()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmSemiGrandRelaxThickness1D.log");

      double L;
      double dG = relaxAndDifferentiate("in/Film/sg_mixed_param",
                                        "in/Film/sg_mixed_w.bf", L);
      if (verbose() > 0) {
         std::cout << "\nMixed:  L = " << L << ",  dG/dL = " << dG;
      }
      TEST_ASSERT(std::abs(dG) < 1.0E-6);

      dG = relaxAndDifferentiate("in/Film/sg_closed_param",
                                 "in/Film/sg_closed_w.bf", L);
      if (verbose() > 0) {
         std::cout << "\nClosed: L = " << L << ",  dG/dL = " << dG;
      }
      TEST_ASSERT(std::abs(dG) < 1.0E-6);
      TEST_ASSERT(std::abs(L - 2.2296630622) < 1.0E-7);
   }

};

TEST_BEGIN(FilmSemiGrandTest)
TEST_ADD(FilmSemiGrandTest, testRelaxThickness1D)
TEST_END(FilmSemiGrandTest)

#endif
