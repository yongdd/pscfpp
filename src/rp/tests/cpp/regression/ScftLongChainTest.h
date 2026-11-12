#ifndef RPC_SCFT_LONG_CHAIN_TEST_H
#define RPC_SCFT_LONG_CHAIN_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/scft/ScftThermo.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Solvent.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/chem/PolymerModel.h>
#include <pscf/backend/Reduce.h>

#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Homogeneous solutions of a homopolymer and a solvent for long chains
* and large solvent molecules, whose partition functions are outside
* the range of a double, and contour steps too small for an int count.
*/
class ScftLongChainTest : public LogFileUnitTest
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
   * Spatial average of an r-grid field.
   */
   template <int D>
   double average(RField<D,CPT> const & field)
   {  return Reduce::sum(field)/double(field.capacity()); }

   /*
   * Solve the SCFT equations for a homopolymer (monomer 0) and a
   * solvent (monomer 1), and compare with Flory-Huggins theory.
   *
   * The chemical potentials are compared in the gauge of the converged
   * fields: mu = ln phi + N <w> for each species of size N.
   */
   void checkFloryHuggins(std::string paramFile, std::string wFile)
   {
      Rp::System<1,CPT> system;
      setupSystem(system, paramFile);
      system.w().readBasis(wFile);
      int error = system.iterate();
      TEST_ASSERT(!error);
      system.scft().compute();

      Rp::Polymer<1,CPT> const & polymer = system.mixture().polymer(0);
      Rp::Solvent<1,CPT> const & solvent = system.mixture().solvent(0);
      double nP;
      if (PolymerModel::isThread()) {
         nP = polymer.length();
      } else {
         nP = (double) polymer.nBead();
      }
      double nS = solvent.size();
      double phiP = polymer.phi();
      double phiS = solvent.phi();
      double chi = system.interaction().chi(0,1);

      double f = phiP*(std::log(phiP) - 1.0)/nP
               + phiS*(std::log(phiS) - 1.0)/nS + chi*phiP*phiS;
      double wP = average(system.w().rgrid(0));
      double wS = average(system.w().rgrid(1));
      double muP = std::log(phiP) + nP*wP;
      double muS = std::log(phiS) + nS*wS;
      double p = -f + phiP*muP/nP + phiS*muS/nS;

      double fScft = system.scft().fHelmholtz();
      double pScft = system.scft().pressure();
      Log::file() << "f = " << fScft << " (FH " << f << ")"
                  << "  p = " << pScft << " (FH " << p << ")"
                  << "\nmu polymer = " << polymer.mu()
                  << " (FH " << muP << ")"
                  << "  mu solvent = " << solvent.mu()
                  << " (FH " << muS << ")" << std::endl;

      TEST_ASSERT(std::isfinite(fScft));
      TEST_ASSERT(std::abs(fScft - f) < 1.0E-8);
      TEST_ASSERT(std::abs(pScft - p) < 1.0E-8);
      TEST_ASSERT(std::abs(polymer.mu() - muP) < 1.0E-6);
      TEST_ASSERT(std::abs(solvent.mu() - muS) < 1.0E-6);
      TEST_ASSERT(std::abs(phiP + phiS - 1.0) < 1.0E-8);
      RField<1,CPT> const & c = system.c().rgrid(0);
      TEST_ASSERT(std::abs(Reduce::max(c) - phiP) < 1.0E-8);
      TEST_ASSERT(std::abs(Reduce::min(c) - phiP) < 1.0E-8);
   }

   /*
   * Thread model, homopolymer length 1800 (ln q about -810).
   */
   void testThread()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftLongChainTestThread.log");
      checkFloryHuggins("in/Scft/long_thread", "in/Scft/long_w.bf");
   }

   /*
   * Bead model, homopolymer of 1700 beads.
   */
   void testBead()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftLongChainTestBead.log");
      checkFloryHuggins("in/Scft/long_bead", "in/Scft/long_w.bf");
   }

   /*
   * Open ensemble for both species, homopolymer length 1790, with
   * mu = ln(0.1) + 0.4*1790 (exp(mu) above the double range). The
   * initial field differs from the solution by 1E-3, because phi of
   * the polymer depends on its field as exp(-1790 w).
   */
   void testOpen()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftLongChainTestOpen.log");
      checkFloryHuggins("in/Scft/long_open", "in/Scft/long_open_w.bf");
   }

   /*
   * Solvent of size 1800 in a monomeric homopolymer. The initial field
   * differs from the solution by 1E-4, because the concentration of
   * the solvent depends on its field as exp(-1800 w).
   */
   void testLargeSolvent()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/ScftLongChainTestLargeSolvent.log");
      checkFloryHuggins("in/Scft/large_solvent",
                        "in/Scft/large_solvent_w.bf");
   }

};

TEST_BEGIN(ScftLongChainTest)
TEST_ADD(ScftLongChainTest, testThread)
TEST_ADD(ScftLongChainTest, testBead)
TEST_ADD(ScftLongChainTest, testOpen)
TEST_ADD(ScftLongChainTest, testLargeSolvent)
TEST_END(ScftLongChainTest)

#endif
