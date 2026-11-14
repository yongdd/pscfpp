/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "HomogeneousComparison.h"

#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <r1/solvers/Block.h>
#include <r1/solvers/Solvent.h>

#include <pscf/interaction/Interaction.h>
#include <pscf/floryHuggins/FhInteraction.h>
#include <pscf/floryHuggins/FhMixture.h>

#include <util/format/Str.h>
#include <util/format/Int.h>
#include <cmath>
#include <util/format/Dbl.h>
#include <util/misc/Log.h>

#include <string>

namespace Pscf {
namespace R1
{

   using namespace Util;

   /*
   * Default constructor.
   */
   HomogeneousComparison::HomogeneousComparison()
    : SystemAccess(),
      fhMixturePtr_(nullptr),
      fhInteractionPtr_(nullptr),
      hasData_(false)
   {}

   /*
   * Constructor.
   */
   HomogeneousComparison::HomogeneousComparison(System& system)
    : SystemAccess(system)
   {
      fhMixturePtr_ = new FhMixture();
      fhInteractionPtr_ = new FhInteraction();
   }

   /*
   * Destructor.
   */
   HomogeneousComparison::~HomogeneousComparison()
   {
      if (fhMixturePtr_) {
         delete fhMixturePtr_;
      }
      if (fhInteractionPtr_) {
         delete fhInteractionPtr_;
      }
   }

   /*
   * Compute properties of a homogeneous reference system.
   *
   * mode == 0 : Composition equals spatial average composition
   * mode == 1:  Chemical potential equal to that of system,
   *             composition guess given at last grid point.
   * mode == 2:  Chemical potential equal to that of system,
   *             composition guess given at first grid point.
   */
   void HomogeneousComparison::compute(int mode)
   {
    
      UTIL_CHECK(fhMixturePtr_);
      UTIL_CHECK(fhInteractionPtr_);
      fhMixture().initialize(mixture()); 
      fhInteraction() = interaction(); 

      int np = mixture().nPolymer();
      int ns = mixture().nSolvent();
      int nt = np + ns;
      UTIL_CHECK(fhMixture().nMolecule() == nt);
      if (!p_.isAllocated()) {
         p_.allocate(nt);
      }
      UTIL_CHECK(p_.capacity() == nt);

      if (mode == 0) {

         // Spatial average volume fractions of the SCFT solution. That
         // of an open species is computed from the fields, so in a
         // mixed or grand-canonical ensemble the total is the spatial
         // average of the total concentration, which is 1 only within
         // the SCFT error tolerance. The homogeneous reference system
         // is incompressible: rescale the fractions to a sum of 1.
         // This keeps absent species (phi = 0) at 0, keeps the relative
         // precision of dilute species, and does not depend on the
         // order of the species. (Closed fractions add up to 1 already.)
         double sum = 0.0;
         for (int i = 0; i < np; ++i) {
            p_[i] = mixture().polymer(i).phi();
            sum += p_[i];
         }
         for (int i = 0; i < ns; ++i) {
            p_[i + np] = mixture().solvent(i).phi();
            sum += p_[i + np];
         }
         Log::file() << "Sum of species volume fractions = "
                     << Dbl(sum, 20, 13) << std::endl;
         if (!(sum > 0.0) || !std::isfinite(sum)) {
            UTIL_THROW("Species volume fractions of the SCFT solution "
                       "do not have a positive finite sum");
         }
         // A converged solution gives |sum - 1| up to about the iterator
         // tolerance, below 1e-6 for the tolerances of the examples
         // (1e-7 to 1e-9). Above 1e-6 (e.g., a loosely converged
         // solution, or COMPUTE without ITERATE), warn that the
         // reference composition differs from the averages. This is
         // not an error: like the other outputs, the comparison does
         // not require a converged solution.
         if (std::abs(sum - 1.0) > 1.0E-6) {
            Log::file() << "Warning: Species volume fractions do not "
                        << "add up to 1 (SCFT solution not converged?).\n"
                        << "The homogeneous reference system (mode 0) "
                        << "uses them rescaled to a sum of 1."
                        << std::endl;
         }
         for (int i = 0; i < nt; ++i) {
            p_[i] /= sum;
         }
         fhMixture().setComposition(p_);
         double xi = 0.0;
         fhMixture().computeMu( fhInteraction(), xi );

      } else 
      if (mode == 1 || mode == 2) {

         // Chemical potentials
         if (!m_.isAllocated()) {
            m_.allocate(np+ns);
         }
         UTIL_CHECK(m_.capacity() == fhMixture().nMolecule());
         for (int i = 0; i < np; ++i) {
            m_[i] = mixture().polymer(i).mu(); 
         }
         for (int i = 0; i < ns; ++i) {
            m_[i + np] = mixture().solvent(i).mu(); 
         }

         // Grid index from which we obtain guess of composition
         int ix; 
         if (mode == 1) {
            ix = domain().nx() - 1;
         } else 
         if (mode == 2) {
            ix = 0;
         }

         // Compute initial guess for array p_ of species volume fractions
         // Initial guess taken from boundary point with index ix
         for (int i = 0; i < np; ++i) {
            // For polymers, use sum of phi values for all blocks 
            p_[i] = 0.0;
            int nb = mixture().polymer(i).nBlock();
            for (int j = 0; j < nb; ++j) {
               Block const & block = mixture().polymer(i).block(j);
               p_[i] += block.cField()[ix];
            }
         }
         for (int i = 0; i < ns; ++i) {
            p_[i + np] = mixture().solvent(i).cField()[ix];
         }

         // A converged SCFT field need not give a valid composition at
         // the boundary point exactly (a fraction can exceed 1 slightly,
         // or be zero or slightly negative): replace values that are not
         // positive by a tiny positive one, clip values above 1, and
         // normalize the guess to a sum of 1. Positive values are kept,
         // however small: computePhi takes Newton steps in ln(phi). A
         // species with mu = -infinity (phi = 0) stays exactly at 0.
         double sum = 0.0;
         for (int i = 0; i < np + ns; ++i) {
            if (std::isinf(m_[i]) && m_[i] < 0.0) {
               p_[i] = 0.0;
            } else {
               if (!(p_[i] > 0.0)) p_[i] = 1.0E-300;
               if (p_[i] > 1.0) p_[i] = 1.0;
            }
            sum += p_[i];
         }
         for (int i = 0; i < np + ns; ++i) {
            p_[i] /= sum;
         }

         double xi = 0.0;
         fhMixture().computePhi(fhInteraction(), m_, p_, xi);

      } else {
         UTIL_THROW("Unknown mode in computeHomogeneous");
      }

      // Compute Helmholtz free energy and pressure
      fhMixture().computeFreeEnergy(fhInteraction());

      hasData_ = true;
   }

   /*
   * Output properties of homogeneous system and free energy difference.
   *
   * Mode 0:      Outputs fHomo (Helmholtz) and difference df
   * Mode 1 or 2: Outputs Helhmoltz, pressure and difference dV
   */
   void HomogeneousComparison::output(int mode, std::ostream& out) const
   {
      int np = mixture().nPolymer();
      int ns = mixture().nSolvent();
      int nt = np + ns;

      // Preconditions
      UTIL_CHECK(hasData_);
      UTIL_CHECK(fhMixture().nMolecule() == nt);

      if (mode == 0) {

         // Output free energies
         double fHomo = fhMixture().fHelmholtz();
         double df = system().fHelmholtz() - fHomo;
         out << std::endl;
         out << "f (homo)    = " << Dbl(fHomo, 18, 11) 
                                 << "   [per monomer volume]" << std::endl;
         out << "delta f     = " << Dbl(df, 18, 11)  
                                 << "   [per monomer volume]" << std::endl;

         // Output species-specific properties
         out << std::endl;
         out << "Species:" << std::endl;
         out << "    i"
             << "      mu(homo)      "
             << "      phi(homo)     "
             << std::endl;
         for (int i = 0; i < nt; ++i) {
            out << Int(i,5)
                << "  " << Dbl(fhMixture().mu(i), 18, 11)
                << "  " << Dbl(fhMixture().phi(i), 18, 11) 
                << std::endl;
         }
         out << std::endl;

      } else
      if (mode == 1 || mode == 2) {

         // Output free energies
         double fHomo = fhMixture().fHelmholtz();
         double pHomo = fhMixture().pressure();
         double pEx   = system().pressure() - pHomo;
         double fEx   = system().fHelmholtz() - fHomo;
         double V     = domain().volume()/mixture().vMonomer();
         double PhiExTot = -1.0*pEx*V;
         double fExTot   = fEx*V; 
         out << std::endl;
         out << "f (homo)   = " << Dbl(fHomo, 18, 11)     
                                << "   [per monomer volume]" << std::endl;
         out << "p (homo)   = " << Dbl(pHomo, 18, 11)     
                                << "   [per monomer volume]" << std::endl;
         out << "delta f    = " << Dbl(fEx, 18, 11)  
                                << "   [per monomer volume]" << std::endl;
         out << "delta p    = " << Dbl(pEx, 18, 11)  
                                << "   [per monomer volume]" << std::endl;
         out << "Phi (ex)   = " << Dbl(PhiExTot, 18, 11) 
                                << "   [total]" << std::endl;
         out << "F (ex)     = " << Dbl(fExTot, 18, 11)  
                                << "   [total]" << std::endl;
         out << "V(tot)/v   = " << Dbl(V, 18, 11) << std::endl;

         // Output species specific properties
         double dV; 
         out << std::endl;
         out << "Species:" << std::endl;
         out << "    i"
             << "      mu            "
             << "      phi(homo)     "
             << "      deltaV        " 
             << std::endl;

         // Polymers
         for (int i = 0; i < np; ++i) {
            dV = mixture().polymer(i).phi() - fhMixture().phi(i);
            dV *= V;
            out << Int(i,5)
                << "  " << Dbl(fhMixture().mu(i), 18, 11)
                << "  " << Dbl(fhMixture().phi(i), 18, 11) 
                << "  " << Dbl(dV, 18, 11)
                << std::endl;
         }

         // Solvents
         int j;
         for (int i = 0; i < ns; ++i) {
            j = i + np;
            dV = mixture().solvent(i).phi() - fhMixture().phi(j);
            dV *= V;
            out << Int(j,5)
                << "  " << Dbl(fhMixture().mu(j), 18, 11)
                << "  " << Dbl(fhMixture().phi(j), 18, 11) 
                << "  " << Dbl(dV, 18, 11)
                << std::endl;
         }

         out << std::endl;
      }

   }

} // namespace R1
} // namespace Pscf
