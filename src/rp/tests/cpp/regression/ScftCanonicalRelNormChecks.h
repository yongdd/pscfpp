#ifndef RP_SCFT_CANONICAL_REL_NORM_CHECKS_H
#define RP_SCFT_CANONICAL_REL_NORM_CHECKS_H

#include <util/tests/LogFileUnitTest.h>
#include <util/param/BracketPolicy.h>
#include <util/misc/Exception.h>
#include <util/containers/DArray.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <pscf/field/RField.h>
#include <pscf/backend/VecOp.h>
#include <pscf/chem/PolymerModel.h>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

namespace Pscf {
namespace Rp {

/*
* The relNorm error of the AM SCFT iterators for a closed system
* (canonical ensemble) must not be divided by a norm that contains a
* field common to all monomer types (a pressure field), nor external
* fields: the residual sees neither, so a large such field made the
* error arbitrarily small at a state that is far from a solution.
* The bead solver keeps the concentrations finite in arbitrarily strong
* nonuniform fields, so that a diverged pressure field was reported as
* converged.
*/
template <class T>
class ScftCanonicalRelNormChecks : public LogFileUnitTest
{

   Util::BracketPolicy::Type bracketPolicy_;

   /*
   * Read a parameter file and set the unit cell.
   */
   void read(System<1,T>& system, std::string const & text, double cell)
   {
      std::istringstream in(text);
      system.readParam(in);
      Util::FSArray<double,6> parameters;
      parameters.append(cell);
      system.setUnitCell(parameters);
   }

   /*
   * Allocate r-grid fields and set them from host values.
   */
   void setValues(Util::DArray<RField<1,T> >& fields,
                  std::vector<double> const & a,
                  std::vector<double> const & b)
   {
      const int n = a.size();
      fields.allocate(2);
      Util::DArray<typename T::Real> host;
      host.allocate(n);
      for (int m = 0; m < 2; ++m) {
         fields[m].allocate(IntVec<1>(n));
         for (int i = 0; i < n; ++i) host[i] = m ? b[i] : a[i];
         VecOp::eqV(fields[m], host);
      }
   }

   /*
   * Iterate; return -1 if an Exception was thrown.
   */
   int iterate(System<1,T>& system)
   {
      int error = 0;
      try {
         error = system.iterate();
      } catch (Util::Exception& e) {
         error = -1;
      }
      return error;
   }

   /*
   * Maximum of |c_A + c_B - 1| over the grid.
   */
   double incompressibilityError(System<1,T>& system)
   {
      const int n = system.domain().mesh().size();
      Util::DArray<typename T::Real> a, b;
      a.allocate(n);
      b.allocate(n);
      VecOp::eqV(a, system.c().rgrid(0));
      VecOp::eqV(b, system.c().rgrid(1));
      double error = 0.0;
      for (int i = 0; i < n; ++i) {
         double value = std::abs(a[i] + b[i] - 1.0);
         if (!(value <= error)) error = value;
      }
      return error;
   }

   /*
   * Diblock A_nA B_nB, chi_AB = chi, phi = 1, mesh 32.
   */
   std::string diblock(std::string const & model,
                       std::string const & blocks, double chi,
                       std::string const & iterator,
                       std::string const & options)
   {
      std::ostringstream text;
      text << "System{ polymerModel " << model
           << " Mixture{ nMonomer 2 monomers[ 1 1 ] nPolymer 1 "
           << "Polymer{ type linear nBlock 2 blocks[ " << blocks
           << " ] phi 1.0 } "
           << (model == "thread" ? "ds 0.01 " : "") << "} "
           << "Interaction{ chi( 1 0 " << chi << " ) } "
           << "Domain{ mesh 32 lattice lamellar groupName P_-1 } "
           << iterator << "{ " << options << " } }";
      return text.str();
   }

public:

   void setUp()
   {
      bracketPolicy_ = Util::BracketPolicy::get();
      Util::BracketPolicy::set(Util::BracketPolicy::Optional);
   }

   void tearDown()
   {
      LogFileUnitTest::tearDown();
      PolymerModel::setModel(PolymerModel::Thread);
      Util::BracketPolicy::set(bracketPolicy_);
   }

   /*
   * A disordered bead diblock (chi 1.5, solution w = 0.75) with a
   * common field xi cos(2 pi x/L), xi = 1e14, added to both w fields,
   * which confines the chains to one plane. The residual is of order
   * 1, but the error was |R|/|w|, about 1e-13. The state may not be
   * accepted (one iteration allowed). (The thread solver gives NaN for
   * such fields.)
   */
   void checkPressureField(std::string const & iterator)
   {
      System<1,T> system;
      read(system, diblock("bead", "0 8 1 8", 1.5, iterator,
                           "epsilon 1e-9 maxItr 1 isFlexible 0"), 4.0);
      const double pi = std::acos(-1.0);
      std::vector<double> w(32);
      for (int i = 0; i < 32; ++i) {
         w[i] = 0.75 + 1.0E14*std::cos(2.0*pi*i/32.0);
      }
      Util::DArray<RField<1,T> > fields;
      setValues(fields, w, w);
      system.w().setRGrid(fields, true);
      TEST_ASSERT(iterate(system) == 1);
   }

   /*
   * A/B homopolymer blend (phi 0.5 each, chi 0.1, mesh 8) with an
   * external field H on monomer A at grid points 3-5, starting from
   * w = h + 0.05. The error was |R|/|w|, about 1e-10 for the bead model
   * with H = 1e9 and 1e-4 for the thread model with H = 1e3, while the
   * total concentration differs from 1 by up to 0.5. Neither state may
   * be accepted (one iteration allowed).
   */
   void checkExternalField(std::string const & iterator)
   {
      std::string models[] = {"bead", "thread"};
      for (int k = 0; k < 2; ++k) {
         const bool bead = (k == 0);
         const double field = bead ? 1.0E9 : 1.0E3;
         std::ostringstream text;
         text << "System{ polymerModel " << models[k]
              << " Mixture{ nMonomer 2 monomers[ 1 1 ] nPolymer 2 "
              << "Polymer{ type linear nBlock 1 blocks[ 0 "
              << (bead ? "10" : "1.0") << " ] phi 0.5 } "
              << "Polymer{ type linear nBlock 1 blocks[ 1 "
              << (bead ? "10" : "1.0") << " ] phi 0.5 } "
              << (bead ? "" : "ds 0.01 ") << "} "
              << "Interaction{ chi( 1 0 0.1 ) } "
              << "Domain{ mesh 8 lattice lamellar groupName P_-1 } "
              << iterator << "{ epsilon " << (bead ? "1e-8" : "1e-3")
              << " maxItr 1 isFlexible 0 } }";
         System<1,T> system;
         read(system, text.str(), 4.0);
         std::vector<double> h(8, 0.0), zero(8, 0.0), wa(8), wb(8);
         for (int i = 3; i <= 5; ++i) h[i] = field;
         for (int i = 0; i < 8; ++i) {
            wa[i] = h[i] + 0.05;
            wb[i] = 0.05;
         }
         Util::DArray<RField<1,T> > hFields, wFields;
         setValues(hFields, h, zero);
         setValues(wFields, wa, wb);
         system.h().setRGrid(hFields, true);
         system.w().setRGrid(wFields, true);
         TEST_ASSERT(iterate(system) == 1);
      }
   }

   /*
   * The bead diblock above, from wA = 1.5(0.5 - 0.4 cos 2 pi x/L) and
   * wB = 1.5(0.5 + 0.4 cos 2 pi x/L), with default AM parameters. The
   * Anderson iteration diverges along the pressure field, and was
   * reported as converged after about 80 iterations with sum c between
   * 0 and 16 and f of order -1e20. The solve must now either fail or
   * end at an incompressible state.
   */
   void checkDivergedIteration(std::string const & iterator)
   {
      System<1,T> system;
      read(system, diblock("bead", "0 8 1 8", 1.5, iterator,
                           "epsilon 1e-9 maxItr 200 isFlexible 0"), 4.0);
      const double pi = std::acos(-1.0);
      std::vector<double> wa(32), wb(32);
      for (int i = 0; i < 32; ++i) {
         const double c = 0.4*std::cos(2.0*pi*i/32.0);
         wa[i] = 1.5*(0.5 - c);
         wb[i] = 1.5*(0.5 + c);
      }
      Util::DArray<RField<1,T> > fields;
      setValues(fields, wa, wb);
      system.w().setRGrid(fields, true);
      const int error = iterate(system);
      TEST_ASSERT(error == 0 || error == 1);
      if (error == 0) {
         TEST_ASSERT(incompressibilityError(system) < 1.0E-6);
      }
   }

};

}
}
#endif
