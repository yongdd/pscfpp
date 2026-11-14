#ifndef RP_SCFT_RANGE_FAILURE_CHECKS_H
#define RP_SCFT_RANGE_FAILURE_CHECKS_H

#include <util/tests/LogFileUnitTest.h>
#include <util/param/BracketPolicy.h>
#include <util/misc/Exception.h>
#include <util/misc/FileMaster.h>
#include <util/containers/DArray.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <pscf/field/RField.h>
#include <pscf/backend/VecOp.h>
#include <pscf/chem/PolymerModel.h>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>

namespace Pscf {
namespace Rp {

/*
* A solver result that is not representable during an SCFT iteration
* must end the solve with failure, as a NaN error does, and not with
* an exception: ITERATE then gives a nonzero status, and SWEEP
* backtracks. A bead homopolymer of 1000 beads in a field of magnitude
* 1e304 or more has a ln Q, or an open-ensemble phi, beyond the double
* range, whatever field shift the solver uses.
*/
template <class T>
class ScftRangeFailureChecks : public LogFileUnitTest
{

   Util::BracketPolicy::Type bracketPolicy_;

   /*
   * Parameter file text: bead homopolymer A (1000 beads) and solvent
   * B (size 1), chi_AB = 2, mesh 8, rigid lamellar cell. Closed: phi
   * 0.5 for both, with a homogeneous solution at w = 0. Mixed: polymer
   * mu = ln(0.6) and solvent phi 0.4, with the homogeneous solution
   * w_A = 0, w_B = 0.4.
   */
   std::string param(std::string const & iterator, bool open,
                     std::string const & mixing = "",
                     std::string const & sweep = "")
   {
      std::ostringstream text;
      text << std::setprecision(17)
           << "System{ polymerModel bead Mixture{ nMonomer 2 "
           << "monomers[ 1 1 ] nPolymer 1 nSolvent 1 "
           << "Polymer{ type linear nBlock 1 blocks[ 0 1000 ] ";
      if (open) {
         text << "mu " << std::log(0.6);
      } else {
         text << "phi 0.5";
      }
      text << " } Solvent{ monomerId 1 size 1 phi "
           << (open ? "0.4" : "0.5") << " } } "
           << "Interaction{ chi( 0 1 2 ) } "
           << "Domain{ mesh 8 lattice lamellar groupName P_-1 } "
           << iterator << "{ epsilon 1e-10 maxItr 20 isFlexible 0 "
           << mixing << " } " << sweep << " }";
      return text.str();
   }

   /*
   * Read the parameter file text and set the unit cell.
   */
   void read(System<1,T>& system, std::string const & text)
   {
      system.fileMaster().setInputPrefix(this->filePrefix());
      system.fileMaster().setOutputPrefix(this->filePrefix());
      std::istringstream in(text);
      system.readParam(in);
      Util::FSArray<double,6> cell;
      cell.append(4.0);
      system.setUnitCell(cell);
   }

   /*
   * Set homogeneous (symmetric) w fields.
   */
   void setFields(System<1,T>& system, double wA, double wB)
   {
      Util::DArray<RField<1,T> > fields;
      fields.allocate(2);
      Util::DArray<typename T::Real> host;
      host.allocate(8);
      for (int m = 0; m < 2; ++m) {
         fields[m].allocate(IntVec<1>(8));
         for (int i = 0; i < 8; ++i) host[i] = m ? wB : wA;
         VecOp::eqV(fields[m], host);
      }
      system.w().setRGrid(fields, true);
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
   * Are all values of both c fields equal to a given value?
   */
   bool hasUniformC(System<1,T>& system, double value)
   {
      if (!system.c().hasData()) return false;
      Util::DArray<typename T::Real> host;
      host.allocate(8);
      for (int m = 0; m < 2; ++m) {
         VecOp::eqV(host, system.c().rgrid(m));
         for (int i = 0; i < 8; ++i) {
            if (!(std::abs(host[i] - value) < 1.0E-10)) return false;
         }
      }
      return true;
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
   * Range failure of the initial MDE solution (w_A = 1e306, ln Q
   * about -1e309), and of the MDE solution after the first AM step
   * (a mixing parameter lambda of 1e307 makes w_A about -2e306 in a
   * mixed ensemble, so that ln Q and phi overflow).
   */
   void checkIterate(std::string const & iterator)
   {
      // Initial state: failure, then repeated failure and recovery
      {
         System<1,T> system;
         read(system, param(iterator, false));
         setFields(system, 1.0E306, 0.0);
         TEST_ASSERT(iterate(system) == 1);
         TEST_ASSERT(!system.c().hasData());
         TEST_ASSERT(iterate(system) == 1);
         setFields(system, 0.0, 0.0);
         TEST_ASSERT(iterate(system) == 0);
         TEST_ASSERT(hasUniformC(system, 0.5));
      }

      // After the first AM step
      {
         System<1,T> system;
         read(system, param(iterator, true,
                            "lambda 1e307 useLambdaRamp 0"));
         setFields(system, 0.0, 0.0);
         TEST_ASSERT(iterate(system) == 1);
         TEST_ASSERT(!system.c().hasData());
         TEST_ASSERT(iterate(system) == 1);
         setFields(system, 0.0, 0.4);
         TEST_ASSERT(iterate(system) == 0);
      }
   }

   /*
   * ITERATE with a range failure gives status 1, and the commands
   * after it are not executed, as after a NaN error.
   */
   void checkCommand()
   {
      System<1,T> system;
      read(system, param("AmIteratorBasis", false));
      setFields(system, 1.0E306, 0.0);
      std::string name = "out/ScftRangeFailureCommand_w.bf";
      std::remove((this->filePrefix() + name).c_str());
      std::istringstream commands("ITERATE\nWRITE_W_BASIS " + name
                                  + "\nFINISH\n");
      int status = 0;
      try {
         status = system.readCommands(commands);
      } catch (Util::Exception& e) {
         status = -1;
      }
      TEST_ASSERT(status == 1);
      std::ifstream file((this->filePrefix() + name).c_str());
      TEST_ASSERT(!file.is_open());
   }

   /*
   * A SWEEP attempt that fails with a range failure backtracks. Here
   * every attempt with a nonzero initial error fails (lambda 1e307),
   * so that the sweep stops after the allowed number of step
   * reductions, with its own exception.
   */
   void checkSweep()
   {
      System<1,T> system;
      std::string sweep = "LinearSweep{ ns 1 "
                          "baseFileName out/ScftRangeFailureSweep_ "
                          "nParameter 1 parameters[ chi 0 1 +0.1 ] }";
      read(system, param("AmIteratorBasis", true,
                         "lambda 1e307 useLambdaRamp 0", sweep));
      setFields(system, 0.0, 0.4);
      std::string message;
      try {
         system.sweep();
      } catch (Util::Exception& e) {
         message = e.message();
      }
      TEST_ASSERT(message.find("Sweep decreased ds too many times")
                  != std::string::npos);
   }

};

}
}
#endif
