#ifndef RPC_SCFT_SCALE_STRESS_TEST_H
#define RPC_SCFT_SCALE_STRESS_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/scft/iterator/Iterator.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>

#include <util/misc/Exception.h>
#include <util/param/Label.h>
#include <util/tests/LogFileUnitTest.h>

#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests of the scaleStress parameter of the AM iterators.
*/
class ScftScaleStressTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * A flexible AmIteratorBasis or AmIteratorGrid rejects a scaleStress
   * that is zero or negative. With scaleStress 0, the cell did not
   * enter the residual and was restored as 0/0: ITERATE reported
   * convergence with the cell unchanged, or gave NaN cell parameters.
   * Positive values, and any value for a rigid unit cell (where
   * scaleStress is not used), are accepted.
   */
   void testScaleStress()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testScaleStress.log");

      // Rejected
      TEST_ASSERT(readThrows("AmIteratorBasis", 1, "0.0"));
      TEST_ASSERT(readThrows("AmIteratorBasis", 1, "-10.0"));
      TEST_ASSERT(readThrows("AmIteratorGrid", 1, "0.0"));
      TEST_ASSERT(readThrows("AmIteratorGrid", 1, "-10.0"));

      // Accepted
      TEST_ASSERT(!readThrows("AmIteratorBasis", 1, "1.0E-6"));
      TEST_ASSERT(!readThrows("AmIteratorGrid", 1, "1.0E+6"));
      TEST_ASSERT(!readThrows("AmIteratorBasis", 0, "0.0"));
      TEST_ASSERT(!readThrows("AmIteratorGrid", 0, "-1.0"));

      // A small positive scaleStress still optimizes the cell
      Rp::System<1,CPT> system;
      TEST_ASSERT(!readParamThrows(system,
                                   param("AmIteratorBasis", 1, "0.1")));
      system.fileMaster().setInputPrefix(filePrefix());
      system.w().readBasis("in/Scft/lam_w.bf");
      int error = system.iterate();
      TEST_ASSERT(!error);
      TEST_ASSERT(system.iterator().isFlexible());
      TEST_ASSERT(system.mixture().hasStress());
      TEST_ASSERT(std::abs(system.mixture().stress(0)) < 1.0E-6);
   }

   // Utility functions

   /*
   * Parameter file of a lamellar diblock melt with an AM iterator.
   */
   std::string param(std::string const & iteratorName, int isFlexible,
                     std::string const & scaleStress)
   {
      std::ostringstream out;
      out << "System{\n"
          << "  Mixture{\n"
          << "    nMonomer  2\n"
          << "    monomers[\n"
          << "              1.0\n"
          << "              1.0\n"
          << "    ]\n"
          << "    nPolymer  1\n"
          << "    Polymer{\n"
          << "      type    linear\n"
          << "      nBlock  2\n"
          << "      blocks[\n"
          << "              0     0.56\n"
          << "              1     0.44\n"
          << "      ]\n"
          << "      phi     1.0\n"
          << "    }\n"
          << "    ds   0.01\n"
          << "  }\n"
          << "  Interaction{\n"
          << "    chi(\n"
          << "         1   0   12.0\n"
          << "    )\n"
          << "  }\n"
          << "  Domain{\n"
          << "    mesh            40\n"
          << "    lattice         lamellar\n"
          << "    groupName       P_-1\n"
          << "  }\n"
          << "  " << iteratorName << "{\n"
          << "    epsilon       1.0e-9\n"
          << "    maxItr        300\n"
          << "    maxHist       10\n"
          << "    isFlexible    " << isFlexible << "\n"
          << "    scaleStress   " << scaleStress << "\n"
          << "  }\n"
          << "}\n";
      return out.str();
   }

   /*
   * Read a parameter file given as a string. Returns true if a
   * Util::Exception was thrown.
   */
   bool readParamThrows(Rp::System<1,CPT>& system,
                        std::string const & text)
   {
      std::istringstream in(text);
      try {
         system.readParam(in);
      } catch (Exception& e) {
         Label::clear(); // Discard any label left by the failed read
         return true;
      }
      return false;
   }

   /*
   * Does reading a parameter file with this iterator block throw?
   */
   bool readThrows(std::string const & iteratorName, int isFlexible,
                   std::string const & scaleStress)
   {
      Rp::System<1,CPT> system;
      return readParamThrows(system,
                             param(iteratorName, isFlexible, scaleStress));
   }

};

TEST_BEGIN(ScftScaleStressTest)
TEST_ADD(ScftScaleStressTest, testScaleStress)
TEST_END(ScftScaleStressTest)

#endif
