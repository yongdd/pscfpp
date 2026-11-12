#ifndef R1_CORE_MIXTURE_TEST_H
#define R1_CORE_MIXTURE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <r1/solvers/Solvent.h>
#include <util/param/ParamComponent.h>
#include <util/tests/LogFileUnitTest.h>

#include <cmath>
#include <fstream>
#include <string>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

/*
* Input checks of a mixture: in/Core/mix_valid is accepted, and each of
* the other in/Core/mix_* files differs from it by one invalid value,
* which must be rejected with an Exception when the file is read.
*/
class CoreMixtureTest : public LogFileUnitTest 
{

public:

   void setUp()
   {  ParamComponent::setEcho(true); }

   void tearDown()
   {
      closeLogFile();
      ParamComponent::setEcho(false);
   }

   void testReadParametersValid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/Core/MixtureReadParametersValid.log");

      std::ifstream in;
      openInputFile("in/Core/mix_valid", in);
      Mixture mix;
      mix.readParam(in);
      TEST_ASSERT(mix.polymer(0).phi() == 0.6);
      TEST_ASSERT(mix.solvent(0).phi() == 0.4);
      TEST_ASSERT(eq(mix.vMonomer(), 0.05));
      TEST_ASSERT(eq(mix.solvent(0).size(), 0.05));
   }

   // Solvent monomerId 2 with nMonomer 2 (it is used as an index)
   void testReadParametersSolventMonomerId()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/Core/MixtureReadParametersSolventMonomerId.log");
      TEST_ASSERT(readParamThrows("in/Core/mix_solventId"));
   }

   // Block monomerId -1 (it is used as an index)
   void testReadParametersBlockMonomerId()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/Core/MixtureReadParametersBlockMonomerId.log");
      TEST_ASSERT(readParamThrows("in/Core/mix_blockId"));
   }

   // Closed ensemble, phi = 0.75, 0.5 and -0.25 (sum 1)
   void testReadParametersPhiNegative()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/Core/MixtureReadParametersPhiNegative.log");
      TEST_ASSERT(readParamThrows("in/Core/mix_phiNegative"));
   }

   // Closed ensemble, phi = 0.6 and 0.5 (sum 1.1)
   void testReadParametersPhiSum()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/Core/MixtureReadParametersPhiSum.log");
      TEST_ASSERT(readParamThrows("in/Core/mix_phiSum"));
   }

   // Solvent size 0
   void testReadParametersSolventSize()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/Core/MixtureReadParametersSolventSize.log");
      TEST_ASSERT(readParamThrows("in/Core/mix_size"));
   }

   // vMonomer 0
   void testReadParametersVMonomer()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/Core/MixtureReadParametersVMonomer.log");
      TEST_ASSERT(readParamThrows("in/Core/mix_vMonomer"));
   }

   /*
   * Closed ensemble, phi = 0.6 and 0.400001: a sum within 1e-5 of 1
   * (e.g., fractions written with few digits) is accepted and rescaled
   * to a sum of 1, keeping the ratios of the volume fractions.
   */
   void testReadParametersPhiRescale()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/Core/MixtureReadParametersPhiRescale.log");

      std::ifstream in;
      openInputFile("in/Core/mix_phiRescale", in);
      Mixture mix;
      mix.readParam(in);

      double phiP = mix.polymer(0).phi();
      double phiS = mix.solvent(0).phi();
      TEST_ASSERT(std::abs(phiP + phiS - 1.0) < 1.0E-14);
      TEST_ASSERT(std::abs(phiP/phiS - 0.6/0.400001) < 1.0E-12);
      TEST_ASSERT(std::abs(phiP - 0.6/1.000001) < 1.0E-14);
   }

private:

   /*
   * Read a Mixture from a file. Return true if this throws an Exception.
   */
   bool readParamThrows(std::string fileName)
   {
      std::ifstream in;
      openInputFile(fileName, in);
      Mixture mix;
      bool thrown = false;
      try {
         mix.readParam(in);
      } catch (Util::Exception&) {
         thrown = true;
      }
      in.close();
      return thrown;
   }

};

TEST_BEGIN(CoreMixtureTest)
TEST_ADD(CoreMixtureTest, testReadParametersValid)
TEST_ADD(CoreMixtureTest, testReadParametersSolventMonomerId)
TEST_ADD(CoreMixtureTest, testReadParametersBlockMonomerId)
TEST_ADD(CoreMixtureTest, testReadParametersPhiNegative)
TEST_ADD(CoreMixtureTest, testReadParametersPhiSum)
TEST_ADD(CoreMixtureTest, testReadParametersSolventSize)
TEST_ADD(CoreMixtureTest, testReadParametersVMonomer)
TEST_ADD(CoreMixtureTest, testReadParametersPhiRescale)
TEST_END(CoreMixtureTest)

#endif
