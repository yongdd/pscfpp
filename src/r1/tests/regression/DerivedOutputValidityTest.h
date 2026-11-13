#ifndef R1_DERIVED_OUTPUT_VALIDITY_TEST_H
#define R1_DERIVED_OUTPUT_VALIDITY_TEST_H

#include "R1TestBase.h"
#include <test/UnitTestRunner.h>
#include <r1/field/FieldIo.h>
#include <r1/field/Domain.h>
#include <limits>

class DerivedOutputValidityTest : public R1TestBase
{
public:
   void testReplacementGuards()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/DerivedOutputValidityReplacement.log");
      System sys;
      readSystem(sys, "in/R1/planar_nr2.prm");
      const std::string outputs[] = {
         "WRITE_C out/derivedoutputvalidity.c\n", "WRITE_BLOCK_C out/derivedoutputvalidity.block\n",
         "WRITE_THERMO out/derivedoutputvalidity.thermo\n",
         "COMPARE_HOMOGENEOUS 0 out/derivedoutputvalidity.homo\n",
         "WRITE_Q_SLICE out/derivedoutputvalidity.slice 0 0 0 0\n",
         "WRITE_Q_TAIL out/derivedoutputvalidity.tail 0 0 0\n",
         "WRITE_Q_VERTEX 0 0 out/derivedoutputvalidity.vertex\n",
         "WRITE_Q out/derivedoutputvalidity.q 0 0 0\n", "WRITE_Q_ALL out/derivedoutputvalidity.all\n"
      };
      runCommands(sys, "READ_W in/R1/planar.w\nCOMPUTE\nFINISH\n");
      // Successful derived output before replacement is a control for every
      // command; none of these calls should require SCFT convergence.
      for (int i = 0; i < 9; ++i) runCommands(sys, outputs[i]);
      const std::string oldC = fileContents("out/derivedoutputvalidity.c");
      runCommands(sys, "READ_W in/R1/planar.w\nFINISH\n");
      for (int i = 0; i < 9; ++i) {
         bool thrown = false;
         try { runCommands(sys, outputs[i]); }
         catch (Exception&) { thrown = true; }
         TEST_ASSERT(thrown);
      }
      TEST_ASSERT(fileContents("out/derivedoutputvalidity.c") == oldC);
      // Recompute in a later command call: validity persists across calls.
      runCommands(sys, "COMPUTE\nFINISH\n");
      for (int i = 0; i < 9; ++i) runCommands(sys, outputs[i]);
      TEST_ASSERT(fileContents("out/derivedoutputvalidity.c") == oldC);
   }

private:

};

TEST_BEGIN(DerivedOutputValidityTest)
TEST_ADD(DerivedOutputValidityTest, testReplacementGuards)
TEST_END(DerivedOutputValidityTest)

#endif
