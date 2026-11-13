#ifndef RPC_ENVIRONMENT_STRESS_INVALIDATION_TEST_H
#define RPC_ENVIRONMENT_STRESS_INVALIDATION_TEST_H

#include "WallFieldScaleTest.h"

class EnvironmentStressInvalidationTest : public WallFieldScaleTestBase
{
public:
   void testFieldReplacementStress()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/EnvironmentStressInvalidationEnvironmentFieldReplacement.log");
      Rp::System<1,CPT> reused;
      initialize(reused, 1.0, 1.0);
      setFields(reused, 0.0);
      reused.compute(true);
      double firstStress = reused.environment().stress(0);

      int count = 0;
      for (double amplitude : {3.0, -2.0, 0.0}) {
         setFields(reused, amplitude);
         TEST_ASSERT(!reused.environment().hasStress());
         TEST_ASSERT(!reused.environment().needsUpdate());
         TEST_ASSERT(!reused.c().hasData());
         if (count % 2 == 0) {
            reused.compute(true);
         } else {
            reused.compute();
            TEST_ASSERT(!reused.environment().hasStress());
            reused.computeStress();
         }
         Rp::System<1,CPT> fresh;
         initialize(fresh, 1.0, 1.0);
         setFields(fresh, amplitude);
         fresh.compute(true);
         double expected = fresh.environment().stress(0);
         TEST_ASSERT(std::abs(reused.environment().stress(0)-expected) < 1e-11);
         TEST_ASSERT(std::abs(reused.mixture().stress(0)
                            -fresh.mixture().stress(0)) < 1e-11);
         for (int m = 0; m < 2; ++m) {
            std::vector<double> a = values(reused.c().rgrid(m));
            std::vector<double> b = values(fresh.c().rgrid(m));
            for (int i = 0; i < 64; ++i) {
               TEST_ASSERT(std::abs(a[i]-b[i]) < 1e-12);
            }
         }
         if (count == 0) TEST_ASSERT(std::abs(expected-firstStress) > 0.1);

         // Stress-only clearing is repeatable and leaves generated fields.
         reused.environment().clearStress();
         reused.environment().clearStress();
         TEST_ASSERT(!reused.environment().hasStress());
         TEST_ASSERT(!reused.environment().needsUpdate());
         reused.computeStress();
         TEST_ASSERT(std::abs(reused.environment().stress(0)-expected) < 1e-11);

         // A full reset retains its separate field-regeneration contract.
         reused.environment().reset();
         TEST_ASSERT(reused.environment().needsUpdate());
         reused.environment().generate();
         reused.computeStress();
         TEST_ASSERT(std::abs(reused.environment().stress(0)-expected) < 1e-11);
         ++count;
      }
   }
};

TEST_BEGIN(EnvironmentStressInvalidationTest)
TEST_ADD(EnvironmentStressInvalidationTest, testFieldReplacementStress)
TEST_END(EnvironmentStressInvalidationTest)

#endif
