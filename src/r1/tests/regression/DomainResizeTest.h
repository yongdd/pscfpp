#ifndef R1_DOMAIN_RESIZE_TEST_H
#define R1_DOMAIN_RESIZE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <r1/field/Domain.h>

using namespace Util;
using namespace Pscf::R1;

class DomainResizeTest : public UnitTest
{
public:
   void testReuseAfterResize()
   {
      printMethod(TEST_FUNC);
      Domain reused;
      int sizes[] = {9, 5, 17, 17, 3};
      for (int mode = 0; mode < 3; ++mode) {
         for (int n : sizes) {
            Domain fresh;
            if (mode == 0) {
               reused.setPlanarParameters(0.2, 2.0, n);
               fresh.setPlanarParameters(0.2, 2.0, n);
            } else if (mode == 1) {
               reused.setCylinderParameters(2.0, n);
               fresh.setCylinderParameters(2.0, n);
            } else {
               reused.setSphereParameters(2.0, n);
               fresh.setSphereParameters(2.0, n);
            }
            Domain::Field f, g;
            f.allocate(n);
            g.allocate(n);
            for (int i = 0; i < n; ++i) {
               f[i] = 1.5;
               g[i] = 2.0;
            }
            TEST_ASSERT(eq(reused.innerProduct(f, g), 3.0));
            for (int i = 0; i < n; ++i) {
               f[i] = 0.2 + double(i)/double(n - 1);
               g[i] = 1.0 - 0.3*f[i];
            }
            TEST_ASSERT(eq(reused.innerProduct(f, g),
                           fresh.innerProduct(f, g)));
         }
      }
   }
};

TEST_BEGIN(DomainResizeTest)
TEST_ADD(DomainResizeTest, testReuseAfterResize)
TEST_END(DomainResizeTest)
#endif
