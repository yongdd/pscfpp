#ifndef R1_FIELD_SHIFT_TEST_H
#define R1_FIELD_SHIFT_TEST_H

#include "R1TestBase.h"

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/field/Domain.h>
#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <r1/solvers/Block.h>

#include <cmath>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

/*
* Field shift of the r1 polymer solver (Polymer::compute): each block is
* solved in w - s, s = min(w) - ln<exp(-N(w - min(w)))>/N. The shift sets
* the Crank-Nicolson contour error, and must keep the propagators in
* range for long chains in strongly nonuniform fields.
*/
class R1FieldShiftTest : public R1TestBase
{

public:

   /*
   * The documented shift of a field, for a polymer of length n.
   */
   double rule(System& sys, DArray<double> const & w, double n)
   {
      int nx = w.capacity();
      double wMin = w[0];
      for (int i = 1; i < nx; ++i) {
         if (w[i] < wMin) wMin = w[i];
      }
      DArray<double> e;
      e.allocate(nx);
      for (int i = 0; i < nx; ++i) e[i] = std::exp(-n*(w[i] - wMin));
      return wMin - std::log(sys.domain().spatialAverage(e))/n;
   }

   /*
   * Shipped spherical and cylindrical micelles (canonical, NR, ds 0.005):
   * the converged free energy lies within 1e-3 of its ds -> 0 limit, as
   * for the unshifted solver. The volume average shift gave -1.7e-3
   * (spherical) and -6.0e-3 (cylindrical) [v06].
   */
   void testRadialCanonical()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1FieldShiftRadial.log");

      char const * params[2] = {"in/R1/shift_sph_c1.prm",
                                "in/R1/shift_cyl.prm"};
      char const * fields[2] = {"READ_W in/R1/shift_sph.w\nFINISH\n",
                                "READ_W in/R1/shift_cyl.w\nFINISH\n"};
      // Richardson limits of ds 0.000625 and 0.0003125 (same nx)
      double limits[2] = {-0.9009112, -0.6923114};
      for (int k = 0; k < 2; ++k) {
         System sys;
         readSystem(sys, params[k]);
         runCommands(sys, fields[k]);
         TEST_ASSERT(sys.iterate() == 0);
         Log::file() << params[k] << ": f = " << sys.fHelmholtz()
                     << "  error = " << sys.fHelmholtz() - limits[k]
                     << std::endl;
         TEST_ASSERT(std::abs(sys.fHelmholtz() - limits[k]) < 1.0E-3);

         // Every block uses the documented shift
         Polymer const & polymer = sys.mixture().polymer(0);
         for (int j = 0; j < polymer.nBlock(); ++j) {
            Block const & block = polymer.block(j);
            double s = rule(sys, sys.wField(block.monomerId()),
                            polymer.length());
            TEST_ASSERT(std::abs(block.wShift() - s) < 1.0E-10);
         }
      }
   }

   /*
   * Shipped grand-canonical spherical micelle (examples/1d/spherical/1,
   * run.g1): AM converges again, and the concentration at the centre of
   * the shipped field is within 3e-3 of its ds -> 0 limit 1.001578
   * (1.2949 with the volume average shift) [v06].
   */
   void testRadialGrandCanonical()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1FieldShiftGrand.log");

      System sys;
      readSystem(sys, "in/R1/shift_sph_g1.prm");
      runCommands(sys, "READ_W in/R1/shift_sph.w\nFINISH\n");
      sys.compute();
      double c0 = sys.cField(0)[0] + sys.cField(1)[0];
      Log::file() << "sum c(r=0) = " << c0 << std::endl;
      TEST_ASSERT(std::abs(c0 - 1.001578) < 3.0E-3);

      TEST_ASSERT(sys.iterate() == 0);
      Log::file() << "f = " << sys.fHelmholtz() << std::endl;
      TEST_ASSERT(std::abs(sys.fHelmholtz() + 0.8956555) < 6.0E-3);
   }

   /*
   * Long chains in a nonuniform field (planar, nx 2, w = (0,2)): the
   * average shift overflowed (NaN). N = 1000 keeps the documented shift;
   * for N = 3000 it underflows (ln Q' about -837), and the shift is
   * raised so that Q' = 1 (before, the average shift again overflowed).
   * Reference: exact contour integrals of exponentials of the same 2x2
   * spatial operator.
   */
   void testLongChain()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1FieldShiftLong.log");

      char const * params[2] = {"in/R1/shift_long.prm",
                                "in/R1/shift_long3000.prm"};
      double lnQ[2] = {-279.6591572320251, -838.1407171197723};
      double c0[2] = {1.9484670702844968, 1.948611222128508};
      double c1[2] = {0.051532929715503044, 0.0513887778714918};
      double tol[2] = {5.0E-4, 1.5E-3};
      for (int k = 0; k < 2; ++k) {
         System sys;
         readSystem(sys, params[k]);
         sys.wField(0)[0] = 0.0;
         sys.wField(0)[1] = 2.0;
         sys.compute();
         Polymer const & polymer = sys.mixture().polymer(0);
         Log::file() << "lnQ = " << polymer.lnQ() << "  c = "
                     << sys.cField(0)[0] << "  " << sys.cField(0)[1]
                     << std::endl;
         TEST_ASSERT(std::isfinite(polymer.lnQ()));
         TEST_ASSERT(std::abs(polymer.lnQ() - lnQ[k]) < tol[k]);
         TEST_ASSERT(std::abs(sys.cField(0)[0] - c0[k]) < 5.0E-8);
         TEST_ASSERT(std::abs(sys.cField(0)[1] - c1[k]) < 5.0E-8);
         double s = rule(sys, sys.wField(0), polymer.length());
         if (k == 0) {
            TEST_ASSERT(std::abs(polymer.block(0).wShift() - s) < 1.0E-12);
         } else {
            // Q' is one up to the change of the Crank-Nicolson error
            TEST_ASSERT(polymer.block(0).wShift() > s);
            TEST_ASSERT(std::abs(std::log(polymer.qShifted())) < 1.0E-2);
         }
      }
   }

   /*
   * Two homopolymers in a step field 0|1500 (planar, nx 41, ds 0.01):
   * the average shift gave a negative Crank-Nicolson factor and mu
   * 674.8. Reference: Richardson limit ds -> 0 at the same nx.
   */
   void testStepField()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1FieldShiftStep.log");

      System sys;
      readSystem(sys, "in/R1/shift_step.prm");
      int nx = sys.domain().nx();
      for (int i = 0; i < nx; ++i) {
         sys.wField(0)[i] = (i <= nx/2) ? 0.0 : 1500.0;
         sys.wField(1)[i] = 1500.0 - sys.wField(0)[i];
      }
      sys.compute();
      double mu = sys.mixture().polymer(0).mu();
      Log::file() << "mu = " << mu << std::endl;
      TEST_ASSERT(std::abs(mu - 0.1991448112) < 5.0E-6);
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(std::isfinite(sys.cField(0)[i]));
         TEST_ASSERT(std::isfinite(sys.cField(1)[i]));
      }
   }

   /*
   * Diblock 1000 + 1000 (planar, nx 17, ds 0.1) with w_A = 0 and w_B
   * alternating 0|6: the documented shift underflows. Only the B shift
   * is raised (the A block is in a uniform field), so that Q' = 1 and
   * the A propagators stay bounded; the concentrations are finite and
   * conserve the amount of each block (NaN before).
   */
   void testDiblock()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1FieldShiftDiblock.log");

      System sys;
      readSystem(sys, "in/R1/shift_diblock.prm");
      int nx = sys.domain().nx();
      for (int i = 0; i < nx; ++i) {
         sys.wField(0)[i] = 0.0;
         sys.wField(1)[i] = (i % 2) ? 6.0 : 0.0;
      }
      sys.compute();
      Polymer const & polymer = sys.mixture().polymer(0);
      Log::file() << "lnQ = " << polymer.lnQ() << "  qShifted = "
                  << polymer.qShifted() << std::endl;
      TEST_ASSERT(std::isfinite(polymer.lnQ()));
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(std::isfinite(sys.cField(0)[i]));
         TEST_ASSERT(std::isfinite(sys.cField(1)[i]));
      }
      TEST_ASSERT(std::abs(sys.domain().spatialAverage(sys.cField(0)) - 0.5)
                  < 1.0E-10);
      TEST_ASSERT(std::abs(sys.domain().spatialAverage(sys.cField(1)) - 0.5)
                  < 1.0E-10);
      double n = polymer.length();
      TEST_ASSERT(std::abs(polymer.block(0).wShift()
                           - rule(sys, sys.wField(0), n)) < 1.0E-5);
      TEST_ASSERT(polymer.block(1).wShift() > rule(sys, sys.wField(1), n));
      TEST_ASSERT(std::abs(std::log(polymer.qShifted())) < 1.0E-2);
   }

   /*
   * Grand-canonical diblock with mu = 750 in a uniform field: phi =
   * exp(750) is not representable. Each iterator (AM, NR, BinaryRelax)
   * must report a failure to converge, not end the program with the
   * exception of the solver (SIGABRT in the first version of the
   * raised shifts; a NaN failure before).
   */
   void testIterateFailure()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1FieldShiftIterateFailure.log");

      char const * params[3] = {"in/R1/shift_gc_am.prm",
                                "in/R1/shift_gc_nr.prm",
                                "in/R1/shift_gc_br.prm"};
      for (int k = 0; k < 3; ++k) {
         System sys;
         readSystem(sys, params[k]);
         runCommands(sys, "READ_W in/R1/shift_gc.w\nFINISH\n");
         bool threw = false;
         int error = 0;
         try {
            error = sys.iterate();
         } catch (Exception& e) {
            threw = true;
         }
         TEST_ASSERT(!threw);
         TEST_ASSERT(error != 0);
      }
   }

   /*
   * Chain of total length 2400 in a field 8 with a single point of field
   * 0 (planar, nx 1025, ds 0.2), as one block or as 120 blocks of 20
   * (alternating monomer types in the same field). Calibrating each block
   * from a uniform head over-raises it by up to ln(1/f) per junction for
   * a chain confined to a small fraction f of the domain, which
   * overflowed for 120 blocks; the common raise then makes the chain the
   * same discrete operator as one block.
   */
   void testManyBlocks()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/R1FieldShiftManyBlocks.log");

      char const * params[2] = {"in/R1/shift_many1.prm",
                                "in/R1/shift_many120.prm"};
      double f[2], mu[2], c[2];
      for (int k = 0; k < 2; ++k) {
         System sys;
         readSystem(sys, params[k]);
         int nx = sys.domain().nx();
         for (int i = 0; i < nx; ++i) {
            double w = (i == nx/2) ? 0.0 : 8.0;
            sys.wField(0)[i] = w;
            sys.wField(1)[i] = w;
         }
         bool threw = false;
         try {
            sys.compute();
         } catch (Exception& e) {
            threw = true;
         }
         TEST_ASSERT(!threw);
         f[k] = sys.fHelmholtz();
         mu[k] = sys.mixture().polymer(0).mu();
         c[k] = sys.cField(0)[nx/2] + sys.cField(1)[nx/2];
         Log::file() << params[k] << ": f = " << f[k] << "  mu = "
                     << mu[k] << "  c = " << c[k] << std::endl;
         TEST_ASSERT(std::isfinite(f[k]) && std::isfinite(mu[k]));
      }
      TEST_ASSERT(std::abs(f[1] - f[0]) < 1.0E-8*std::abs(f[0]));
      TEST_ASSERT(std::abs(mu[1] - mu[0]) < 1.0E-8*std::abs(mu[0]));
      TEST_ASSERT(std::abs(c[1] - c[0]) < 1.0E-8*std::abs(c[0]));
   }

};

TEST_BEGIN(R1FieldShiftTest)
TEST_ADD(R1FieldShiftTest, testRadialCanonical)
TEST_ADD(R1FieldShiftTest, testRadialGrandCanonical)
TEST_ADD(R1FieldShiftTest, testLongChain)
TEST_ADD(R1FieldShiftTest, testStepField)
TEST_ADD(R1FieldShiftTest, testDiblock)
TEST_ADD(R1FieldShiftTest, testIterateFailure)
TEST_ADD(R1FieldShiftTest, testManyBlocks)
TEST_END(R1FieldShiftTest)

#endif
