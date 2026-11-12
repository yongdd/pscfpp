#ifndef R1_SYSTEM_TEST_H
#define R1_SYSTEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/field/Domain.h>
#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <r1/iterator/Iterator.h>
#include <r1/field/FieldIo.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/chem/Vertex.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

class SystemTest : public UnitTest 
{

private:

   std::ofstream logFile_;

public:

   void setUp()
   {}

   void tearDown()
   {
      if (logFile_.is_open()) {
         logFile_.close();
      }
      ParamComponent::setEcho(false);
      setVerbose(0);
   }

   void openLogFile(char const * filename)
   {
      openOutputFile(filename, logFile_);
      Log::setFile(logFile_);
      ParamComponent::setEcho(true);
   }

   /*
   * Solve the SCFT equations with the iterator of the system. Return
   * true if the iterator converged, false if it failed or threw.
   */
   bool solve(System& sys)
   {
      int error = 1;
      try {
         error = sys.iterator().solve();
      } catch (Exception& e) {
         Log::file() << "Exception: " << e.message() << std::endl;
         return false;
      }
      return (error == 0);
   }

   /*
   * Execute commands. Check that every command was executed:
   * readCommands stops early, without an error, after an ITERATE that
   * failed or at an unknown command.
   */
   void runCommands(System& sys, std::istream& in)
   {
      bool hasThrown = false;
      try {
         sys.readCommands(in);
      } catch (Exception& e) {
         Log::file() << "Exception: " << e.message() << std::endl;
         hasThrown = true;
      }
      TEST_ASSERT(!hasThrown);

      // Nothing may be left to read after FINISH, the last command
      std::string rest;
      in >> rest;
      TEST_ASSERT(in.fail());
   }

   /*
   * Execute a command file, as runCommands.
   */
   void runCommandFile(System& sys, char const * filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      runCommands(sys, in);
      in.close();
   }

   /*
   * Check that the partition function Q of each polymer is the same
   * from every propagator (product of its head and the tail of its
   * partner), from every vertex (product of the incoming tails), and
   * as stored by the polymer, which normalizes the concentrations.
   */
   void checkQ(System& sys)
   {
      Mixture& mix = sys.mixture();
      Domain& domain = sys.domain();
      int nx = domain.nx();
      DArray<double> product;
      product.allocate(nx);
      double q0, q;
      int i, j;
      for (int ip = 0; ip < mix.nPolymer(); ++ip) {
         Polymer& polymer = mix.polymer(ip);
         polymer.propagator(0, 0).computeQ(q0);
         Log::file() << "Polymer " << ip << ", Q = " << q0 << "\n";
         TEST_ASSERT(q0 > 0.0);
         for (int ib = 0; ib < polymer.nBlock(); ++ib) {
            for (int id = 0; id < 2; ++id) {
               polymer.propagator(ib, id).computeQ(q);
               TEST_ASSERT(std::abs(q/q0 - 1.0) < 1.0E-6);
            }
         }
         for (int iv = 0; iv < polymer.nVertex(); ++iv) {
            Vertex const & vertex = polymer.vertex(iv);
            for (i = 0; i < nx; ++i) {
               product[i] = 1.0;
            }
            for (j = 0; j < vertex.size(); ++j) {
               Pair<int> const & pId = vertex.inPropagatorId(j);
               DArray<double> const & tail
                                 = polymer.propagator(pId[0], pId[1]).tail();
               for (i = 0; i < nx; ++i) {
                  product[i] *= tail[i];
               }
            }
            q = domain.spatialAverage(product);
            TEST_ASSERT(std::abs(q/q0 - 1.0) < 1.0E-6);
         }
         TEST_ASSERT(std::abs(polymer.q()/q0 - 1.0) < 1.0E-6);
      }
   }

   /*
   * Check the spatial average of each monomer concentration: the sum
   * over polymers of phi times the fraction of the chain length in
   * blocks of that monomer type (for any w fields, in a mixture of
   * closed polymer species without solvent).
   */
   void checkSumRules(System& sys)
   {
      Mixture& mix = sys.mixture();
      UTIL_CHECK(mix.nSolvent() == 0);
      int nm = mix.nMonomer();
      DArray<double> expected;
      expected.allocate(nm);
      int im;
      for (im = 0; im < nm; ++im) {
         expected[im] = 0.0;
      }
      for (int ip = 0; ip < mix.nPolymer(); ++ip) {
         Polymer& polymer = mix.polymer(ip);
         for (int ib = 0; ib < polymer.nBlock(); ++ib) {
            Block& block = polymer.block(ib);
            expected[block.monomerId()]
                       += polymer.phi()*block.length()/polymer.length();
         }
      }
      double average;
      for (im = 0; im < nm; ++im) {
         average = sys.domain().spatialAverage(sys.cField(im));
         Log::file() << "Average c[" << im << "] = " << average
                     << ", expected " << expected[im] << "\n";
         TEST_ASSERT(std::abs(average - expected[im]) < 1.0E-6);
      }
   }

   /*
   * Check that the current w fields solve the SCFT equations: compute
   * the concentrations again, then check incompressibility (the sum of
   * the c fields is 1) and the exchange equations (w_i - sum_j chi_ij
   * c_j is the same for every monomer type i) at every grid point,
   * against an absolute tolerance, and the checks of checkQ and
   * checkSumRules. NR stops when the largest residual is below its
   * epsilon; the tolerances used below for NR are 10 or 100 times
   * epsilon. AM stops when the norm of its residual r is below epsilon
   * (1e-7) times the norm of all w fields (about 300 in the planar
   * tests), and for two monomer types the errors checked here are at
   * most 2 max|r_i|, so at most about 6e-5; the AM tolerance is 1e-4.
   * (A tolerance proportional to the norm of w would also accept a
   * diverged state with huge w, which AM can report as converged.)
   */
   void checkScft(System& sys, double tolerance)
   {
      sys.mixture().compute(sys.wFields(), sys.cFields());
      DMatrix<double> const & chi = sys.interaction().chi();
      int nm = sys.mixture().nMonomer();
      int nx = sys.domain().nx();
      double incompressMax = 0.0;
      double exchangeMax = 0.0;
      double sum, r0, r, value;
      int i, j, k;
      for (k = 0; k < nx; ++k) {
         sum = 0.0;
         for (i = 0; i < nm; ++i) {
            sum += sys.cField(i)[k];
         }
         value = std::abs(sum - 1.0);
         if (!(value <= incompressMax)) incompressMax = value;
         r0 = 0.0;
         for (i = 0; i < nm; ++i) {
            r = sys.wField(i)[k];
            for (j = 0; j < nm; ++j) {
               r -= chi(i, j)*sys.cField(j)[k];
            }
            if (i == 0) {
               r0 = r;
            } else {
               value = std::abs(r - r0);
               if (!(value <= exchangeMax)) exchangeMax = value;
            }
         }
      }
      Log::file() << "Tolerance                   = " << tolerance
                  << "\n";
      Log::file() << "Max incompressibility error = " << incompressMax
                  << "\n";
      Log::file() << "Max exchange error          = " << exchangeMax
                  << "\n";
      TEST_ASSERT(incompressMax < tolerance);
      TEST_ASSERT(exchangeMax < tolerance);
      checkQ(sys);
      checkSumRules(sys);
   }

   /*
   * Return the maximum difference of the c fields of two systems.
   */
   double maxDifferenceC(System& a, System& b)
   {
      int nm = a.mixture().nMonomer();
      int nx = a.domain().nx();
      UTIL_CHECK(b.mixture().nMonomer() == nm);
      UTIL_CHECK(b.domain().nx() == nx);
      double value, maxValue = 0.0;
      for (int im = 0; im < nm; ++im) {
         for (int i = 0; i < nx; ++i) {
            value = std::abs(a.cField(im)[i] - b.cField(im)[i]);
            if (!(value <= maxValue)) maxValue = value;
         }
      }
      Log::file() << "Max difference of c fields = " << maxValue << "\n";
      return maxValue;
   }

   /*
   * Set the initial w fields of the planar tests for chi = 20.
   */
   void setPlanarGuess(System& sys)
   {
      int nx = sys.domain().nx();
      double cs;
      double chi = 20.0;
      for (int i = 0; i < nx; ++i) {
         cs = cos(Constants::Pi*double(i)/double(nx-1));
         sys.wField(0)[i] = chi*(-0.5*cs + 0.25*cs*cs);
         sys.wField(1)[i] = chi*(+0.5*cs + 0.25*cs*cs);
      }
      double shift = sys.wField(1)[nx-1];
      for (int i = 0; i < nx; ++i) {
         sys.wField(0)[i] -= shift;
         sys.wField(1)[i] -= shift;
      }
   }
  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      System sys;
   }

   void testReadParameters()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestReadParameters.log");

      std::ifstream in;
      openInputFile("in/planar_nr1.prm", in);

      System sys;
      sys.readParam(in);

      Log::file() << "\n";
      sys.writeParam(Log::file());

      TEST_ASSERT(sys.domain().mode() == Planar);
   }

   void testSolveMdePlanar()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestSolveMdePlanar.log");

      std::ifstream in;
      openInputFile("in/planar_nr1.prm", in);

      System sys;
      sys.readParam(in);

      Log::file() << "\n";
      sys.writeParam(Log::file());

      Mixture& mix = sys.mixture();
      Domain& domain = sys.domain();

      double nx = (double)domain.nx();
      double cs;
      for (int i = 0; i < nx; ++i) {
         cs = cos(2.0*Constants::Pi*(double(i)+0.5)/double(nx-1));
         sys.wField(0)[i] = 0.5 + cs;
         sys.wField(1)[i] = 0.7 - cs;
      }
      mix.compute(sys.wFields(), sys.cFields());

      // Test if same Q is obtained from different methods
      double q00, q01, q10, q11;
      mix.polymer(0).propagator(0, 0).computeQ(q00);
      mix.polymer(0).propagator(0, 1).computeQ(q01);
      mix.polymer(0).propagator(1, 0).computeQ(q10);
      mix.polymer(0).propagator(1, 1).computeQ(q11);

      setVerbose(1);
      if (verbose() > 0) {
         Log::file() << q00 << "\n";
         Log::file() << q01 << "\n";
         Log::file() << q10 << "\n";
         Log::file() << q11 << "\n";
      }
      TEST_ASSERT(std::abs((q01 - q00)/q00) < 1.0E-6);
      TEST_ASSERT(std::abs((q10 - q00)/q00) < 1.0E-6);
      TEST_ASSERT(std::abs((q11 - q00)/q00) < 1.0E-6);
      checkQ(sys);

      // Test spatial integral of block concentration
      double sum0 = domain.spatialAverage(sys.cField(0));
      double sum1 = domain.spatialAverage(sys.cField(1));
      Log::file() << "Volume fraction of block 0 = " << sum0 << "\n";
      Log::file() << "Volume fraction of block 1 = " << sum1 << "\n";
      checkSumRules(sys);

      TEST_ASSERT(eq(mix.polymer(0).length(), 5.0));
   }


   void testSolveMdeSpherical()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestSolveMdeSpherical.log");

      std::ifstream in;
      openInputFile("in/spherical_nr1.prm", in);

      System sys;
      sys.readParam(in);

      TEST_ASSERT( !sys.domain().isShell() );

      Log::file() << "\n";
      sys.writeParam(Log::file());

      Mixture& mix = sys.mixture();
      Domain& domain = sys.domain();

      double nx = (double)domain.nx();
      double cs;
      for (int i = 0; i < nx; ++i) {
         cs = cos(2.0*Constants::Pi*double(i)/double(nx-1));
         sys.wField(0)[i] = -cs;
         sys.wField(1)[i] = +cs;
      }
      mix.compute(sys.wFields(), sys.cFields());

      // Test if same Q is obtained from different methods
      double q00, q01, q10, q11;
      mix.polymer(0).propagator(0, 0).computeQ(q00);
      mix.polymer(0).propagator(0, 1).computeQ(q01);
      mix.polymer(0).propagator(1, 0).computeQ(q10);
      mix.polymer(0).propagator(1, 1).computeQ(q11);
      if (verbose() > 1) {
         Log::file() << q00 << "\n";
         Log::file() << q01 << "\n";
         Log::file() << q10 << "\n";
         Log::file() << q11 << "\n";
      }
      TEST_ASSERT(std::abs((q01 - q00)/q00) < 1.0E-6);
      TEST_ASSERT(std::abs((q10 - q00)/q00) < 1.0E-6);
      TEST_ASSERT(std::abs((q11 - q00)/q00) < 1.0E-6);
      checkQ(sys);

      // Test spatial integral of block concentration
      double sum0 = domain.spatialAverage(sys.cField(0));
      double sum1 = domain.spatialAverage(sys.cField(1));
      Log::file() << "Volume fraction of block 0 = " << sum0 << "\n";
      Log::file() << "Volume fraction of block 1 = " << sum1 << "\n";
      checkSumRules(sys);
   }

   /*
   * Test NR iterator on input field generated within function.
   */
   void testIteratorPlanarNr1()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorPlanarNr1.log");

      std::ifstream in;
      openInputFile("in/planar_nr2.prm", in);

      System sys;
      sys.readParam(in);
      FieldIo fieldIo;
      fieldIo.associate(sys.domain(), sys.fileMaster());

      Log::file() << "\n";
      sys.writeParam(Log::file());

      Mixture& mix = sys.mixture();
      setPlanarGuess(sys);

      // Compute initial state
      mix.compute(sys.wFields(), sys.cFields());

      std::ofstream out;
      openOutputFile("out/initialPlanarNr1.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/initialPlanarNr1.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

      TEST_ASSERT(solve(sys));

      openOutputFile("out/finalPlanarNr1.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/finalPlanarNr1.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

      checkScft(sys, 1.0E-5);

   }

   /*
   * Test NR iterator on input w field read from file.
   */
   void testIteratorPlanarNr2()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorPlanarNr2.log");

      std::ifstream in;
      openInputFile("in/planar_nr2.prm", in);

      System sys;
      sys.readParam(in);
      in.close();

      FieldIo fieldIo;
      fieldIo.associate(sys.domain(), sys.fileMaster());

      Log::file() << "\n";
      // sys.writeParam(Log::file());

      openInputFile("in/planar.w", in);
      fieldIo.readFields(sys.wFields(), in);
      in.close();

      TEST_ASSERT(solve(sys));

      std::ofstream out;
      openOutputFile("out/planarNr2.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

      checkScft(sys, 1.0E-5);

   }

   /*
   * Test NR iterator controlled by a command file.
   */
   void testIteratorPlanarNr3()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorPlanarNr3.log");

      System sys;
      std::ifstream in;
      Log::file() << "\n";

      openInputFile("in/planar_nr2.prm", in);
      sys.readParam(in);
      in.close();

      // Set System filemaster prefixes to unit test file prefix
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      runCommandFile(sys, "in/planar_nr.cmd");
      checkScft(sys, 1.0E-5);
   }

   /*
   * Test NR iterator on input field generated within function.
   */
   void testIteratorSphericalNr1()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorSphericalNr1.log");

      std::ifstream in;
      openInputFile("in/spherical_nr1.prm", in);

      System sys;
      sys.readParam(in);

      Log::file() << "\n";
      sys.writeParam(Log::file());

      Mixture& mix = sys.mixture();
      Domain& domain = sys.domain();
      FieldIo fieldIo;
      fieldIo.associate(sys.domain(), sys.fileMaster());

      // Create initial chemical potential fields
      double nx = (double)domain.nx();
      double cs;
      double chi = 80.0;
      for (int i = 0; i < nx; ++i) {
         cs = cos(Constants::Pi*double(i)/double(nx-1));
         sys.wField(0)[i] = -chi*cs/2.0;
         sys.wField(1)[i] = +chi*cs/2.0;
      }
      double shift = sys.wField(1)[nx-1];
      for (int i = 0; i < nx; ++i) {
         sys.wField(0)[i] -= shift;
         sys.wField(1)[i] -= shift;
      }

      // Solve MDE for initial fields
      mix.compute(sys.wFields(), sys.cFields());

      Log::file() << "Average fraction 0 = " 
                << domain.spatialAverage(sys.cField(0)) << "\n";
      Log::file() << "Average fraction 1 = " 
                << domain.spatialAverage(sys.cField(1)) << "\n";

      std::ofstream out;
      openOutputFile("out/initialSphericalNr1.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/initialSphericalNr1.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

      TEST_ASSERT(solve(sys));

      openOutputFile("out/finalSphericalNr1.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/finalSphericalNr1.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

      checkScft(sys, 1.0E-4);

   }

   void testIteratorSphericalNr3()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorSphericalNr3.log");

      System sys;
      std::ifstream in;
      Log::file() << "\n";

      openInputFile("in/spherical_nr2.prm", in);
      sys.readParam(in);
      in.close();

      TEST_ASSERT(!sys.domain().isShell());

      // Set System filemaster prefixes to unit test file prefix
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      runCommandFile(sys, "in/spherical_nr.cmd");
      checkScft(sys, 1.0E-6);
   }

   void testIteratorPlanarAm1()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorPlanarAm.log");

      std::ifstream in;
      openInputFile("in/planar_am2.prm", in);

      System sys;
      sys.readParam(in);
      FieldIo fieldIo;
      fieldIo.associate(sys.domain(), sys.fileMaster());

      Log::file() << "\n";
      sys.writeParam(Log::file());

      Mixture& mix = sys.mixture();
      setPlanarGuess(sys);

      // Compute initial state
      mix.compute(sys.wFields(), sys.cFields());

      std::ofstream out;
      openOutputFile("out/initialPlanarAm.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/initialPlanarAm.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

      TEST_ASSERT(solve(sys));

      openOutputFile("out/finalPlanarAm.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/finalPlanarAm.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

      checkScft(sys, 1.0E-4);

      // The NR iterator must find the same solution of the same problem
      // from the same initial fields (free energy and concentrations,
      // which do not depend on the gauge of w)
      System sysNr;
      in.close();
      openInputFile("in/planar_nr2.prm", in);
      sysNr.readParam(in);
      in.close();
      setPlanarGuess(sysNr);
      TEST_ASSERT(solve(sysNr));
      sysNr.computeFreeEnergy();
      sys.computeFreeEnergy();
      Log::file() << "fHelmholtz (AM, NR) = " << sys.fHelmholtz() << "  "
                  << sysNr.fHelmholtz() << "\n";
      TEST_ASSERT(std::abs(sys.fHelmholtz() - sysNr.fHelmholtz()) < 1.0E-6);
      TEST_ASSERT(maxDifferenceC(sys, sysNr) < 1.0E-5);

   }

   /*
   * Test AM iterator controlled by a command file.
   */
   void testIteratorPlanarAm3()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorPlanarAm3.log");

      System sys;
      std::ifstream in;
      Log::file() << "\n";

      openInputFile("in/planar_am2.prm", in);
      sys.readParam(in);
      in.close();

      // Set System filemaster prefixes to unit test file prefix
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      runCommandFile(sys, "in/planar_am.cmd");
      checkScft(sys, 1.0E-4);
   }

   void testSweepSpherical()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestSweepSpherical.log");

      System sys;
      std::ifstream in;
      Log::file() << "\n";

      openInputFile("in/spherical3_nr.prm", in);
      sys.readParam(in);
      in.close();
      sys.writeParam(Log::file());

      TEST_ASSERT( !sys.domain().isShell() );

      // Set System filemaster prefixes to unit test file prefix
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      runCommandFile(sys, "in/sphericalSweep.cmd");

      // The sweep ends at the final composition, in a converged state
      TEST_ASSERT(std::abs(sys.mixture().polymer(0).phi() - 0.135) < 1.0E-10);
      TEST_ASSERT(std::abs(sys.mixture().polymer(1).phi() - 0.865) < 1.0E-10);
      checkScft(sys, 1.0E-6);
      sys.computeFreeEnergy();

      // ITERATE at the final composition, from the initial fields of
      // the sweep, must find the same state
      System ref;
      openInputFile("in/spherical3_nr.prm", in);
      ref.readParam(in);
      in.close();
      ref.fileMaster().setInputPrefix(filePrefix());
      ref.fileMaster().setOutputPrefix(filePrefix());
      ref.mixture().polymer(0).setPhi(0.135);
      ref.mixture().polymer(1).setPhi(0.865);
      std::istringstream commands("READ_W in/spherical.w\nITERATE\nFINISH\n");
      runCommands(ref, commands);
      Log::file() << "fHelmholtz (SWEEP, ITERATE) = " << sys.fHelmholtz()
                  << "  " << ref.fHelmholtz() << "\n";
      TEST_ASSERT(std::abs(sys.fHelmholtz() - ref.fHelmholtz()) < 1.0E-6);
      TEST_ASSERT(maxDifferenceC(sys, ref) < 1.0E-5);
   }

};

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testConstructor)
TEST_ADD(SystemTest, testReadParameters)
TEST_ADD(SystemTest, testSolveMdePlanar)
TEST_ADD(SystemTest, testSolveMdeSpherical)
TEST_ADD(SystemTest, testIteratorPlanarNr1)
TEST_ADD(SystemTest, testIteratorPlanarNr2)
TEST_ADD(SystemTest, testIteratorPlanarNr3)
TEST_ADD(SystemTest, testIteratorSphericalNr1)
TEST_ADD(SystemTest, testIteratorSphericalNr3)
TEST_ADD(SystemTest, testIteratorPlanarAm1)
TEST_ADD(SystemTest, testIteratorPlanarAm3)
TEST_ADD(SystemTest, testSweepSpherical)
TEST_END(SystemTest)

#endif
