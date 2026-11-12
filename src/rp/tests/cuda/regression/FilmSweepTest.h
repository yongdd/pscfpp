#ifndef RPG_FILM_SWEEP_TEST_H
#define RPG_FILM_SWEEP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/Mask.h>
#include <rp/environment/MixAndMatchEnv.h>
#include <rp/environment/FieldGenerator.h>
#include <rp/environment/FilmFieldGenExt.h>

#include <pscf/crystal/BFieldComparison.h>

#include <util/containers/DArray.h>
#include <util/misc/Exception.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <string>

using namespace Util;
using namespace Pscf;

/*
* Regression tests for sweeps of the wall parameters of thin films.
*/
class FilmSweepTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Read a parameter file to initialize a System.
   */
   template <int D>
   void setupSystem(Rp::System<D,CUT>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

   /*
   * Read a parameter file section to create a FilmFieldGenExt object.
   */
   template <int D>
   void createFilmFieldGenExt(Rp::FilmFieldGenExt<D,CUT>& ext,
                              std::string fname)
   {
      std::ifstream in;
      openInputFile(fname, in);
      ext.readParameters(in);
      in.close();
   }

   /*
   * Sweep the wall parameters chi_bottom and chi_top of a 1D film with
   * a rigid unit cell. At the end of the sweep, the wall fields h must
   * be those generated for the final parameter values.
   */
   void testSweepWallChi1D()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmSweepWallChi1D.log");

      Rp::System<1,CUT> system;
      setupSystem(system, "in/Film/sweep1D_param");
      system.w().readBasis("in/Film/sweep1D_w.bf");
      system.sweep();

      // Check the final values of the swept parameters
      Rp::MixAndMatchEnv& env
            = dynamic_cast<Rp::MixAndMatchEnv&>(system.environment());
      DArray<int> ids;
      ids.allocate(1);
      ids[0] = 0;
      Rp::FieldGenerator const & ext = env.fieldGenerator2();
      TEST_ASSERT(eq(ext.getParameter("chi_bottom", ids), 3.0));
      TEST_ASSERT(eq(ext.getParameter("chi_top", ids), 3.0));

      // Wall fields at the end of the sweep
      TEST_ASSERT(system.h().hasData());
      DArray< DArray<double> > hSweep;
      hSweep = system.h().basis();

      // Regenerate the fields for the final parameters, and compare
      env.reset();
      env.generate();
      BFieldComparison comparison(0);
      comparison.compare(system.h().basis(), hSweep);
      Log::file() << "\nmax change of h on regeneration = "
                  << comparison.maxDiff() << std::endl;
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   /*
   * Monomer ids of the sweep parameters chi_bottom and chi_top must be
   * in range [0, nMonomer).
   */
   void testParameterIdRange()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmParameterIdRange.log");

      Rp::System<1,CUT> system;
      setupSystem(system, "in/Film/sweep1D_param");
      Rp::FilmFieldGenExt<1,CUT> ext(system);
      createFilmFieldGenExt(ext, "in/Film/ext_asym");

      DArray<int> ids;
      ids.allocate(1);
      std::string names[2] = {"chi_bottom", "chi_top"};
      int badIds[2] = {-1, 2};
      bool thrown;
      for (int i = 0; i < 2; ++i) {

         // Valid ids
         ids[0] = 1;
         ext.setParameter(names[i], ids, 3.0);
         TEST_ASSERT(eq(ext.getParameter(names[i], ids), 3.0));

         // Invalid ids
         for (int j = 0; j < 2; ++j) {
            ids[0] = badIds[j];
            thrown = false;
            try {
               ext.getParameter(names[i], ids);
            } catch (Exception& e) {
               thrown = true;
            }
            TEST_ASSERT(thrown);
         }
         ids[0] = 2;
         thrown = false;
         try {
            ext.setParameter(names[i], ids, 1.0);
         } catch (Exception& e) {
            thrown = true;
         }
         TEST_ASSERT(thrown);
      }
   }

};

TEST_BEGIN(FilmSweepTest)
TEST_ADD(FilmSweepTest, testSweepWallChi1D)
TEST_ADD(FilmSweepTest, testParameterIdRange)
TEST_END(FilmSweepTest)

#endif
