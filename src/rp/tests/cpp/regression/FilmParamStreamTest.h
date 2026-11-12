#ifndef RPC_FILM_PARAM_STREAM_TEST_H
#define RPC_FILM_PARAM_STREAM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/environment/MixAndMatchEnv.h>
#include <rp/environment/FieldGenerator.h>
#include <rp/environment/FilmFieldGenMaskBase.h>
#include <rp/environment/FilmFieldGenExtBase.h>

#include <util/param/Label.h>
#include <util/misc/Exception.h>
#include <util/tests/LogFileUnitTest.h>

#include <istream>
#include <sstream>
#include <streambuf>
#include <string>

using namespace Util;
using namespace Pscf;

/*
* Regression tests for reading a FilmEnvironment parameter block from a
* stream that cannot be rewound (a pipe or a FIFO).
*/
class FilmParamStreamTest : public LogFileUnitTest
{

public:

   /*
   * Stream buffer that gives a text and cannot seek, like a pipe.
   */
   class ForwardOnlyStreamBuf : public std::streambuf
   {
   public:

      explicit ForwardOnlyStreamBuf(std::string const & text)
       : text_(text)
      {
         char* begin = &text_[0];
         setg(begin, begin, begin + text_.size());
      }

   private:

      std::string text_;

   };

   void setUp()
   {  setVerbose(0); }

   /*
   * Parameter file of a 1D film with an Iterator.
   *
   * \param fBulk  line with the parameter fBulk, or ""
   */
   std::string filmParam(std::string const & fBulk)
   {
      std::ostringstream param;
      param << "System{\n"
            << "  Mixture{\n"
            << "    nMonomer  2\n"
            << "    monomers[  1.0  1.0  ]\n"
            << "    nPolymer  1\n"
            << "    Polymer{\n"
            << "      type    linear\n"
            << "      nBlock  2\n"
            << "      blocks[\n"
            << "              0  0.40\n"
            << "              1  0.60\n"
            << "      ]\n"
            << "      phi     1.0\n"
            << "    }\n"
            << "    ds   0.01\n"
            << "  }\n"
            << "  Interaction{\n"
            << "    chi(\n"
            << "         1   0   8.0\n"
            << "    )\n"
            << "  }\n"
            << "  Domain{\n"
            << "    mesh           32\n"
            << "    lattice        lamellar\n"
            << "    groupName      P_1\n"
            << "  }\n"
            << "  FilmEnvironment{\n"
            << "    normalVecId        0\n"
            << "    interfaceThickness 0.3\n"
            << "    excludedThickness  0.5\n"
            << fBulk
            << "    chiBottom[  5.0  0.0  ]\n"
            << "    chiTop[  2.0  10.0  ]\n"
            << "  }\n"
            << "  AmIteratorBasis{\n"
            << "    epsilon      1.0e-8\n"
            << "    maxItr       100\n"
            << "    maxHist      25\n"
            << "    isFlexible   0\n"
            << "  }\n"
            << "}\n";
      return param.str();
   }

   /*
   * Read a parameter file from a stream; returns true if an Exception
   * was thrown.
   */
   bool readParamThrows(Rp::System<1,CPT>& system, std::istream& in)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      try {
         system.readParam(in);
      } catch (Exception& e) {
         Label::clear();
         return true;
      }
      return false;
   }

   /*
   * Read the same parameter file from a seekable stream and from a
   * stream that cannot seek, and compare the two systems.
   */
   void checkFilmParam(std::string const & fBulk)
   {
      std::string text = filmParam(fBulk);

      // Seekable stream
      Rp::System<1,CPT> ref;
      std::istringstream refIn(text);
      TEST_ASSERT(!readParamThrows(ref, refIn));

      // Stream that cannot seek
      Rp::System<1,CPT> system;
      ForwardOnlyStreamBuf buf(text);
      std::istream in(&buf);
      TEST_ASSERT(in.tellg() == std::streampos(-1));
      TEST_ASSERT(!readParamThrows(system, in));

      // The blocks after the Environment were read
      TEST_ASSERT(system.hasEnvironment());
      TEST_ASSERT(system.hasIterator());

      // Parameters of the mask and wall field generators
      Rp::MixAndMatchEnv const & env
         = dynamic_cast<Rp::MixAndMatchEnv const &>(system.environment());
      Rp::FilmFieldGenMaskBase<1> const & mask
         = dynamic_cast<Rp::FilmFieldGenMaskBase<1> const &>
                                                   (env.fieldGenerator1());
      Rp::FilmFieldGenExtBase<1> const & ext
         = dynamic_cast<Rp::FilmFieldGenExtBase<1> const &>
                                                   (env.fieldGenerator2());
      TEST_ASSERT(mask.normalVecId() == 0);
      TEST_ASSERT(eq(mask.interfaceThickness(), 0.3));
      TEST_ASSERT(eq(mask.excludedThickness(), 0.5));
      TEST_ASSERT(mask.hasFBulk() == !fBulk.empty());
      TEST_ASSERT(ext.normalVecId() == 0);
      TEST_ASSERT(eq(ext.interfaceThickness(), 0.3));
      TEST_ASSERT(eq(ext.excludedThickness(), 0.5));
      TEST_ASSERT(eq(ext.chiBottom(0), 5.0));
      TEST_ASSERT(eq(ext.chiBottom(1), 0.0));
      TEST_ASSERT(eq(ext.chiTop(0), 2.0));
      TEST_ASSERT(eq(ext.chiTop(1), 10.0));

      // The parameter files written by the two systems are identical
      std::ostringstream refOut, out;
      ref.writeParamNoSweep(refOut);
      system.writeParamNoSweep(out);
      TEST_ASSERT(out.str() == refOut.str());
   }

   /*
   * Without fBulk: the mask generator leaves the label chiBottom in the
   * Label buffer when it looks for the optional fBulk.
   */
   void testNoFBulk()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmParamStreamNoFBulk.log");
      checkFilmParam("");
   }

   /*
   * With fBulk.
   */
   void testFBulk()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmParamStreamFBulk.log");
      checkFilmParam("    fBulk              2.5\n");
   }

};

TEST_BEGIN(FilmParamStreamTest)
TEST_ADD(FilmParamStreamTest, testNoFBulk)
TEST_ADD(FilmParamStreamTest, testFBulk)
TEST_END(FilmParamStreamTest)

#endif
