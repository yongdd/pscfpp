#ifndef PSCF_CORE_BASIS_DATA_TEST_H
#define PSCF_CORE_BASIS_DATA_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/fieldIo/rFieldIo.h>
#include <pscf/crystal/Basis.h>
#include <pscf/crystal/SpaceGroup.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>

#include <util/containers/DArray.h>
#include <util/containers/FSArray.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

using namespace Util;
using namespace Pscf;

/*
* Data section of basis field files (writeBasisData, readBasisData):
* the Miller indices and the star size must be separate columns, also
* when an index fills its field width, and lines whose values fit must
* be written as before.
*/
class CoreBasisDataTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   /*
   * Make a basis with group p_1 for a rectangular 2D mesh.
   */
   void makeBasis2D(IntVec<2> const & d, Basis<2>& basis,
                    Mesh<2>& mesh, UnitCell<2>& cell)
   {
      FSArray<double, 6> p;
      p.append(0.5);
      p.append(0.5*d[1]/d[0]);
      cell.set(UnitCell<2>::Rectangular, p);
      mesh.setDimensions(d);
      SpaceGroup<2> group;
      readGroup("p_1", group);
      basis.makeBasis(mesh, cell, group);
      TEST_ASSERT(basis.isValid());
   }

   /*
   * Line of the data section in the format written before the fix.
   */
   std::string oldLine(DArray< DArray<double> > const & fields, int ib,
                       Basis<2>::Star const & star)
   {
      std::ostringstream out;
      for (int j = 0; j < fields.capacity(); ++j) {
         out << Dbl(fields[j][ib], 20, 10);
      }
      out << "   ";
      for (int j = 0; j < 2; ++j) {
         out << Int(star.waveMin[j], 5);
      }
      out << Int(star.size, 5);
      return out.str();
   }

   /*
   * Write fields of a 2D basis, check every line and read them back.
   */
   void checkWriteRead(IntVec<2> const & d, int& nWide)
   {
      Basis<2> basis;
      Mesh<2> mesh;
      UnitCell<2> cell;
      makeBasis2D(d, basis, mesh, cell);
      const int nBasis = basis.nBasis();
      const int nMonomer = 2;

      DArray< DArray<double> > fields;
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         fields[i].allocate(nBasis);
         for (int ib = 0; ib < nBasis; ++ib) {
            fields[i][ib] = (i + 1)*std::cos(0.1*ib);
         }
      }

      std::ostringstream out;
      writeBasisData(out, fields, basis);

      // Check the columns of every line
      std::istringstream in(out.str());
      std::string line;
      nWide = 0;
      int ib = 0;
      for (int i = 0; i < basis.nStar(); ++i) {
         Basis<2>::Star const & star = basis.star(i);
         if (star.cancel) continue;
         TEST_ASSERT(std::getline(in, line));

         // Split into words: nMonomer values, 2 indices and the size
         std::istringstream ls(line);
         std::vector<std::string> words;
         std::string word;
         while (ls >> word) {
            words.push_back(word);
         }
         TEST_ASSERT((int)words.size() == nMonomer + 3);
         for (int j = 0; j < 2; ++j) {
            std::size_t pos = 0;
            int index = std::stoi(words[nMonomer + j], &pos);
            TEST_ASSERT(pos == words[nMonomer + j].size());
            TEST_ASSERT(index == star.waveMin[j]);
         }
         TEST_ASSERT(std::stoi(words[nMonomer + 2]) == star.size);

         // Unchanged where every integer fits in 4 characters
         bool fits = (star.size < 1000);
         for (int j = 0; j < 2; ++j) {
            int index = star.waveMin[j];
            if (index > 9999 || index < -999) fits = false;
         }
         if (fits) {
            TEST_ASSERT(line == oldLine(fields, ib, star));
         } else {
            ++nWide;
         }
         ++ib;
      }
      TEST_ASSERT(ib == nBasis);

      // Read back
      DArray< DArray<double> > fieldsIn;
      fieldsIn.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         fieldsIn[i].allocate(nBasis);
      }
      std::istringstream in2(out.str());
      readBasisData(in2, fieldsIn, cell, mesh, basis, nBasis);
      for (int i = 0; i < nMonomer; ++i) {
         for (int ib = 0; ib < nBasis; ++ib) {
            double diff = std::abs(fieldsIn[i][ib] - fields[i][ib]);
            TEST_ASSERT(diff < 1.0E-9);
         }
      }
   }

   /*
   * Mesh of 2 x 2400 points: second indices of about -1200 to 1200,
   * of which those of -1000 or less fill 5 characters.
   */
   void testNegativeIndices()
   {
      printMethod(TEST_FUNC);
      IntVec<2> d;
      d[0] = 2;
      d[1] = 2400;
      int nWide;
      checkWriteRead(d, nWide);
      TEST_ASSERT(nWide >= 2*199);
   }

   /*
   * Mesh of 2 x 20002 points: second indices of about -10000 to
   * 10000, of which 10000 or more and -1000 or less fill 5 characters
   * or more.
   */
   void testLargeIndices()
   {
      printMethod(TEST_FUNC);
      IntVec<2> d;
      d[0] = 2;
      d[1] = 20002;
      int nWide;
      checkWriteRead(d, nWide);
      TEST_ASSERT(nWide > 0);
   }

};

TEST_BEGIN(CoreBasisDataTest)
TEST_ADD(CoreBasisDataTest, testNegativeIndices)
TEST_ADD(CoreBasisDataTest, testLargeIndices)
TEST_END(CoreBasisDataTest)

#endif
