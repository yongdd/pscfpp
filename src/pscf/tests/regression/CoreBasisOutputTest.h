#ifndef PSCF_CORE_BASIS_OUTPUT_TEST_H
#define PSCF_CORE_BASIS_OUTPUT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/crystal/Basis.h>
#include <pscf/crystal/SpaceGroup.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>

#include <util/containers/FSArray.h>

#include <sstream>
#include <string>
#include <vector>

using namespace Util;
using namespace Pscf;

/*
* Basis::outputWaves and Basis::outputStars (WRITE_WAVES, WRITE_STARS):
* every integer column must stay separated from its neighbors, also
* when a value fills its field width.
*/
class CoreBasisOutputTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   /*
   * Split a line into integer tokens (fails if a token is not an int).
   */
   std::vector<long> intTokens(std::string const & line, int nInt)
   {
      std::istringstream ss(line);
      std::vector<long> tokens;
      std::string token;
      while (ss >> token) {
         if ((int)tokens.size() < nInt) {
            std::size_t pos = 0;
            long value = std::stol(token, &pos);
            TEST_ASSERT(pos == token.size());
            tokens.push_back(value);
         } else {
            tokens.push_back(0);
         }
      }
      return tokens;
   }

   /*
   * Make a basis for a 1D mesh with group P_-1.
   */
   void makeBasis1D(int n, Basis<1>& basis, Mesh<1>& mesh,
                    UnitCell<1>& cell)
   {
      FSArray<double, 6> p;
      p.append(1.6*n/1000.0);
      cell.set(UnitCell<1>::Lamellar, p);
      IntVec<1> d;
      d[0] = n;
      mesh.setDimensions(d);
      SpaceGroup<1> group;
      readGroup("P_-1", group);
      basis.makeBasis(mesh, cell, group);
      TEST_ASSERT(basis.isValid());
   }

   /*
   * Check the waves written by outputWaves: one line per wave, with
   * k, i, the D minimum-image indices and starId as separate columns.
   */
   template <int D>
   void checkWaves(Basis<D> const & basis)
   {
      std::ostringstream out;
      basis.outputWaves(out);
      std::istringstream in(out.str());
      std::string line;
      std::getline(in, line);
      TEST_ASSERT(line == "N_wave");
      std::getline(in, line);
      int nWave = std::stoi(line);
      TEST_ASSERT(nWave == basis.nBasisWave());
      int nCol = 2 + D + 1 + 2;
      int k = 0;
      for (int i = 0; i < basis.nWave(); ++i) {
         if (basis.star(basis.wave(i).starId).cancel) continue;
         std::getline(in, line);
         std::vector<long> t = intTokens(line, 3 + D);
         TEST_ASSERT((int)t.size() == nCol);
         TEST_ASSERT(t[0] == k);
         TEST_ASSERT(t[1] == i);
         for (int j = 0; j < D; ++j) {
            TEST_ASSERT(t[2+j] == basis.wave(i).indicesMin[j]);
         }
         TEST_ASSERT(t[2+D] == basis.wave(i).starId);
         ++k;
      }
      TEST_ASSERT(k == nWave);
   }

   /*
   * Check the stars written by outputStars: basisId, starId, size,
   * beginId, endId, invertFlag and the D indices of waveMin.
   */
   template <int D>
   void checkStars(Basis<D> const & basis)
   {
      std::ostringstream out;
      basis.outputStars(out);
      std::istringstream in(out.str());
      std::string line;
      std::getline(in, line);
      TEST_ASSERT(line == "N_basis");
      std::getline(in, line);
      int nBasis = std::stoi(line);
      TEST_ASSERT(nBasis == basis.nBasis());
      int nCol = 6 + D;
      int ib = 0;
      for (int i = 0; i < basis.nStar(); ++i) {
         typename Basis<D>::Star const & star = basis.star(i);
         if (star.cancel) continue;
         std::getline(in, line);
         std::vector<long> t = intTokens(line, nCol);
         TEST_ASSERT((int)t.size() == nCol);
         TEST_ASSERT(t[0] == ib);
         TEST_ASSERT(t[0] == star.basisId);
         TEST_ASSERT(t[1] == i);
         TEST_ASSERT(t[2] == star.size);
         TEST_ASSERT(t[3] == star.beginId);
         TEST_ASSERT(t[4] == star.endId);
         TEST_ASSERT(t[5] == star.invertFlag);
         for (int j = 0; j < D; ++j) {
            TEST_ASSERT(t[6+j] == star.waveMin[j]);
         }
         ++ib;
      }
      TEST_ASSERT(ib == nBasis);
   }

   /*
   * 1D mesh of 2001 points: the wave index -1000 fills its field.
   */
   void testWaves1D()
   {
      printMethod(TEST_FUNC);
      Basis<1> basis;
      Mesh<1> mesh;
      UnitCell<1> cell;
      makeBasis1D(2001, basis, mesh, cell);
      checkWaves<1>(basis);
      checkStars<1>(basis);
   }

   /*
   * 1D mesh of 200002 points: star and basis ids of 100000 or more,
   * and wave indices of 100000 or more.
   */
   void testStars1DLarge()
   {
      printMethod(TEST_FUNC);
      Basis<1> basis;
      Mesh<1> mesh;
      UnitCell<1> cell;
      makeBasis1D(200002, basis, mesh, cell);
      TEST_ASSERT(basis.nStar() > 100000);
      checkStars<1>(basis);
      checkWaves<1>(basis);
   }

};

TEST_BEGIN(CoreBasisOutputTest)
TEST_ADD(CoreBasisOutputTest, testWaves1D)
TEST_ADD(CoreBasisOutputTest, testStars1DLarge)
TEST_END(CoreBasisOutputTest)

#endif
