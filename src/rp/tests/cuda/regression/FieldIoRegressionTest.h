#ifndef RPG_FIELD_IO_REGRESSION_TEST_H
#define RPG_FIELD_IO_REGRESSION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>

#include <pscf/field/cuda/RField.h>
#include <pscf/field/cuda/RFieldDft.h>
#include <pscf/field/cuda/RFieldComparison.h>
#include <pscf/crystal/BFieldComparison.h>
#include <pscf/crystal/Basis.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/backend/cuda/HostRecvArray.h>
#include <pscf/backend/cuda/HostSendArray.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/DArray.h>
#include <util/misc/FileMaster.h>
#include <util/misc/Exception.h>
#include <util/tests/LogFileUnitTest.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Regression tests for FieldIo: reading basis files written on another
* mesh or for another unit cell, and file format conversions in systems
* without a space group.
*/
class FieldIoRegressionTest : public LogFileUnitTest
{

   FileMaster fileMaster_;
   int nMonomer_;

public:

   void setUp()
   {
      setVerbose(0);
      nMonomer_ = 2;
      fileMaster_.setRootPrefix(UnitTest::filePrefix());
      openLogFile("out/FieldIo_regression.log");
   }

   /*
   * Open and read file header to initialize Rp::Domain<D,CUT> system.
   */
   template <int D>
   void readHeader(std::string filename, Rp::Domain<D,CUT>& domain)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.readRGridFieldHeader(in, nMonomer_);
      in.close();
   }

   // Allocate an array of fields in symmetry adapated format
   void allocateFields(int nMonomer, int nStar,
                            DArray< DArray<double> >& fields)
   {
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         fields[i].allocate(nStar);
      }
   }

   // Allocate an array of r-grid fields
   template <int D>
   void allocateFields(int nMonomer, IntVec<D> dimensions,
                            DArray< RField<D,CUT> >& fields)
   {
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         fields[i].allocate(dimensions);
      }
   }

   // Allocate an array of k-grid fields
   template <int D>
   void allocateFields(int nMonomer, IntVec<D> dimensions,
                            DArray< RFieldDft<D,CUT> >& fields)
   {
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         fields[i].allocate(dimensions);
      }
   }

   template <int D>
   void readFields(std::string filename, Rp::Domain<D,CUT>& domain,
                   DArray< RField<D,CUT> >& fields)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.fieldIo().readFieldsRGrid(in, fields, domain.unitCell());
      in.close();
   }

   /*
   * Initialize a Domain from field file header text, which ends with
   * the mesh dimensions.
   */
   template <int D>
   void readHeaderString(std::string const & header,
                         Rp::Domain<D,CUT>& domain)
   {
      std::istringstream in(header);
      int nMonomer;
      domain.readRGridFieldHeader(in, nMonomer);
   }

   /*
   * Read basis fields from a file, without letting an exception escape.
   *
   * The unit cell of the file header is read into a separate object.
   * Returns false if the reader threw an Exception.
   */
   template <int D>
   bool tryReadFieldsBasis(std::string filename,
                           Rp::Domain<D,CUT>& domain,
                           DArray< DArray<double> >& fields)
   {
      UnitCell<D> unitCell;
      std::ifstream in;
      openInputFile(filename, in);
      bool success = true;
      try {
         domain.fieldIo().readFieldsBasis(in, fields, unitCell);
      } catch (Exception& e) {
         success = false;
      }
      in.close();
      return success;
   }

   /*
   * Compare basis fields defined on two different meshes in k-space.
   *
   * The fields are compared at every wavevector that both meshes
   * represent without aliasing, i.e., with |k_i| < N_i/2 for the
   * smaller of the two mesh dimensions N_i in each direction. Returns
   * the maximum absolute difference of the DFT coefficients.
   */
   template <int D>
   double compareCommonWaves(Rp::Domain<D,CUT>& domain1,
                             DArray< DArray<double> > const & bf1,
                             Rp::Domain<D,CUT>& domain2,
                             DArray< DArray<double> > const & bf2)
   {
      int nMonomer = bf1.capacity();
      UTIL_CHECK(bf2.capacity() == nMonomer);
      IntVec<D> n1 = domain1.mesh().dimensions();
      IntVec<D> n2 = domain2.mesh().dimensions();

      DArray< RFieldDft<D,CUT> > kf1, kf2;
      allocateFields(nMonomer, n1, kf1);
      allocateFields(nMonomer, n2, kf2);
      domain1.fieldIo().convertBasisToKGrid(bf1, kf1);
      domain2.fieldIo().convertBasisToKGrid(bf2, kf2);

      // Copy the DFT coefficients to host memory
      DArray< HostRecvArray<cudaComplex,CUT> > h1, h2;
      h1.allocate(nMonomer);
      h2.allocate(nMonomer);
      for (int j = 0; j < nMonomer; ++j) {
         h1[j] = kf1[j];
         h2[j] = kf2[j];
      }

      Mesh<D> dftMesh2(kf2[0].dftDimensions());
      MeshIterator<D> iter(kf1[0].dftDimensions());
      IntVec<D> p, q;
      int i, j, k, nMin, rank1, rank2;
      bool isCommon;
      double diff;
      double maxDiff = 0.0;
      int nCommon = 0;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         p = iter.position();
         isCommon = true;
         for (i = 0; i < D; ++i) {
            k = p[i];
            if (2*k > n1[i]) k -= n1[i];  // signed wavevector component
            nMin = (n1[i] < n2[i]) ? n1[i] : n2[i];
            if (2*std::abs(k) >= nMin) isCommon = false;
            q[i] = (k >= 0) ? k : k + n2[i];
         }
         if (!isCommon) continue;
         ++nCommon;
         rank1 = iter.rank();
         rank2 = dftMesh2.rank(q);
         for (j = 0; j < nMonomer; ++j) {
            diff = std::abs(h1[j][rank1].x - h2[j][rank2].x);
            if (diff > maxDiff) maxDiff = diff;
            diff = std::abs(h1[j][rank1].y - h2[j][rank2].y);
            if (diff > maxDiff) maxDiff = diff;
         }
      }
      UTIL_CHECK(nCommon > 1);
      return maxDiff;
   }

   /*
   * Does a text file contain a string?
   */
   bool fileContains(std::string filename, std::string const & text)
   {
      std::ifstream in;
      openInputFile(filename, in);
      std::string line;
      bool found = false;
      while (std::getline(in, line)) {
         if (line.find(text) != std::string::npos) {
            found = true;
         }
      }
      in.close();
      return found;
   }

   /*
   * Read a basis file written on a coarser mesh, for a group without
   * inversion symmetry (I_41_3_2, file mesh 36^3, read on mesh 48^3).
   *
   * Waves with a component 18 = 36/2 belong to closed stars on the
   * mesh of the file, but to pairs of open stars on the finer mesh.
   */
   void testReadBasisCoarseToFine_altG()
   {
      printMethod(TEST_FUNC);
      nMonomer_ = 3;

      // Domain of the file (mesh 36^3)
      Rp::Domain<3,CUT> coarse;
      coarse.setFileMaster(fileMaster_);
      readHeader("../field/in/w_altG.rf", coarse);
      DArray< DArray<double> > bf_c;
      allocateFields(nMonomer_, coarse.basis().nBasis(), bf_c);
      TEST_ASSERT(tryReadFieldsBasis("../field/in/w_altG.bf", coarse, bf_c));

      // Domain with a finer mesh (48^3)
      Rp::Domain<3,CUT> fine;
      fine.setFileMaster(fileMaster_);
      readHeaderString("format 1 0\n dim 3\n crystal_system cubic\n"
                       " N_cell_param 1\n cell_param 3.9\n"
                       " group_name I_41_3_2\n N_monomer 3\n"
                       " mesh 48 48 48\n", fine);
      TEST_ASSERT(fine.basis().nBasis() > coarse.basis().nBasis());

      // Read the coarse-mesh file on the fine mesh
      DArray< DArray<double> > bf_f;
      allocateFields(nMonomer_, fine.basis().nBasis(), bf_f);
      bool success = tryReadFieldsBasis("../field/in/w_altG.bf", fine, bf_f);
      if (verbose() > 0) {
         std::cout << "\nRead succeeded = " << success;
      }
      TEST_ASSERT(success);

      // Compare all waves that the coarse mesh represents exactly
      double maxDiff = compareCommonWaves(coarse, bf_c, fine, bf_f);
      if (verbose() > 0) {
         std::cout << "\nMax k-space difference = " << maxDiff << "\n";
      }
      TEST_ASSERT(maxDiff < 1.0E-10);
   }

   /*
   * Read a basis file written on a finer mesh (double gyroid, group
   * I_a_-3_d, file mesh 32^3, read on mesh 24^3). Some stars of the
   * coarse mesh that contain a wave with a component 12 = 24/2 are
   * cancelled there, but not on the mesh of the file.
   */
   void testReadBasisFineToCoarse_gyr()
   {
      printMethod(TEST_FUNC);

      // Domain of the file (mesh 32^3)
      Rp::Domain<3,CUT> fine;
      fine.setFileMaster(fileMaster_);
      readHeaderString("format 1 0\n dim 3\n crystal_system cubic\n"
                       " N_cell_param 1\n cell_param 3.5985168166\n"
                       " group_name I_a_-3_d\n N_monomer 2\n"
                       " mesh 32 32 32\n", fine);
      DArray< DArray<double> > bf_f;
      allocateFields(nMonomer_, fine.basis().nBasis(), bf_f);
      TEST_ASSERT(tryReadFieldsBasis("in/FieldIo/w_gyr.bf", fine, bf_f));

      // Domain with a coarser mesh (24^3)
      Rp::Domain<3,CUT> coarse;
      coarse.setFileMaster(fileMaster_);
      readHeaderString("format 1 0\n dim 3\n crystal_system cubic\n"
                       " N_cell_param 1\n cell_param 3.5985168166\n"
                       " group_name I_a_-3_d\n N_monomer 2\n"
                       " mesh 24 24 24\n", coarse);
      TEST_ASSERT(coarse.basis().nBasis() < fine.basis().nBasis());

      // Read the fine-mesh file on the coarse mesh
      DArray< DArray<double> > bf_c;
      allocateFields(nMonomer_, coarse.basis().nBasis(), bf_c);
      bool success = tryReadFieldsBasis("in/FieldIo/w_gyr.bf", coarse, bf_c);
      if (verbose() > 0) {
         std::cout << "\nRead succeeded = " << success;
      }
      TEST_ASSERT(success);

      // Compare all waves that the coarse mesh represents exactly
      double maxDiff = compareCommonWaves(fine, bf_f, coarse, bf_c);
      if (verbose() > 0) {
         std::cout << "\nMax k-space difference = " << maxDiff << "\n";
      }
      TEST_ASSERT(maxDiff < 1.0E-10);
   }

   /*
   * Read a basis file written on a finer mesh (bcc, I_m_-3_m, file
   * mesh 32^3, read on mesh 16^3).
   *
   * The file has more lines than the coarse basis has stars, and lines
   * of waves that the coarse mesh cannot represent come before some
   * stars of the coarse basis. Every coarse star must still be read.
   */
   void testReadBasisFineToCoarse_bcc()
   {
      printMethod(TEST_FUNC);

      // Domain of the file (mesh 32^3)
      Rp::Domain<3,CUT> fine;
      fine.setFileMaster(fileMaster_);
      readHeader("../field/in/w_bcc.rf", fine);
      DArray< DArray<double> > bf_f;
      allocateFields(nMonomer_, fine.basis().nBasis(), bf_f);
      TEST_ASSERT(tryReadFieldsBasis("../field/in/w_bcc.bf", fine, bf_f));

      // Domain with a coarser mesh (16^3)
      Rp::Domain<3,CUT> coarse;
      coarse.setFileMaster(fileMaster_);
      readHeaderString("format 1 0\n dim 3\n crystal_system cubic\n"
                       " N_cell_param 1\n cell_param 1.7593562142\n"
                       " group_name I_m_-3_m\n N_monomer 2\n"
                       " mesh 16 16 16\n", coarse);
      TEST_ASSERT(coarse.basis().nBasis() < fine.basis().nBasis());

      // Read the fine-mesh file on the coarse mesh
      DArray< DArray<double> > bf_c;
      allocateFields(nMonomer_, coarse.basis().nBasis(), bf_c);
      TEST_ASSERT(tryReadFieldsBasis("../field/in/w_bcc.bf", coarse, bf_c));

      // Compare all waves that the coarse mesh represents exactly
      double maxDiff = compareCommonWaves(fine, bf_f, coarse, bf_c);
      if (verbose() > 0) {
         std::cout << "\nMax k-space difference = " << maxDiff << "\n";
      }
      TEST_ASSERT(maxDiff < 1.0E-10);
   }

   /*
   * Read a basis file written on a finer mesh (I_41_3_2, file mesh
   * 36^3, read on mesh 24^3). The reader skips the lines of waves that
   * are not minimum images on the coarse mesh, and must report them.
   */
   void testReadBasisFineToCoarse_altG()
   {
      printMethod(TEST_FUNC);
      nMonomer_ = 3;

      // Domain of the file (mesh 36^3)
      Rp::Domain<3,CUT> fine;
      fine.setFileMaster(fileMaster_);
      readHeader("../field/in/w_altG.rf", fine);
      DArray< DArray<double> > bf_f;
      allocateFields(nMonomer_, fine.basis().nBasis(), bf_f);
      TEST_ASSERT(tryReadFieldsBasis("../field/in/w_altG.bf", fine, bf_f));

      // Domain with a coarser mesh (24^3)
      Rp::Domain<3,CUT> coarse;
      coarse.setFileMaster(fileMaster_);
      readHeaderString("format 1 0\n dim 3\n crystal_system cubic\n"
                       " N_cell_param 1\n cell_param 3.9\n"
                       " group_name I_41_3_2\n N_monomer 3\n"
                       " mesh 24 24 24\n", coarse);
      TEST_ASSERT(coarse.basis().nBasis() < fine.basis().nBasis());

      // Read the fine-mesh file on the coarse mesh, with a separate log
      closeLogFile();
      openLogFile("out/FieldIo_readBasisFineToCoarse_altG.log");
      DArray< DArray<double> > bf_c;
      allocateFields(nMonomer_, coarse.basis().nBasis(), bf_c);
      bool success = tryReadFieldsBasis("../field/in/w_altG.bf", coarse, bf_c);
      closeLogFile();
      TEST_ASSERT(success);

      // Compare all waves that the coarse mesh represents exactly
      double maxDiff = compareCommonWaves(fine, bf_f, coarse, bf_c);
      if (verbose() > 0) {
         std::cout << "\nMax k-space difference = " << maxDiff << "\n";
      }
      TEST_ASSERT(maxDiff < 1.0E-10);

      // The skipped lines must be reported in the log
      TEST_ASSERT(fileContains("out/FieldIo_readBasisFineToCoarse_altG.log",
                               "not a minimum image"));
   }

   /*
   * Read back a basis file whose header cell is rounded, at a tie of
   * minimum images (rhombic cell with gamma = 2 pi/3 in full double
   * precision, group p_2, mesh 24 x 24).
   *
   * The file header stores cell parameters with 11 significant digits.
   * The minimum images of some waves differ between the exact cell
   * and the rounded cell of the header.
   */
   void testReadBasisHeaderCell_rhombic()
   {
      printMethod(TEST_FUNC);

      Rp::Domain<2,CUT> domain;
      domain.setFileMaster(fileMaster_);
      readHeaderString("format 1 0\n dim 2\n crystal_system rhombic\n"
                       " N_cell_param 2\n"
                       " cell_param 1.6908669157 2.0943951023931953\n"
                       " group_name p_2\n N_monomer 2\n"
                       " mesh 24 24\n", domain);
      int nBasis = domain.basis().nBasis();
      TEST_ASSERT(nBasis > 0);

      // Fields with a nonzero component for every star
      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, nBasis, bf_0);
      for (int i = 0; i < nMonomer_; ++i) {
         for (int j = 0; j < nBasis; ++j) {
            bf_0[i][j] = (0.5 + std::cos(1.3*j + i))/(1.0 + 0.1*j);
         }
      }

      // Write, then read back with a separate unit cell object
      domain.fieldIo().writeFieldsBasis("out/FieldIo_w_rhombic_tie.bf", bf_0,
                                        domain.unitCell());
      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, nBasis, bf_1);
      TEST_ASSERT(tryReadFieldsBasis("out/FieldIo_w_rhombic_tie.bf",
                                     domain, bf_1));

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      if (verbose() > 0) {
         std::cout << "\nMax difference = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-9);
   }

   /*
   * File conversions between r-grid and k-grid formats, and rescaling
   * of an r-grid file, in a domain without a space group. The input
   * file ../field/in/w_hex.rf has a group_name (p_6_m_m) in its header.
   */
   void testConvertFilesWithoutGroup_hex()
   {
      printMethod(TEST_FUNC);

      // Hexagonal domain with the mesh of ../field/in/w_hex.rf, but no group
      Rp::Domain<2,CUT> domain;
      domain.setFileMaster(fileMaster_);
      readHeaderString("format 1 0\n dim 2\n crystal_system hexagonal\n"
                       " N_cell_param 1\n cell_param 1.6651597308\n"
                       " N_monomer 2\n mesh 32 32\n", domain);
      TEST_ASSERT(!domain.hasGroup());
      domain.fieldIo().setNMonomer(nMonomer_);
      IntVec<2> dimensions = domain.mesh().dimensions();

      DArray< RField<2,CUT> > rf_0;
      allocateFields(nMonomer_, dimensions, rf_0);
      readFields("../field/in/w_hex.rf", domain, rf_0);

      // r-grid -> k-grid -> r-grid
      bool success = true;
      try {
         domain.fieldIo().convertRGridToKGrid("../field/in/w_hex.rf",
                                              "out/FieldIo_w_hex_nogroup.kf");
         domain.fieldIo().convertKGridToRGrid("out/FieldIo_w_hex_nogroup.kf",
                                              "out/FieldIo_w_hex_nogroup.rf");
      } catch (Exception& e) {
         success = false;
      }
      TEST_ASSERT(success);

      DArray< RField<2,CUT> > rf_1;
      allocateFields(nMonomer_, dimensions, rf_1);
      readFields("out/FieldIo_w_hex_nogroup.rf", domain, rf_1);
      RFieldComparison<2,CUT> comparison;
      comparison.compare(rf_0, rf_1);
      if (verbose() > 0) {
         std::cout << "\nMax difference = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

      // Rescale the r-grid file (its header names a group)
      success = true;
      try {
         domain.fieldIo().scaleFieldsRGrid("../field/in/w_hex.rf",
                                           "out/FieldIo_w_hex_nogroup_x2.rf",
                                           2.0);
      } catch (Exception& e) {
         success = false;
      }
      TEST_ASSERT(success);
      readFields("out/FieldIo_w_hex_nogroup_x2.rf", domain, rf_1);
      domain.fieldIo().scaleFieldsRGrid(rf_0, 2.0);
      comparison.compare(rf_0, rf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

      // Output files of a domain without a group carry no group name
      TEST_ASSERT(!fileContains("out/FieldIo_w_hex_nogroup.kf", "group_name"));
      TEST_ASSERT(!fileContains("out/FieldIo_w_hex_nogroup.rf", "group_name"));
      TEST_ASSERT(!fileContains("out/FieldIo_w_hex_nogroup_x2.rf",
                                "group_name"));
   }

   /*
   * Conversion of an r-grid file without a group_name to k-grid format
   * in a domain with a group: the output must not claim the group.
   */
   void testConvertRGridToKGridAsymmetric_hex()
   {
      printMethod(TEST_FUNC);

      Rp::Domain<2,CUT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("../field/in/w_hex.rf", domain);
      TEST_ASSERT(domain.hasGroup());
      domain.fieldIo().setNMonomer(nMonomer_);
      IntVec<2> dimensions = domain.mesh().dimensions();

      // Write an r-grid field that has no group_name in its header
      DArray< RField<2,CUT> > rf_0;
      allocateFields(nMonomer_, dimensions, rf_0);
      readFields("../field/in/w_hex.rf", domain, rf_0);
      int n = domain.mesh().size();
      HostRecvArray<cudaReal,CUT> h_rf(n);
      HostSendArray<cudaReal,CUT> h_rf2(n);
      for (int i = 0; i < nMonomer_; ++i) {
         h_rf = rf_0[i];
         for (int j = 0; j < n; ++j) {
            h_rf2[j] = h_rf[j] + 0.1*std::sin(0.37*j + i);
         }
         rf_0[i] = h_rf2;
      }
      domain.fieldIo().writeFieldsRGrid("out/FieldIo_w_hex_asym.rf", rf_0,
                                        domain.unitCell(), false);
      TEST_ASSERT(!fileContains("out/FieldIo_w_hex_asym.rf", "group_name"));

      domain.fieldIo().convertRGridToKGrid("out/FieldIo_w_hex_asym.rf",
                                           "out/FieldIo_w_hex_asym.kf");
      TEST_ASSERT(!fileContains("out/FieldIo_w_hex_asym.kf", "group_name"));
   }

};

TEST_BEGIN(FieldIoRegressionTest)
TEST_ADD(FieldIoRegressionTest, testReadBasisCoarseToFine_altG)
TEST_ADD(FieldIoRegressionTest, testReadBasisFineToCoarse_gyr)
TEST_ADD(FieldIoRegressionTest, testReadBasisFineToCoarse_bcc)
TEST_ADD(FieldIoRegressionTest, testReadBasisFineToCoarse_altG)
TEST_ADD(FieldIoRegressionTest, testReadBasisHeaderCell_rhombic)
TEST_ADD(FieldIoRegressionTest, testConvertFilesWithoutGroup_hex)
TEST_ADD(FieldIoRegressionTest, testConvertRGridToKGridAsymmetric_hex)
TEST_END(FieldIoRegressionTest)

#endif
