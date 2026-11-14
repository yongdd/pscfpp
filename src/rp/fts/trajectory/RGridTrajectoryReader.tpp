#ifndef RP_RGRID_TRAJECTORY_READER_TPP
#define RP_RGRID_TRAJECTORY_READER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "RGridTrajectoryReader.h"

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <pscf/crystal/UnitCell.h>
#include <util/misc/ioUtil.h>
#include <util/misc/FileMaster.h>

namespace Pscf {
namespace Rp {

   template <int D, class T> class FieldIo;

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   RGridTrajectoryReader<D,T>::RGridTrajectoryReader(
                                         System<D,T>& system)
    : TrajectoryReader<D,T>(system),
      isAllocated_(false),
      hasFrameCell_(false)
   {}

   /*
   * Allocate required memory.
   */
   template <int D, class T>
   void RGridTrajectoryReader<D,T>::allocate()
   {
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer > 0);
      meshDimensions_ = system().domain().mesh().dimensions();
      if (!isAllocated_){
         wField_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            wField_[i].allocate(meshDimensions_);
         }
         isAllocated_ = true;
      }
   }

   /*
   * Open file and allocate memory if necessary.
   */
   template <int D, class T>
   void RGridTrajectoryReader<D,T>::open(std::string filename)
   {
      // The -i input prefix applies to the trajectory file, as to other
      // input files. A name that already includes the prefix (as earlier
      // versions required) is still accepted: it is opened as given when
      // no file exists at the prefixed path. Both paths are opened as
      // plain streams: a failed FileMaster open would write an error
      // message to the log even when the second path is found.
      std::string path = system().fileMaster().inputFilePath(filename);
      inputfile_.open(path.c_str());
      if (inputfile_.fail() && path != filename) {
         inputfile_.clear();
         inputfile_.open(filename.c_str());
      }
      if (inputfile_.fail()) {
         std::string message = "Error opening trajectory file. Filename: ";
         message += path;
         if (path != filename) {
            message += " (or, as given, " + filename + ")";
         }
         UTIL_THROW(message.c_str());
      }
      allocate();
   }

   /*
   * Read trajectory file header section.
   */
   template <int D, class T>
   void RGridTrajectoryReader<D,T>::readHeader()
   {
      int nMonomer = system().mixture().nMonomer();
      FieldIo<D,T> const & fieldIo = system().domain().fieldIo();
      UnitCell<D> tmpUnitCell;
      bool hasSymmetry;
      fieldIo.readFieldHeader(inputfile_, nMonomer, tmpUnitCell,
                              hasSymmetry);
      // readFieldHeader returns the N_monomer of the file
      if (nMonomer != system().mixture().nMonomer()) {
         Log::file() << "N_monomer in trajectory file = " << nMonomer
                     << ", in the system = "
                     << system().mixture().nMonomer() << "\n";
         UTIL_THROW("Trajectory file has a different number of monomers");
      }
      system().setUnitCell(tmpUnitCell);
      headerParameters_ = tmpUnitCell.parameters();
      hasFrameCell_ = false;
      Log::file() << "Read Header" << "\n";
   }

   /*
   * Read frame, return false if end-of-file.
   */
   template <int D, class T>
   bool RGridTrajectoryReader<D,T>::readFrame()
   {
      // Preconditions
      if (!isAllocated_) {
         UTIL_THROW("Memory not allocated in RGridTrajectoryReader");
      }

      bool notEnd;
      std::stringstream line;

      // Attempt to read first line, check for end of file
      notEnd = getNextLine(inputfile_, line);
      if (!notEnd) {
         return false;
      }

      // Read line containing time step
      checkString(line, "i");
      checkString(line, "=");
      #if 0
      std::string value;
      line >> value;
      int step;
      step = std::stoi(value);
      Log::file() << "step " << step << "\n";
      #endif

      // Read an optional "cell_param" line (written by TrajectoryWriter
      // when the unit cell of a frame differs from the header cell, e.g.,
      // under a cell_param ramp), then the mesh dimensions line.
      notEnd = getNextLine(inputfile_, line);
      UTIL_CHECK(notEnd);
      std::string label;
      line >> label;
      if (label == "cell_param") {
         int nParameter = system().domain().unitCell().nParameter();
         FSArray<double, 6> parameters;
         double value;
         for (int i = 0; i < nParameter; ++i) {
            line >> value;
            UTIL_CHECK(!line.fail());
            parameters.append(value);
         }
         system().setUnitCell(parameters);
         hasFrameCell_ = true;
         notEnd = getNextLine(inputfile_, line);
         UTIL_CHECK(notEnd);
         line >> label;
      } else if (hasFrameCell_) {
         // A frame without the line uses the header cell
         system().setUnitCell(headerParameters_);
         hasFrameCell_ = false;
      }
      if (label != "mesh") {
         Log::file() << "Expected \"mesh\", found \"" << label << "\"\n";
         UTIL_THROW("Bad frame in r-grid trajectory file");
      }

      // Read the mesh dimensions of the frame and compare them with the
      // system mesh (a permuted mesh of the same size would otherwise
      // be read point by point in the wrong order)
      notEnd = getNextLine(inputfile_, line);
      UTIL_CHECK(notEnd);
      IntVec<D> frameDimensions;
      for (int i = 0; i < D; ++i) {
         line >> frameDimensions[i];
      }
      if (line.fail() || frameDimensions != meshDimensions_) {
         Log::file() << "Mesh dimensions of frame: " << line.str()
                     << "\nMesh dimensions of system: " << meshDimensions_
                     << "\n";
         UTIL_THROW("Trajectory frame has other mesh dimensions");
      }

      // Read a w-field configuration in r-grid format
      int nMonomer = system().mixture().nMonomer();
      FieldIo<D,T> const & fieldIo = system().domain().fieldIo();
      fieldIo.readFieldsRGridData(inputfile_, wField_, nMonomer);

      // Update system r-grid w fields
      system().w().setRGrid(wField_);

      return true;
   }

   /*
   * Close trajectory file.
   */
   template <int D, class T>
   void RGridTrajectoryReader<D,T>::close()
   {  inputfile_.close(); }

}
}
#endif
