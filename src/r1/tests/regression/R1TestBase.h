#ifndef R1_TEST_BASE_H
#define R1_TEST_BASE_H

#include <test/UnitTest.h>

#include <r1/system/System.h>
#include <util/param/ParamComponent.h>
#include <util/misc/Log.h>
#include <util/global.h>

#include <fstream>
#include <sstream>
#include <string>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

/*
* Base class of the r1 regression tests: log file and System helpers.
*/
class R1TestBase : public UnitTest
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
   * Read a System parameter file and set the file prefixes of the
   * System to the unit test file prefix.
   */
   void readSystem(System& sys, char const * filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      sys.readParam(in);
      in.close();
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());
   }

   /*
   * Execute the commands in a string.
   */
   void runCommands(System& sys, std::string const & commands)
   {
      std::istringstream in(commands);
      sys.readCommands(in);
   }

   /*
   * Return the contents of a file.
   */
   std::string fileContents(char const * filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      std::stringstream buffer;
      buffer << in.rdbuf();
      return buffer.str();
   }

   /*
   * Return the value of a parameter in parameter file format text.
   */
   double paramValue(std::string const & text, std::string const & label)
   {
      std::size_t pos = text.find(label);
      TEST_ASSERT(pos != std::string::npos);
      std::istringstream in(text.substr(pos + label.size()));
      double value;
      in >> value;
      TEST_ASSERT(!in.fail());
      return value;
   }

};

#endif
