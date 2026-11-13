#ifndef RP_TESTS_COMPRESSED_INITIAL_STATE_H
#define RP_TESTS_COMPRESSED_INITIAL_STATE_H

#include <rp/fts/compressor/Compressor.h>

namespace Pscf {
namespace Rp {

// Failed-move/step tests intentionally allow only one compressor iteration
// per call. Prepare a valid starting saddle with repeated calls to that same
// compressor, preserving its tolerance, iteration limit, input field files
// and physical parameters. A failed startup may no longer enter the step loop.
template <class SystemT, class SimulatorT>
void prepareCompressedInitialState(SystemT& system, SimulatorT& simulator)
{
   system.compute();
   int error = 1;
   for (int i = 0; i < 1000 && error; ++i) {
      error = simulator.compressor().compress();
   }
   UTIL_CHECK(error == 0);
}

}
}
#endif
