/*! raster_times.h
 *  Makes a list of functions that do clustering so that Python
 *  can pick which one to run. Intended to make it easier to
 *  wrap and run sets of candidate functions.
 */

#ifndef _RASTER_TIMES_H_
#define _RASTER_TIMES_H_ 1

#include "timing_harness.hpp"

namespace raster_stats {
std::vector<timing_harness> raster_times();
}

#endif
