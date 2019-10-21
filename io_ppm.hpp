/*! This file writes the PPM file format, an easy text-only format.
 */
#ifndef _IO_PPM_H_
#define _IO_PPM_H_ 1


#include "raster.hpp"

namespace raster_stats {
	void write_ppm(const landscape_t& raster, const char* filename);
}

#endif // _IO_PPM_H_
