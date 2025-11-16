/*! io_geotiff.h
 *  This reads the data from a TIFF file using the geotiff library.
 */
#ifndef _IO_GEOTIFF_H_
#define _IO_GEOTIFF_H_

#include <utility>
#include <memory>
#include <boost/array.hpp>
#include "raster.hpp"

namespace raster_stats {
    boost::array<size_t,2> tiff_dimensions(const char* filename);
    void tiff_data_format(const char* filename);
    std::shared_ptr<landscape_t> read_tiff(const char* filename);
    std::shared_ptr<landscape_t> resize_replicate(
                                       std::shared_ptr<landscape_t> praster,
                                       boost::array<landscape_t::size_type,2> ns);
    std::shared_ptr<landscape_t> multi_value(boost::array<landscape_t::size_type,2> ns,
                                             boost::array<landscape_t::value_type,2> vals);
}


#endif // _IO_GEOTIFF_H_
