#include <utility>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <memory>
#include <boost/assert.hpp>
#include <boost/array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "tiffio.h"
#include "xtiffio.h"
//#include "geotiff/xtiffio.h"
#include "io_geotiff.hpp"

using namespace std;
using namespace boost::numeric::ublas;

namespace raster_stats {


boost::array<size_t,2> tiff_dimensions(const char* filename)
{
    uint32 width, height;
    
    TIFF* raster = XTIFFOpen(filename,"r");
    if (NULL == raster) {
        throw std::runtime_error("Could not open TIFF.");
    }
    int read_width = TIFFGetField(raster, TIFFTAG_IMAGEWIDTH, &width);
    assert(read_width == 1);
    int read_height = TIFFGetField(raster, TIFFTAG_IMAGELENGTH, &height);
    assert(read_height == 1);
    if (read_width != 1 || read_height != 1) {
		throw std::runtime_error("Could not read TIFF width and height.");
    }
    
    XTIFFClose(raster);
    
    return {{height,width}};
}


template<class T>
T tiff_get_field(TIFF* raster, ttag_t tag, const std::string& name)
{
    T item;
    int read_item = TIFFGetField(raster, tag, &item);
    //assert(read_item == 1);
    if (1 != read_item) {
		std::cout << name << " could not read" << std::endl;
        //throw std::runtime_error("Could not read tiff tag.");
    } else {
		std::cout << name << " = " << item << std::endl;
	}
    return item;
}



/*! tiff_data_format, pulls out formatting information as found in
 *  http://partners.adobe.com/public/developer/en/tiff/TIFF6.pdf.
 */
void tiff_data_format(const char* filename)
{
    
    TIFF* raster = XTIFFOpen(filename,"r");
    if (NULL == raster) {
        throw std::runtime_error("Could not open TIFF.");
    }
    
	// This could be 3 for RGB, and each would have a different samples per pixel.
	uint16 samples_per_pixel = tiff_get_field<uint16>(raster,TIFFTAG_SAMPLESPERPIXEL,
				"samples per pixel");
	
    uint16 bits_per_sample = tiff_get_field<uint16>(raster,TIFFTAG_BITSPERSAMPLE,
			"bits per sample");
	uint16 bits_per_sample2 = tiff_get_field<uint16>(raster,258,
			"bits per sample");
    
    uint16 data_type = tiff_get_field<uint16>(raster,TIFFTAG_DATATYPE, "data type");
    
	// endianness. default to 1 for little-endian.
    uint16 fill_order = tiff_get_field<uint16>(raster,TIFFTAG_FILLORDER, "fill order");
    
    uint32 image_depth = tiff_get_field<uint32>(raster,TIFFTAG_IMAGEDEPTH, "image depth");
    uint32 strip_byte_counts = tiff_get_field<uint32>(raster,TIFFTAG_STRIPBYTECOUNTS, "strip byte counts");
    
	const uint16 PlanarConfiguration = 284;
	// 1 = packed, 2 = separate planes.
    uint16 planar_configuration = tiff_get_field<uint16>(raster,PlanarConfiguration,
		"packed or planes");

    uint16 min_value = tiff_get_field<uint16>(raster,TIFFTAG_MINSAMPLEVALUE,
		"min sample value");
    uint16 max_value = tiff_get_field<uint16>(raster,TIFFTAG_MAXSAMPLEVALUE,
		"max sample value");
	// default orientation: 0th row is top of image. 0th column is lhs.
    uint16 orientation = tiff_get_field<uint16>(raster,TIFFTAG_ORIENTATION, "orientation");
    uint16 sample_format = tiff_get_field<uint16>(raster,TIFFTAG_SAMPLEFORMAT, "sample format");
    // 1 = Byte, 2 = ASCII (7-bit), 3 = SHORT (2-byte), 4 = LONG (4-byte),
	// 5 = RATIONAL (2 LONGS, a numerator and denominator), 6 = SBYTE (signed byte),
	// 7 = UNDEFINED (1 byte), 8 = SSHORT (signed 2-byte), 9 = SLONG (signed 4-byte),
	// 10 = SRATIONAL (two signed longs), 11 = FLOAT (single precision 4-byte IEEE),
	// 12 = DOUBLE (double precision 8-byte IEEE)
	//TIFFDataType tiff_data_type = data_type;
    //int width_bytes = TIFFDataWidth(tiff_data_type);
    //assert(width_bytes != 0);
    //if (width_bytes == 0) {
    //    throw std::runtime_error("Could not get TIFF data width in bytes.")
    //}
    
    XTIFFClose(raster);
}


std::shared_ptr<landscape_t> read_tiff(const char* filename)
{
    uint32 width=0, height=0;
    TIFF* raster = XTIFFOpen(filename,"r");
	if ( 0 == raster ) {
	 	throw std::runtime_error("Could not open TIFF.");
	}

    int read_width = TIFFGetField(raster, TIFFTAG_IMAGEWIDTH, &width);
    assert(read_width == 1);
    int read_height = TIFFGetField(raster, TIFFTAG_IMAGELENGTH, &height);
    assert(read_height == 1);
    if (read_width != 1 || read_height != 1) {
		throw std::runtime_error("Could not read TIFF width and height.");
    }
	
	// Create a boost matrix large enough to hold the scan lines.
	std::shared_ptr<landscape_t> image(new landscape_t(height,width));
	landscape_t& rimage = *image;
	
	uint32 scanline_size = TIFFScanlineSize(raster);

	// Create a separate buffer to hold a single scan line b/c this can be longer than matrix width.
	uint8* line_buffer = static_cast<uint8*>(_TIFFmalloc(TIFFScanlineSize(raster)));
	
	// Read scan lines and copy the used length into the matrix, from the top left.
	const int pixel_size = 1;
	for (size_t row_idx = 0; row_idx < height; row_idx++) {
		TIFFReadScanline(raster, line_buffer, row_idx);
		//::memcpy(&image->data()+(height-row_idx-1)*width, line_buffer, pixel_size*width);
		for (size_t col_idx=0; col_idx<width; col_idx++) {
			rimage(height-row_idx-1,col_idx) = line_buffer[col_idx];
		}
	}

	_TIFFfree(line_buffer);

    XTIFFClose(raster);
	return image;
}



/*! This creates a new matrix of the given size using copies
 *  of the given matrix. It copies the matrix in blocks using
 *  BLAS functions.
 */
std::shared_ptr<landscape_t> resize_replicate(std::shared_ptr<landscape_t> praster,
                               boost::array<landscape_t::size_type,2> ns)
{
  const landscape_t& raster(*praster);

  std::shared_ptr<landscape_t> pmorph(new landscape_t(ns[0],ns[1]));
  landscape_t& morph(*pmorph);

  const boost::array<landscape_t::size_type,2> os = {{ raster.size1(), raster.size2() }};

  typedef landscape_t::size_type size_type;

  // Copy the data in blocks.
  size_type iblocks = ns[0]/os[0];
  size_type jblocks = ns[1]/os[1];

  for (size_type i=0; i<iblocks; i++) {
    for (size_type j=0; j<jblocks; j++) {
		noalias(project(morph,
			range(i*os[0],(i+1)*os[0]),
			range(j*os[1],(j+1)*os[1])
			))=raster;
    }
	// The end of the row has a partial block.
    if (0 != ns[1] % os[1]) {
      auto col_lrange = range(jblocks*os[1],ns[1]);
      auto col_rrange = range(0,ns[1]%os[1]);
      BOOST_ASSERT(col_lrange.size() == col_rrange.size());

      auto end_block_view = project(morph,range(i*os[0],(i+1)*os[0]),col_lrange);
      auto end_raster_view = project(raster,range(0,os[0]),col_rrange);
      noalias(end_block_view) = end_raster_view;
    }
  }
  // Rows at the end may not have been copied yet.
  if (0 != ns[0] % os[0]) {
    const auto row_lrange = range(iblocks*os[0],ns[0]);
    const auto row_rrange = range(0,ns[0]%os[0]);
    BOOST_ASSERT(row_lrange.size() == row_rrange.size());

    for (int j=0; j<jblocks; j++) {
      const auto col_lrange = range(j*os[1],(j+1)*os[1]);
      const auto col_rrange = range(0,os[1]);
      BOOST_ASSERT(col_lrange.size() == col_rrange.size());

      auto chopped_block = project(morph,row_lrange,col_lrange);
      auto chopped_raster = project(raster,row_rrange,col_rrange);
      noalias(chopped_block) = chopped_raster;
    }

    if (0 != ns[1] % os[1]) {
      auto col_lrange = range(jblocks*os[1],ns[1]);
      auto col_rrange = range(0,ns[1]%os[1]);
      BOOST_ASSERT(col_lrange.size() == col_rrange.size());

      auto corner_block = project(morph,row_lrange,col_lrange);
      auto corner_raster = project(raster,row_rrange,col_rrange);
      noalias(corner_block) = corner_raster;
    }
  }
  return pmorph;
}



/*! Fills a portion of a raster of the given size with the given range of values.
 *  b is the bounds within the matrix.
 *  vals is the half-open range of values to assign, so (2,3] means assign 2
 *  to everything.
 */
size_t color_range(landscape_t& raster, boost::array<landscape_t::size_type,4> b,
                   boost::array<landscape_t::value_type,2> vals)
{
    typedef boost::array<landscape_t::size_type,4> dim_t;
    typedef boost::array<landscape_t::value_type,2> val_t;
 
    size_t color_cnt=0;
    BOOST_ASSERT(vals[1]>vals[0]);
    if (vals[1]-vals[0]==1 || (b[1]-b[0]==1 && b[3]-b[2]==1)) {
        auto tile=project(raster,range(b[0],b[1]),range(b[2],b[3]));
        noalias(tile) = scalar_matrix<landscape_t::value_type>(
                                    b[1]-b[0],b[3]-b[2],vals[0]
                                    );
        color_cnt += 1;
    } else {
        landscape_t::value_type midval = (vals[0]+vals[1])/2;
        BOOST_ASSERT(b[1]>b[0]);
        BOOST_ASSERT(b[3]>b[2]);
        if (b[1]-b[0] > b[3]-b[2]) {
            auto mid=(b[1]+b[0])/2;
            dim_t long_low = {b[0],mid,b[2],b[3]};
            val_t low = {vals[0],midval};
            color_cnt += color_range(raster,long_low,low);
            dim_t long_high = {mid,b[1],b[2],b[3]};
            val_t high = {midval,vals[1]};
            color_cnt += color_range(raster,long_high,high);
        } else {
            auto mid=(b[3]+b[2])/2;
            dim_t wide_left={b[0],b[1],b[2],mid};
            val_t left={vals[0],midval};
            color_cnt += color_range(raster,wide_left,left);
            dim_t wide_right={b[0],b[1],mid,b[3]};
            val_t right={midval,vals[1]};
            color_cnt += color_range(raster,wide_right,right);
        }
    }
    return color_cnt;
}



/*! Creates a raster of the given size with the given range of values.
 *  ns is the total size of the raster, and vals is the range of values.
 *  Both are half-open intervals.
 */
std::shared_ptr<landscape_t> multi_value(boost::array<landscape_t::size_type,2> ns,
                                           boost::array<landscape_t::value_type,2> vals)
{
    typedef landscape_t::value_type value_type;
    std::shared_ptr<landscape_t> pmorph(new landscape_t(ns[0],ns[1]));

    typedef boost::array<landscape_t::size_type,4> region_t;
    region_t whole = {0,ns[0],0,ns[1]};
    size_t color_cnt = color_range(*pmorph, whole, vals);
    BOOST_ASSERT( color_cnt == vals[1] - vals[0] );
    if ( color_cnt < vals[1]-vals[0] ) {
        throw runtime_error("The requested matrix was too small to hold all the values.");
    }
    return pmorph;
}





} // namespace
