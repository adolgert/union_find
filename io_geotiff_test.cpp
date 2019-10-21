#define BOOST_TEST_MODULE io_geotiff
#include <boost/test/included/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/array.hpp>
#include "raster.hpp"
#include "io_geotiff.hpp"

using namespace std;
using namespace boost::unit_test;
using namespace raster_stats;

#define SMALL_TIFF "feep.tif"

BOOST_AUTO_TEST_SUITE( geotiff )

BOOST_AUTO_TEST_CASE( read_dims )
{
  auto dims = tiff_dimensions(SMALL_TIFF);
  BOOST_CHECK_EQUAL(dims[0],7);
  BOOST_CHECK_EQUAL(dims[1],24);
}

BOOST_AUTO_TEST_CASE( read_vals )
{
  tiff_data_format(SMALL_TIFF);
  auto landscape = read_tiff(SMALL_TIFF);
  BOOST_CHECK_EQUAL((*landscape)(0,0),0);
  BOOST_CHECK_EQUAL((*landscape)(1,1),3);
  BOOST_CHECK_EQUAL((*landscape)(5,1),3);
  BOOST_CHECK_EQUAL((*landscape)(1,19),17); // It's flipped.
  BOOST_CHECK_EQUAL((*landscape)(5,19),253); // It's flipped.
}

BOOST_AUTO_TEST_CASE( resize_replicate_test )
{
  auto landscape = read_tiff(SMALL_TIFF);
  typedef boost::array<size_t,2> dim_t;
  dim_t size = {{2,2}};
  auto smaller = resize_replicate(landscape, size);
  BOOST_CHECK_EQUAL((*smaller)(0,0),0);
  BOOST_CHECK_EQUAL((*smaller)(1,1),3);

  dim_t larger = {{7,48}};
  auto larger_right = resize_replicate(landscape, larger);
  BOOST_CHECK_EQUAL((*larger_right)(1,1),3);
  BOOST_CHECK_EQUAL((*larger_right)(1,19),17); // It's flipped.

  dim_t lright = {{7,50}};
  auto llarger_right = resize_replicate(landscape, lright);
  BOOST_CHECK_EQUAL((*llarger_right)(1,49),3);

  dim_t ldown = {{14,24}};
  auto larger_down = resize_replicate(landscape,ldown);
  BOOST_CHECK_EQUAL((*larger_down)(8,1),3);
}

BOOST_AUTO_TEST_CASE( generate_single_raster )
{
  typedef boost::array<size_t,2> dim_t;
  dim_t xy = {{100,100}};
  boost::array<unsigned char,2> vals = {{3,4}};
  auto landscape = multi_value(xy,vals);
  BOOST_CHECK_EQUAL((*landscape)(0,0),3);
  BOOST_CHECK_EQUAL((*landscape)(99,99),3);
  BOOST_CHECK_EQUAL((*landscape)(25,25),3);
  for (size_t i=0; i<100; i++) {
    for (size_t j=0; j<100; j++) {
      BOOST_CHECK_EQUAL(landscape->operator()(i,j),3);
    }
  }
}



BOOST_AUTO_TEST_CASE( generate_raster )
{
  typedef boost::array<size_t,2> dim_t;
  dim_t xy = {{100,100}};
  boost::array<unsigned char,2> vals = {{3,36}};
  auto landscape = multi_value(xy,vals);
    std::set<landscape_t::value_type> types;
    for (size_t i=0; i<100; i++) {
        for (size_t j=0; j<100; j++) {
            types.insert((*landscape)(i,j));
        }
    }
    BOOST_CHECK_EQUAL(types.size(),33);
}


BOOST_AUTO_TEST_SUITE_END()
