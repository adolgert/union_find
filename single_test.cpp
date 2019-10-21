#include <map>
#include <boost/test/included/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/version.hpp>
#include "python2.7/patchlevel.h"
#include "array_store.hpp"
#include "single.hpp"
#include "grid2d.hpp"

using namespace std;
using namespace boost::unit_test;
using namespace raster_stats;


void test_sample()
{
	
	BOOST_CHECK(1,1);
}



bool init_function( )
{
  BOOST_TEST_MESSAGE("Using Boost version " << BOOST_VERSION);

  auto& master = framework::master_test_suite();
  master.add( BOOST_TEST_CASE( test_clusters_generic ) );
  return true;
}


struct rs_version {
	static std::string svn_version;
};

std::string rs_version::svn_version=RASTER_STATS_VERSION;


int main(int argc, char* argv[])
{
	cout << "Boost version " << BOOST_VERSION << endl << std::flush;
	cout << "Python version " << PY_MAJOR_VERSION << "." << PY_MINOR_VERSION <<
	    "." << PY_MICRO_VERSION << endl;
    cout << "TIFFLib version " << TIFFLIB_VERSION << endl;
	cout << "TBB version " << TBB_VERSION_MAJOR << "." << TBB_VERSION_MINOR << endl;
#ifdef __GNUC__
	cout << "Gnu C++ version " << __GNUC__ << "." << __GNUC_MINOR__ <<
		"." << __GNUC_PATCHLEVEL__ << endl;
    cout << "Raster stats version " << RASTER_STATS_VERSION << endl;
#endif // __GNUC__
	
    return ::boost::unit_test::unit_test_main(&init_function, argc, argv);
}
