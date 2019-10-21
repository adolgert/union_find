#include <map>
#include <boost/test/included/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/version.hpp>
#include "tbb/task_scheduler_init.h"
#include "tbb/tbb_stddef.h"
#include "python2.7/patchlevel.h"
#include "raster.hpp"
#include "raster_version.hpp"
#include "cluster_generic.hpp"
#include "io_geotiff.hpp"
#include "tiffvers.h"
#include "cluster_tbb.hpp"


using namespace std;
using namespace boost::unit_test;
using namespace raster_stats;


void test_clusters_generic()
{
    std::cout << "test_clusters_generic: enter" << std::endl;
    const int thread_cnt = 1;
    tbb::task_scheduler_init init(thread_cnt);
    boost::shared_ptr<landscape_t> raster = multi_value({{100,100}},{{0,1}});
	cluster_raster<landscape_t> cr;
	cr(*raster);
	size_t cluster_cnt = count(cr);
	BOOST_CHECK_EQUAL(cluster_cnt,1);
}


void test_clusters_basis()
{
	// This array is (0,10) in the i and (20,40) in the j.
	boost::array<size_t,4> bounds = {{0,10,20,40}};
	array_basis basis(bounds,32);
	BOOST_CHECK_EQUAL(*basis.begin(),0);
	BOOST_CHECK_EQUAL(*basis.end(),  200);
	array_basis half(basis,tbb::split());
	BOOST_CHECK_EQUAL(*basis.begin(),0);
	BOOST_CHECK_EQUAL(*basis.end(),  200);
	BOOST_CHECK_EQUAL(*half.begin(), 10);
	BOOST_CHECK_EQUAL(*half.end(),   210);

    auto iter = basis.begin();
    auto stop = basis.end();
    for (size_t i=0; i<10; ++i) {
        for (size_t j=0; j<10; ++j) {
            BOOST_CHECK_EQUAL(*iter,i*20+j);
            ++iter;
        }
    }
    BOOST_CHECK_EQUAL(*iter,*stop);
    auto halfiter = half.begin();
    auto halfstop=half.end();
    size_t halfcount = 0;
    while (halfiter!=halfstop) {
        ++halfiter;
        ++halfcount;
    }
    BOOST_CHECK_EQUAL(halfcount,100);
}



void test_clusters_adjacent()
{
	boost::array<size_t,4> bounds = {{0,10,20,40}};
	array_basis basis(bounds,32);
	array_basis half(basis,tbb::split());
    adjacent_iterator::loc_type loc={{5,30}};
    adjacent_iterator begin(bounds,bounds,loc,-1);
    adjacent_iterator end(bounds,bounds,loc,4);
    size_t count=0;
    while (begin!=end) {
        count++;
        begin++;
    }
    BOOST_CHECK_EQUAL(count,4);

    std::map<size_t,size_t> neighbors;
	auto basis_iter = basis.begin();
	auto basis_end = basis.end();
	while (basis_iter != basis_end) {
		boost::array<adjacent_iterator,2> adj = basis_iter.adjacent();
        size_t adj_cnt=0;
		while (adj[0]!=adj[1]) {
			++adj[0];
            ++adj_cnt;
	    }
        if (neighbors.find(adj_cnt)==neighbors.end()) {
            neighbors[adj_cnt]=0;
        }
        neighbors[adj_cnt]++;
		++basis_iter;
	}
	// The basis is 10 x 10. That means
    // 8x8 - 4 neighbors
    // 4*8 - 3 neighbors
    // 4 - 2 neighbors
    // 64 + 32 + 4 = 100.
    BOOST_CHECK_EQUAL(neighbors[4],64);
    BOOST_CHECK_EQUAL(neighbors[3],32);
    BOOST_CHECK_EQUAL(neighbors[2],4);
}



void test_edge_adjacent()
{
	boost::array<size_t,4> bounds = {{0,10,20,40}};
	array_basis basis(bounds,32);
	array_basis half(basis,tbb::split());
    // basis is 0,10,20,30
    // half is  0,10,30,40
    auto edge_begin=edge_iterator(basis,half,false);
    auto edge_end=edge_iterator(basis,half,true);
    size_t cnt=0;
    for ( ; edge_begin!=edge_end; edge_begin++) {
        boost::array<size_t,2> const& edge = *edge_begin;
        cnt++;
    }
    BOOST_CHECK_EQUAL(cnt,10);
}



bool pt_in_rect(const boost::array<size_t,2>& pt,
                const boost::array<size_t,4>& rect)
{
    return (pt[0]>=rect[0] && pt[0]<rect[1] && pt[1]>=rect[2] && pt[1]<rect[3]);
}



size_t compare_adjacent(const boost::array<size_t,4>& a_bounds,
         const boost::array<size_t,4>& b_bounds)
{
	array_basis a(a_bounds,32);
	array_basis b(b_bounds,32);
    auto edge_begin=edge_iterator(a,b,false);
    auto edge_end=edge_iterator(a,b,true);
    size_t cnt=0;
    for ( ; edge_begin!=edge_end; edge_begin++) {
        boost::array<size_t,2> const& edge = *edge_begin;
        boost::array<boost::array<size_t,2>,2> pt=edge_begin.coords();
        bool la=pt_in_rect(pt[0],a_bounds);
        bool lb=pt_in_rect(pt[0],b_bounds);
        bool ra=pt_in_rect(pt[1],a_bounds);
        bool rb=pt_in_rect(pt[1],b_bounds);
        BOOST_CHECK_EQUAL((la&&rb) || (lb&&ra),true);
        BOOST_CHECK_EQUAL(!(la&&ra) && !(lb&&rb),true);
        cnt++;
    }
    return cnt;
}



void test_adjacent_types()
{
    size_t cnt=0;
    typedef boost::array<size_t,4> bounds_t;
	bounds_t a_bounds = {{ 0,10,20,30}};
	bounds_t b_bounds = {{ 0,10,30,40}};
    bounds_t c_bounds = {{10,20,20,30}};
    bounds_t d_bounds = {{10,20,31,41}};
	bounds_t e_bounds = {{ 0,10,31,41}};
    cnt=compare_adjacent(a_bounds,b_bounds);
    BOOST_CHECK_EQUAL(cnt,10);
    cnt=compare_adjacent(b_bounds,a_bounds);
    BOOST_CHECK_EQUAL(cnt,10);
    cnt=compare_adjacent(a_bounds,c_bounds);
    BOOST_CHECK_EQUAL(cnt,10);
    cnt=compare_adjacent(c_bounds,a_bounds);
    BOOST_CHECK_EQUAL(cnt,10);
    cnt=compare_adjacent(a_bounds,d_bounds);
    BOOST_CHECK_EQUAL(cnt,0);    
    cnt=compare_adjacent(a_bounds,e_bounds);
    BOOST_CHECK_EQUAL(cnt,0);    
}




void test_clusters_tiny_tbb0()
{
    const int thread_cnt = 1;
    tbb::task_scheduler_init init(thread_cnt);
    boost::shared_ptr<landscape_t> raster = resize_replicate(read_tiff("34418039.tif"),{{10,10}});
    boost::shared_ptr<cluster_t> clusters = clusters_tbb0(*raster);
}



void test_clusters_tbb0()
{
    const int thread_cnt = 1;
    tbb::task_scheduler_init init(thread_cnt);
    boost::shared_ptr<landscape_t> raster = resize_replicate(read_tiff("34418039.tif"),{{200,200}});
    boost::shared_ptr<cluster_t> clusters = clusters_tbb0(*raster);
}


void known_single_thread_tbb0()
{
    const int thread_cnt = 1;
    tbb::task_scheduler_init init(thread_cnt);
    boost::shared_ptr<landscape_t> raster = multi_value({{100,100}},{{0,1}});
    boost::shared_ptr<cluster_t> clusters = clusters_tbb0(*raster);
    BOOST_CHECK_EQUAL(clusters->size(),1);
}



void known_single_tbb0()
{
    const int thread_cnt = 6;
    tbb::task_scheduler_init init(thread_cnt);
    boost::shared_ptr<landscape_t> raster = multi_value({{100,100}},{{0,1}});
    boost::shared_ptr<cluster_t> clusters = clusters_tbb0(*raster);
    BOOST_CHECK_EQUAL(clusters->size(),1);
}



void known_many_tbb0()
{
    const int thread_cnt = 6;
    tbb::task_scheduler_init init(thread_cnt);
    boost::shared_ptr<landscape_t> raster = multi_value({{100,100}},{{0,25}});
    boost::shared_ptr<cluster_t> clusters = clusters_tbb0(*raster);
    BOOST_CHECK_EQUAL(clusters->size(),25);
}



bool init_function( )
{
  BOOST_TEST_MESSAGE("Using Boost version " << BOOST_VERSION);

  auto& master = framework::master_test_suite();
  master.add( BOOST_TEST_CASE( test_clusters_generic ) );
  master.add( BOOST_TEST_CASE( test_clusters_basis ) );
  master.add( BOOST_TEST_CASE( test_clusters_adjacent ) );
  master.add( BOOST_TEST_CASE( test_edge_adjacent ) );
  master.add( BOOST_TEST_CASE( test_adjacent_types ) );
  master.add( BOOST_TEST_CASE( known_single_thread_tbb0 ) );
  master.add( BOOST_TEST_CASE( test_clusters_tiny_tbb0 ) );
  master.add( BOOST_TEST_CASE( test_clusters_tbb0 ) );
  master.add( BOOST_TEST_CASE( known_single_tbb0 ) );
  master.add( BOOST_TEST_CASE( known_many_tbb0 ) );
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
