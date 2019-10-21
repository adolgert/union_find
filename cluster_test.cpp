#include <memory>
#include <map>
#include <boost/test/included/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include "raster.hpp"
#include "io_geotiff.hpp"
#include "cluster.hpp"
#include "unique_values.hpp"
#include "array_store.hpp"
#include "grid2d.hpp"
#include "single.hpp"
#include "array_init.hpp"
#include "gather_clusters.hpp"
#include "morton.hpp"


using namespace std;
using namespace boost::unit_test;
using namespace raster_stats;


void test_transform()
{
	transform_ij tij(512);
    {
        transform_ij::key_type coord;
        coord[0]=0;
        coord[1]=0;
        BOOST_CHECK_EQUAL(tij(coord),0);
        coord[1]=100;
        BOOST_CHECK_EQUAL(tij(coord),100);
        coord[0]=2;
        BOOST_CHECK_EQUAL(tij(coord),512*2+100);
    }
	transform_map<transform_ij,unsigned char> data(tij,512*256);
    {
        transform_ij::key_type coord;
        coord[0]=100;
        coord[1]=100;
        data[coord]=42;
        unsigned char v=get(data,coord);
        BOOST_CHECK_EQUAL(v,42);
    }
}



template<class TR>
void write_read_transform(TR tij)
{
    // You could allocate 512*256, but the Morton order
    // uses a square allocation whether you like it or not.
    transform_map<TR,unsigned char> data(tij,512*512);

    for (size_t wi=0; wi<256; wi++) {
        for (size_t wj=0; wj<512; wj++) {
            typename TR::key_type coord;
            coord[0]=wi;
            coord[1]=wj;
            data[coord]=(wi*3571+wj*2663) % 256;
        }
    }

    for (size_t ri=0; ri<256; ri++) {
        for (size_t rj=0; rj<512; rj++) {
            typename TR::key_type coord;
            coord[0]=ri;
            coord[1]=rj;
            unsigned char v=get(data,coord);
            BOOST_CHECK_EQUAL(v,(ri*3571+rj*2663) % 256);
        }
    }
}



void test_transform_coverage()
{
	transform_ij tij(512);
	write_read_transform<transform_ij>(tij);
}



void test_blocked_coverage()
{
	transform_ij_blocked tij(512,32);
    write_read_transform<transform_ij_blocked>(tij);
}



void test_full_blocked()
{
	transform_ij_full_blocked tij(512,32);
    write_read_transform<transform_ij_full_blocked>(tij);
}



void test_morton()
{
    transform_morton_ij tij;
    write_read_transform<transform_morton_ij>(tij);
}



void test_bits()
{
    size_t a = alternating_bits<1,1>::value;
    BOOST_CHECK_EQUAL( a, 1);
    a = alternating_bits<8,1>::value;
    BOOST_CHECK_EQUAL( a, 0b0101010101010101);

    a = alternating_bits<1,2>::value;
    BOOST_CHECK_EQUAL( a, 2);
    a = alternating_bits<8,2>::value;
    BOOST_CHECK_EQUAL( a, 0b1010101010101010);

    a=morton_xy<0b10,0>::value;
    BOOST_CHECK_EQUAL( a, 0b0100);
    a=morton_xy<0b110,0>::value;
    BOOST_CHECK_EQUAL( a, 0b010100);

    a=morton_xy<0,0b10>::value;
    BOOST_CHECK_EQUAL( a, 0b1000);
    a=morton_xy<0,0b110>::value;
    BOOST_CHECK_EQUAL( a, 0b101000);

    a=morton_xy<1,0>::value;
    BOOST_CHECK_EQUAL(a,1);

    a=morton_xy<0,1>::value;
    BOOST_CHECK_EQUAL(a,0b10);

    a=morton_xy<-1,0>::value;
    size_t b=alternating_bits<32,1>::value;
    BOOST_CHECK_EQUAL(a,b);

    a=morton_xy<0,-1>::value;
    b=alternating_bits<32,0b10>::value;
    BOOST_CHECK_EQUAL(a,b);
}



void test_bit_shift()
{
    std::vector<boost::array<size_t,2>> a;
    a.push_back({{0b011,0b110 }});
    a.push_back({{0b110,0b011 }});
    //a.push_back({{0b0,0b0}}); // Interesting case.
    a.push_back({{0b111,0b111}});

    auto res=morton_calculations::combine_xy<size_t,size_t>({{ 0b011,0b110 }});
    BOOST_CHECK_EQUAL(res,0b101101);

    auto idx=a.begin();
    while (idx!=a.end()) {
        const boost::array<size_t,2>& val=*idx;
        size_t m = morton_calculations::combine_xy<size_t,size_t>(val);
        const boost::array<size_t,2>& res=morton_calculations::detangle<size_t,size_t>(m);
        BOOST_CHECK_EQUAL(val[0],res[0]);
        BOOST_CHECK_EQUAL(val[1],res[1]);

        size_t n;
        n = morton_calculations::add_interleaved(m,morton_xy<1,0>::value);
        const boost::array<size_t,2>& xp=morton_calculations::detangle<size_t,size_t>(n);
        BOOST_CHECK_EQUAL(val[0]+1,xp[0]);
        BOOST_CHECK_EQUAL(val[1],xp[1]);

        n = morton_calculations::add_interleaved(m,morton_xy<-1,0>::value);
        const boost::array<size_t,2>& xm=morton_calculations::detangle<size_t,size_t>(n);
        BOOST_CHECK_EQUAL(val[0]-1,xm[0]);
        BOOST_CHECK_EQUAL(val[1],xm[1]);

        n = morton_calculations::add_interleaved(m,morton_xy<0,1>::value);
        const boost::array<size_t,2>& yp=morton_calculations::detangle<size_t,size_t>(n);
        BOOST_CHECK_EQUAL(val[0],yp[0]);
        BOOST_CHECK_EQUAL(val[1]+1,yp[1]);

        n = morton_calculations::add_interleaved(m,morton_xy<0,-1>::value);
        const boost::array<size_t,2>& ym=morton_calculations::detangle<size_t,size_t>(n);
        BOOST_CHECK_EQUAL(val[0],ym[0]);
        BOOST_CHECK_EQUAL(val[1]-1,ym[1]);

        ++idx;
    }

}



void test_tiny_grid()
{
    typedef array_basis<size_t> basis_t;
    basis_t::bounds_type bounds;
    bounds[0][0]=0;
    bounds[0][1]=3;
    bounds[1][0]=0;
    bounds[1][1]=3;
    basis_t basis(bounds,1);
    BOOST_CHECK_EQUAL(basis.empty(),false);

    auto iter=make_vertex_iterator(basis);
    size_t neighbor_cnt=0;
    size_t vertex_cnt  =0;
    while (iter[0]!=iter[1]) {
        vertex_cnt++;
        auto adj=make_four_adjacent(basis,*iter[0]);
        while (adj[0]!=adj[1]) {
            neighbor_cnt++;
            adj[0]++;
        }
        iter[0]++;
    }
    BOOST_CHECK_EQUAL(vertex_cnt,9);
    BOOST_CHECK_EQUAL(neighbor_cnt,4*1+4*3+2*4);
}




void test_grid_basis()
{
    typedef array_basis<size_t> basis_t;
    basis_t::bounds_type bounds;
    bounds[0][0]=0;
    bounds[0][1]=256;
    bounds[1][0]=0;
    bounds[1][1]=512;
    basis_t basis(bounds,32);
    BOOST_CHECK_EQUAL(basis.empty(),false);

    auto iter=make_vertex_iterator(basis);
    size_t neighbor_cnt=0;
    size_t vertex_cnt  =0;
    while (iter[0]!=iter[1]) {
        vertex_cnt++;
        auto adj=make_four_adjacent(basis,*iter[0]);
        while (adj[0]!=adj[1]) {
            neighbor_cnt++;
            adj[0]++;
        }
        iter[0]++;
    }
    BOOST_CHECK_EQUAL(vertex_cnt,512*256);
    BOOST_CHECK_EQUAL(neighbor_cnt,4*(512-2)*(256-2)+3*(510+254)*2+2*4);
}



// std::map has third and fourth default template arguments that this hides.
template<typename FROM, typename TO>
struct bound_map : std::map<FROM,TO> {};



void test_generic_single()
{
    typedef unsigned char value_type;
    boost::array<size_t,2> extent;
    extent[0]=100;
    extent[1]=100;

    typedef array_basis<size_t> basis_t;
    basis_t::bounds_type bounds;
    bounds[0][0]=0;
    bounds[0][1]=extent[0];
    bounds[1][0]=0;
    bounds[1][1]=extent[1];
    basis_t basis(bounds,32);
	transform_ij tij(extent[0]);
    typedef transform_map<transform_ij,value_type> map_t;
	map_t data(tij,extent[0]*extent[1]);

    boost::array<value_type,2> limits;
    limits[0]=0;
    limits[1]=100;
    checkerboard_array(data,extent,limits);

    typedef AreEqual<map_t::key_type,map_t> comparison_t;
    comparison_t comparison(data);

    typedef construct_disjoint_set<basis_t::vertex_type,
                                   bound_map,
                                   bound_map> disj_t;
    // Just to prove that disj_t is well-formed, let's
    // instantiate one, but really it's a policy class.
    disj_t example;

    union_find_st<disj_t> ufind;
    ufind(basis,comparison,
          make_vertex_iterator<basis_t>,
          make_four_adjacent<basis_t>);
    
    auto clusters = gather_clusters(ufind.rank_pmap_,ufind.dset_,
                                         extent);
    BOOST_CHECK_EQUAL(clusters->size(),limits[1]-limits[0]);
}




void test_generic_blocked()
{
    typedef unsigned char value_type;
    boost::array<size_t,2> extent;
    extent[0]=100;
    extent[1]=100;

    typedef array_basis<size_t> basis_t;
    basis_t::bounds_type bounds;
    bounds[0][0]=0;
    bounds[0][1]=extent[0];
    bounds[1][0]=0;
    bounds[1][1]=extent[1];
    basis_t basis(bounds,32);
    transform_ij_blocked tij(extent[0],10);
    typedef transform_map<transform_ij_blocked,value_type> map_t;
    map_t data(tij,extent[0]*extent[1]);

    boost::array<value_type,2> limits;
    limits[0]=0;
    limits[1]=100;
    checkerboard_array(data,extent,limits);

    typedef AreEqual<map_t::key_type,map_t> comparison_t;
    comparison_t comparison(data);

    typedef construct_disjoint_set<basis_t::vertex_type,
                                   bound_map,
                                   bound_map> disj_t;
    // Just to prove that disj_t is well-formed, let's
    // instantiate one, but really it's a policy class.
    disj_t example;

    union_find_st<disj_t> ufind;
    ufind(basis,comparison,
          make_vertex_iterator<basis_t>,
          make_four_adjacent<basis_t>);

    auto clusters = gather_clusters(ufind.rank_pmap_,ufind.dset_,
                                         extent);
    BOOST_CHECK_EQUAL(clusters->size(),limits[1]-limits[0]);
}






void test_generic_full()
{
    typedef unsigned char value_type;
    boost::array<size_t,2> extent;
    extent[0]=100;
    extent[1]=100;

    typedef array_basis<size_t> basis_t;
    basis_t::bounds_type bounds;
    bounds[0][0]=0;
    bounds[0][1]=extent[0];
    bounds[1][0]=0;
    bounds[1][1]=extent[1];
    basis_t basis(bounds,32);
    transform_ij_full_blocked tij(extent[0],10);
    typedef transform_map<transform_ij_full_blocked,value_type> map_t;
    map_t data(tij,extent[0]*extent[1]);

    boost::array<value_type,2> limits;
    limits[0]=0;
    limits[1]=100;
    checkerboard_array(data,extent,limits);

    typedef AreEqual<map_t::key_type,map_t> comparison_t;
    comparison_t comparison(data);

    typedef construct_disjoint_set<basis_t::vertex_type,
                                   bound_map,
                                   bound_map> disj_t;
    // Just to prove that disj_t is well-formed, let's
    // instantiate one, but really it's a policy class.
    disj_t example;

    union_find_st<disj_t> ufind;
    ufind(basis,comparison,
          make_vertex_iterator<basis_t>,
          make_four_adjacent<basis_t>);

    auto clusters = gather_clusters(ufind.rank_pmap_,ufind.dset_,
                                         extent);
    BOOST_CHECK_EQUAL(clusters->size(),limits[1]-limits[0]);
}




void test_unique_values()
{
  boost::shared_ptr<landscape_t> raster = read_tiff("34418039.tif");
  std::set<landscape_t::value_type> uniques;
  unique_values(*raster,uniques);
  BOOST_CHECK_EQUAL(uniques.size(),15);
}



void known_single_blank()
{
    boost::shared_ptr<landscape_t> raster = multi_value({{100,100}},{{0,1}});
    cluster_t clusters = find_clusters(*raster);
    BOOST_CHECK_EQUAL(clusters.size(),1);
}


void known_many_blank()
{
    boost::shared_ptr<landscape_t> raster = multi_value({{100,100}},{{0,25}});
    cluster_t clusters = find_clusters(*raster);
    BOOST_CHECK_EQUAL(clusters.size(),25);
}




void known_single_twopass()
{
    boost::shared_ptr<landscape_t> raster = multi_value({{100,100}},{{0,1}});
    cluster_t clusters = find_clusters_twopass(*raster);
    BOOST_CHECK_EQUAL(clusters.size(),1);
}


void known_many_twopass()
{
    boost::shared_ptr<landscape_t> raster = multi_value({{100,100}},{{0,25}});
    cluster_t clusters = find_clusters_twopass(*raster);
    BOOST_CHECK_EQUAL(clusters.size(),25);
}




void known_single_pair()
{
    boost::shared_ptr<landscape_t> raster = multi_value({{100,100}},{{0,1}});
    cluster_loc_t clusters = find_clusters_pair(*raster);
    BOOST_CHECK_EQUAL(clusters.size(),1);
}


void known_many_pair()
{
    boost::shared_ptr<landscape_t> raster = multi_value({{100,100}},{{0,25}});
    cluster_loc_t clusters = find_clusters_pair(*raster);
    BOOST_CHECK_EQUAL(clusters.size(),25);
}


void known_single_remap()
{
    boost::shared_ptr<landscape_t> raster = multi_value({{100,100}},{{0,1}});
    cluster_t clusters = find_clusters_remap(*raster);
    BOOST_CHECK_EQUAL(clusters.size(),1);
}


void known_many_remap()
{
    boost::shared_ptr<landscape_t> raster = multi_value({{100,100}},{{0,25}});
    cluster_t clusters = find_clusters_remap(*raster);
    BOOST_CHECK_EQUAL(clusters.size(),25);
}



bool init_function( )
{
  BOOST_TEST_MESSAGE("Loading from 34418039.tif");

  framework::master_test_suite().add( BOOST_TEST_CASE( test_transform ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( test_transform_coverage));
  framework::master_test_suite().add( BOOST_TEST_CASE( test_blocked_coverage ));
  framework::master_test_suite().add( BOOST_TEST_CASE( test_full_blocked ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( test_morton ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( test_bits ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( test_bit_shift ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( test_tiny_grid ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( test_grid_basis ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( test_generic_single ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( test_generic_blocked ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( test_generic_full ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( test_unique_values ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( known_single_blank ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( known_many_blank ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( known_single_twopass ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( known_many_twopass ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( known_single_pair ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( known_many_pair ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( known_single_remap ) );
  framework::master_test_suite().add( BOOST_TEST_CASE( known_many_remap ) );
  return true;
}


int main(int argc, char* argv[])
{
  return ::boost::unit_test::unit_test_main(&init_function, argc, argv);
}

