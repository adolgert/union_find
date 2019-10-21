/*! main.cpp
 *  This is intended for quick testing without Python.
 */

#include <iostream>
#include <functional>
#include <boost/program_options.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include "io_geotiff.hpp"
#include "io_ppm.hpp"
#include "unique_values.hpp"
#include "cluster.hpp"
#include "timing.hpp"
#include "timing_harness.hpp"
#include "single_timing.hpp"
#include "io_hdf.hpp"

using namespace std;
using namespace raster_stats;
namespace po = boost::program_options;




void do_stuff()
{
    double start=0.9;
    for (size_t i=0; i<9999; i++) {
        start=start*0.9999;
    }
}




/*! Randomly reorder a random-access container with Knuth-Fisher-Yates.
 */
template<class RandomAccessIter,class Generator>
void random_reorder(RandomAccessIter begin, RandomAccessIter end,
                    Generator& rng)
{
    for (RandomAccessIter cursor = end-1; cursor>begin; --cursor) {
        auto gen =
            boost::random::uniform_int_distribution<size_t>(0,
                             distance(begin,cursor));
        swap(*(begin+gen(rng)),*cursor);
    }
}



int main(int argc, char* argv[])
{
    po::options_description desc("Allowed options");
    size_t side_length;
    size_t iterations;
    size_t count;
    size_t depth;
    size_t block;
    desc.add_options()
        ("help","show help message")
        ("size,s", po::value<size_t>(&side_length)->default_value(100),
         "length of a side of the raster")
        ("depth,d", po::value<size_t>(&depth)->default_value(100),
         "number of land use types")
        ("block,b", po::value<size_t>(&block)->default_value(32),
         "size of blocks")
        ("iter,i",po::value<size_t>(&iterations)->default_value(1),
         "number of times to run test during a single timing run")
        ("count,c",po::value<size_t>(&count)->default_value(1),
         "number of times to run sets of iterations of all tests")
        ("tiff", po::value<std::string>(),"filename of a TIFF to read")
        ;

    po::variables_map vm;
    auto parsed_options=po::parse_command_line(argc,argv,desc);
    po::store(parsed_options, vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << endl;
        return 1;
    }

    boost::random::mt19937 rng;
    rng.seed(static_cast<unsigned int>(std::time(0)));


    std::string tiff_filename;
    if (vm.count("tiff")) {
        tiff_filename = vm["tiff"].as<std::string>();
	}
	
	vector<boost::shared_ptr<timing_harness>> tests;

    // Create the first test using simple ij
    typedef array_basis<size_t> basis_t;
    typedef transform_map<transform_ij,unsigned char> map_t;

    auto bd=make_data<basis_t,map_t>(side_length,side_length,block,depth);
    auto basis=bd.get<0>();
    auto data=bd.get<1>();

    auto run=single_run<basis_t,map_t>(basis,data);
    auto timing0=make_timing(run,"single");
    tests.push_back(timing0);

    // Create the second test with a blocked transform.
    typedef transform_map<transform_ij_blocked,unsigned char> mapb_t;

    auto bdb=make_data<basis_t,mapb_t>(side_length,side_length,block,depth);
    auto basisb=bdb.get<0>();
    auto datab=bdb.get<1>();

    auto run1=single_run<basis_t,mapb_t>(basisb,datab);
    auto timing1=make_timing(run1,"blocked");
    tests.push_back(timing1);

    // The third test uses a full blocking.
    typedef transform_map<transform_ij_full_blocked,unsigned char> mapc_t;

    auto bdc=make_data<basis_t,mapc_t>(side_length,side_length,block,depth);
    auto basisc=bdc.get<0>();
    auto datac=bdc.get<1>();

    auto run2=single_run<basis_t,mapc_t>(basisc,datac);
    auto timing2=make_timing(run2,"full_blocked");
    tests.push_back(timing2);


    // Run all tests, randomizing the order for each set of runs.
    // Store results in a list for later.
    std::vector<size_t> order(tests.size());
    std::vector<std::vector<boost::array<size_t,2>>> results(tests.size());
    for (size_t run_cnt=0; run_cnt<count; run_cnt++) {
        for (size_t j=0; j<order.size(); j++) { order[j]=j; }
        random_reorder(begin(order),end(order),rng);
        for (size_t test_idx=0; test_idx<tests.size(); test_idx++) {
            size_t which_test=order[test_idx];
            auto test=tests[which_test];
            auto time = test->time(iterations);
            results[which_test].push_back({{iterations,time}});
        }
    }

    timing_file out_file("out.h5", parsed_options);

    // Save results
    for (size_t save_idx=0; save_idx<tests.size(); save_idx++) {
        const std::vector<boost::array<size_t,2>>& times = results[save_idx];
        std::string test_name(tests[save_idx]->name());
        out_file.store_test(times,test_name);
    }


    // Print results
    for (size_t test_result=0; test_result<tests.size(); test_result++) {
        const auto& test_list=results[test_result];
        cout << tests[test_result]->name() << "\t";
        auto presult = test_list.begin();
        for (; presult!=test_list.end();presult++) {
            cout << " " << (*presult)[0] << " " << (*presult)[1];
        }
        cout << endl;
    }
    return 0;
}
