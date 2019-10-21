#include "cluster.hpp"
#include "unique_values.hpp"
#include "timing.hpp"
#include "raster_times.hpp"

using namespace std;

namespace raster_stats {

vector<timing_harness> raster_times() {
	vector<timing_harness> tests;
	//tests.push_back([](const landscape_t& raster) { unique_values(raster); });
	timing_harness::subject_t st0 = [](const landscape_t& raster) { 
        std::set<landscape_t::value_type> uniques;
        unique_values(raster,uniques); };
	tests.push_back(timing_harness(st0,"unique_values"));
	timing_harness::subject_t st1 = [](const landscape_t& raster) {
        unique_values_direct(raster); };
	tests.push_back(timing_harness(st1,"unique_values_direct"));
	timing_harness::subject_t st2 = [](const landscape_t& raster) { find_clusters(raster); };
	tests.push_back(timing_harness(st2,"find_clusters"));
	return tests;
}


}
