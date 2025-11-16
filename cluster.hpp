#ifndef _CLUSTER_H_
#define _CLUSTER_H_ 1

#include <set>
#include <list>
#include <set>
#include <memory>
#include "raster.hpp"
#include "gather_clusters.hpp"

namespace raster_stats {

typedef unsigned char arr_type;
std::set<arr_type> unique_values_direct(const landscape_t& raster);
cluster_t find_clusters(const landscape_t& raster);
cluster_t find_clusters_twopass(const landscape_t& raster);
std::shared_ptr<cluster_t> find_clusters_pointer(const landscape_t& raster);
cluster_loc_t find_clusters_pair(const landscape_t& raster);
cluster_t find_clusters_remap(const landscape_t& raster);





}

#endif
