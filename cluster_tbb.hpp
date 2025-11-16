#ifndef _CLUSTER_TBB_H_
#define _CLUSTER_TBB_H_ 1

#include <memory>
#include "raster.hpp"

namespace raster_stats {

  std::shared_ptr<cluster_t> clusters_tbb0(const landscape_t& raster);

}



#endif // _CLUSTER_TBB_H_
