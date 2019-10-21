#ifndef _GATHER_CLUSTERS_HPP_
#define _GATHER_CLUSTERS_HPP_ 1

#include <map>
#include <list>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include "raster.hpp"

namespace raster_stats {

template<typename parent_map, typename disjoint_set>
  boost::shared_ptr<cluster_t> gather_clusters(parent_map& parent,
                                        disjoint_set& dset,
                                        size_t icnt, size_t jcnt)
{
  // At this point, each element points to a parent.
  // The return value is a list of lists of elements.
  // The temporary map will associate a parent with the list of its children.
  boost::shared_ptr<cluster_t> clusters(new cluster_t); // a list of lists

  typedef std::map<size_t,cluster_t::iterator> ptl_t;
  ptl_t parent_to_list; // The map from parent to list of children.

  for (size_t pull=0; pull<icnt*jcnt; pull++) {
    size_t parent = dset.find_set(pull);
    ptl_t::iterator plist = parent_to_list.find(parent);
    if (plist==parent_to_list.end()) {
      cluster_t::iterator nlist = clusters->insert(clusters->end(),
                                                   std::list<size_t>());
      parent_to_list[parent]=nlist;
      nlist->push_back(pull);
    } else {
      plist->second->push_back(pull);
    }
  }
  return clusters;
}



template<typename parent_map, typename disjoint_set>
boost::shared_ptr<std::list<std::list<typename parent_map::key_type>>>
gather_clusters(parent_map& parent,
                disjoint_set& dset,
                boost::array<size_t,2> dim)
{
    // At this point, each element points to a parent.
    // The return value is a list of lists of elements.
    // The temporary map will associate a parent with the list of its children.
    typedef typename parent_map::key_type loc_t;
    typedef std::list<std::list<loc_t>> clus_t;
    auto clusters = boost::make_shared<clus_t>(); // a list of lists

    typedef std::map<loc_t,typename clus_t::iterator> ptl_t;
    ptl_t parent_to_list; // The map from parent to list of children.
    
    for (size_t i=0; i<dim[0]; i++) {
        for (size_t j=0; j<dim[1]; j++) {
            loc_t loc;
            loc[0]=i;
            loc[1]=j;
            auto const parent = dset.find_set(loc);
            typename ptl_t::iterator plist = parent_to_list.find(parent);
            if (plist==parent_to_list.end()) {
                auto nlist = clusters->insert(clusters->end(),
                                 std::list<loc_t>());
                parent_to_list[parent]=nlist;
                nlist->push_back(loc);
            } else {
                plist->second->push_back(loc);
            }
        }
    }
    return clusters;
}


}

#endif
