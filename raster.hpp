#ifndef _RASTER_H_
#define _RASTER_H_ 1

#include <map>
#include <utility>
#include <list>
#include <boost/numeric/ublas/matrix.hpp>

namespace raster_stats {

//! The raster image holding the landscape characteristics.
typedef boost::numeric::ublas::matrix<unsigned char> landscape_t;
//! A map from an individual quadrant to a list of neighboring, similar quadrants.
typedef std::list<std::list<size_t> > cluster_t;
//! The (i,j) coordinates of a quadrant of the landscape.
typedef std::pair<size_t,size_t> loc_t;
typedef std::map<loc_t,std::list<loc_t> > cluster_loc_t;

}

#endif
