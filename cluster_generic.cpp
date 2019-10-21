#include <boost/numeric/interval.hpp>
#include "cluster_generic.hpp"


namespace raster_stats {

    edge_iterator::edge_iterator(const array_basis& a, const array_basis& b,
                                 bool end=false) :
        whole_(a.whole_)
    {
        vertical_=false;
        start_   ={{0,0}};
        count_   =0;
        idx_     =0;
    struct rel {
        static bool above(const array_basis&a, const array_basis&b) {
            return a.bounds_[3]==b.bounds_[2]; }
        static bool rightof(const array_basis&a, const array_basis&b) {
            return a.bounds_[1]==b.bounds_[0]; }
    };
    if (rel::above(a,b)) {
        // intervals are closed, while bounds are half-open, so -1.
        boost::numeric::interval<size_t> ai(a.bounds_[0],a.bounds_[1]-1);
        boost::numeric::interval<size_t> bi(b.bounds_[0],b.bounds_[1]-1);
        if (boost::numeric::overlap(ai,bi)) {
            auto intersection = intersect(ai,bi);
            vertical_=true;
            start_={{boost::numeric::lower(intersection),a.bounds_[3]-1}};
            count_=boost::numeric::width(intersection)+1;
        }
    }
    else if (rel::above(b,a)) {
        boost::numeric::interval<size_t> ai(a.bounds_[0],a.bounds_[1]-1);
        boost::numeric::interval<size_t> bi(b.bounds_[0],b.bounds_[1]-1);
        if (boost::numeric::overlap(ai,bi)) {
            auto intersection = intersect(ai,bi);
            vertical_=true;
            start_={{boost::numeric::lower(intersection),b.bounds_[3]-1}};
            count_=boost::numeric::width(intersection)+1;
        }
    }
    else if (rel::rightof(a,b)) {
        boost::numeric::interval<size_t> ai(a.bounds_[2],a.bounds_[3]-1);
        boost::numeric::interval<size_t> bi(b.bounds_[2],b.bounds_[3]-1);
        if (boost::numeric::overlap(ai,bi)) {
            auto intersection = intersect(ai,bi);
            vertical_=false;
            start_={{a.bounds_[1]-1,boost::numeric::lower(intersection)}};
            count_=boost::numeric::width(intersection)+1;
        }
    }
    else if (rel::rightof(b,a)) {
        boost::numeric::interval<size_t> ai(a.bounds_[2],a.bounds_[3]-1);
        boost::numeric::interval<size_t> bi(b.bounds_[2],b.bounds_[3]-1);
        if (boost::numeric::overlap(ai,bi)) {
            auto intersection = intersect(ai,bi);
            vertical_=false;
            start_={{b.bounds_[1]-1,boost::numeric::lower(intersection)}};
            count_=boost::numeric::width(intersection)+1;
        }
    }
    if (end) {
        idx_=count_;
    }
}



}
