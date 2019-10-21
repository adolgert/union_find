#ifndef _ARRAY_INIT_H_
#define _ARRAY_INIT_H_ 1


#include <exception>
#include <boost/array.hpp>



/*! Fills a portion of a raster of the given size with the given range of values.
 *  b is the bounds within the matrix.
 *  vals is the half-open range of values to assign, so (2,3] means assign 2
 *  to everything.
 */
template<class ARRAY,class VT>
size_t checkerboard_range(ARRAY& raster,
                          boost::array<size_t,4> b,
                          boost::array<VT,2> vals)
{
    typedef boost::array<size_t,4> dim_t;
    typedef boost::array<VT,2> val_t;

    size_t color_cnt=0;
    BOOST_ASSERT(vals[1]>vals[0]);
    if (vals[1]-vals[0]==1 || (b[1]-b[0]==1 && b[3]-b[2]==1)) {
        for (size_t i=b[0]; i<b[1]; i++) {
            for (size_t j=b[2]; j<b[3]; j++) {
                boost::array<size_t,2> v;
                v[0]=i;
                v[1]=j;
                raster[v]=vals[0];
            }
        }
        color_cnt += 1;
    } else {
        VT midval = (vals[0]+vals[1])/2;
        BOOST_ASSERT(b[1]>b[0]);
        BOOST_ASSERT(b[3]>b[2]);
        if (b[1]-b[0] > b[3]-b[2]) {
            auto mid=(b[1]+b[0])/2;
            dim_t long_low = {b[0],mid,b[2],b[3]};
            val_t low = {vals[0],midval};
            color_cnt += checkerboard_range(raster,long_low,low);
            dim_t long_high = {mid,b[1],b[2],b[3]};
            val_t high = {midval,vals[1]};
            color_cnt += checkerboard_range(raster,long_high,high);
        } else {
            auto mid=(b[3]+b[2])/2;
            dim_t wide_left={b[0],b[1],b[2],mid};
            val_t left={vals[0],midval};
            color_cnt += checkerboard_range(raster,wide_left,left);
            dim_t wide_right={b[0],b[1],mid,b[3]};
            val_t right={midval,vals[1]};
            color_cnt += checkerboard_range(raster,wide_right,right);
        }
    }
    return color_cnt;
}



/*! Creates a raster of the given size with the given range of values.
 *  ns is the total size of the raster, and vals is the range of values.
 *  Both are half-open intervals.
 */
template<class ARRAY, class VT>
void checkerboard_array(ARRAY& pmorph,
                        boost::array<size_t,2> ns,
                        boost::array<VT,2> limits)
{
    typedef boost::array<size_t,4> region_t;
    region_t whole = {0,ns[0],0,ns[1]};
    size_t color_cnt = checkerboard_range(pmorph, whole, limits);
    BOOST_ASSERT( color_cnt == limits[1] - limits[0] );
    if ( color_cnt < limits[1]-limits[0] ) {
        throw std::runtime_error("The requested matrix was too small "
                                 "to hold all the values.");
    }
}

#endif // _ARRAY_INIT_H_
