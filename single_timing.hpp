#ifndef _SINGLE_TIMING_HPP_
#define _SINGLE_TIMING_HPP_ 1

#include <memory>
#include <map>
#include <boost/assert.hpp>
#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>
#include "array_store.hpp"
#include "grid2d.hpp"
#include "single.hpp"
#include "array_init.hpp"
#include "gather_clusters.hpp"


namespace raster_stats
{    
    
    template<class BASIS,class MAP>
    boost::tuple<std::shared_ptr<BASIS>,std::shared_ptr<MAP>>
        make_data(size_t width, size_t height, size_t block, size_t level_cnt)
        {
            typedef typename MAP::value_type value_type;
            boost::array<size_t,2> extent;
            typename BASIS::bounds_type bounds;
            boost::array<value_type,2> limits;
            extent[0]=width;
            extent[1]=height;

            bounds[0][0]=0;
            bounds[0][1]=extent[0];
            bounds[1][0]=0;
            bounds[1][1]=extent[1];
            auto basis=std::make_shared<BASIS>(bounds,32);

            auto transform=std::make_shared<typename MAP::transform_type>(
                                                         extent[0],block);
            auto data=std::make_shared<MAP>(*transform,extent[0]*extent[1]);

            limits[0]=0;
            limits[1]=level_cnt;
            checkerboard_array(*data,extent,limits);
            return boost::make_tuple<std::shared_ptr<BASIS>,
                                std::shared_ptr<MAP>>(basis,data);
        }



        // std::map has third and fourth default template arguments
        // that this hides.
        template<typename FROM, typename TO>
        struct bound_map : std::map<FROM,TO> {};


    template<class BASIS,class DATA>
    class single_run
    {
        std::shared_ptr<BASIS> basis_;
        std::shared_ptr<DATA> data_;
    public:
        single_run(std::shared_ptr<BASIS> basis,
            std::shared_ptr<DATA> data) : basis_(basis),data_(data) {}

        void operator()() {

            typedef AreEqual<typename DATA::key_type,DATA> comparison_t;
            comparison_t comparison(*data_);


            typedef construct_disjoint_set<typename BASIS::vertex_type,
                                           bound_map,
                                           bound_map> disj_t;

            union_find_st<disj_t> ufind;
            ufind(*basis_,comparison,
                  make_vertex_iterator<BASIS>,
                  make_four_adjacent<BASIS>);
    
            //auto clusters = gather_clusters(ufind.rank_pmap_,ufind.dset_,
            //                                extent);
            //BOOST_ASSERT(clusters->size(),limits[1]-limits[0]);
        }
    };
}

#endif // _SINGLE_TIMING_HPP_
