#ifndef _SINGLE_HPP_
#define _SINGLE_HPP_ 1

#include <functional>
#include <memory>
#include "boost/concept/assert.hpp"
#include "boost/property_map/property_map.hpp"
#include "boost/pending/disjoint_sets.hpp"



namespace raster_stats {


    /*! This binary comparison uses an associative property map
     *  to determine whether the two vertices are equal.
     */
    template<class Vertex,class Property>
    struct AreEqual :
        public std::binary_function<Vertex,Vertex,void>
    {
        const Property& p_;
        AreEqual(const Property& p) : p_(p) {}
        AreEqual(const AreEqual& b) : p_(b.p_) {}
        bool operator()(const Vertex& a, const Vertex& b) const {
            return get(p_,a)==get(p_,b);
        }
    };



    template<class VT,template<typename T,typename A> class RANK,
             template<typename B,typename C> class PARENT>
    struct construct_disjoint_set
    {
        typedef RANK<VT,size_t> rank_t;
        typedef PARENT<VT,VT> parent_t;
        typedef boost::associative_property_map<rank_t> rank_pmap_t;
        typedef boost::associative_property_map<parent_t> parent_pmap_t;
        typedef boost::disjoint_sets<rank_pmap_t,parent_pmap_t> dset_t;

        rank_t        rank_map_;
        parent_t      parent_map_;
        rank_pmap_t   rank_pmap_;
        parent_pmap_t parent_pmap_;
        dset_t        dset_;

        construct_disjoint_set()
            : rank_pmap_(rank_map_), parent_pmap_(parent_map_),
              dset_(rank_pmap_,parent_pmap_) {}
    };


    

    /*! A serial union-find with lots of parameterization.
     *  BASIS is the grid on which the values live.
     *  ADJACENCY finds neighbors of a vertex.
     *  TEST compares two vertices to determine whether they are in the set.
     *  CONSTRUCT chooses how to build the disjoint_set. It's a Policy.
     */
    template<class CONSTRUCT>
	class union_find_st : public CONSTRUCT
	{
    public:

        template<class BASIS, class TEST, class ITER, class ADJACENCY>
		void operator()(const BASIS& basis, const TEST& compare,
                        ITER vertex_iterator, ADJACENCY adjacent)
    	{
	        auto vertices=vertex_iterator(basis);
	        auto vertex_idx=vertices[0];
	        auto vertex_cnt=vertices[1];

            size_t unioned=0;
            size_t total=0;
	        while (vertex_idx!=vertex_cnt) {
	            this->dset_.make_set(*vertex_idx);
	            auto adj = adjacent(basis,*vertex_idx);
                bool b_union=false;
	            while ( adj[0] != adj[1] ) {
	                if (this->parent_map_.find(*adj[0])!=
                            this->parent_map_.end()) {
	                    if (compare(*adj[0],*vertex_idx)) {
	                        this->dset_.union_set(*adj[0],*vertex_idx);
                            b_union=true;
	                    }
                    }
                    
	                ++adj[0];
	            }
                if (b_union) ++unioned;
                ++total;
	            ++vertex_idx;
	        }
            std::cout << "Number unioned " << unioned << std::endl;
            std::cout << "Not unioned " << total-unioned << std::endl;
	    }
	};

}


#endif // _SINGLE_HPP_
