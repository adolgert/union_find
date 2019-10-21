#ifndef _ARRAY_STORE_H_
#define _ARRAY_STORE_H_ 1

#include <vector>
#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/shared_ptr.hpp>
#include <morton.hpp>




namespace raster_stats {
    /*! This property map reorders access to underlying array.
     *  It is taken, whole-cloth, from boost::property_map, but this
     *  one accesses multi-arrays when given a token containing indices.
     *
     *  C is the container type. The value_type of this map comes from
     *  the array type.
     */
    template<typename T, typename C, typename IndexMap = boost::identity_property_map>
    class array_property_map
        : public boost::put_get_helper< 
              typename std::iterator_traits< 
                  typename C::iterator >::reference,
        array_property_map<T, C, IndexMap> >
    {
    public:
        typedef typename boost::property_traits<IndexMap>::key_type  key_type;
        typedef T value_type;
        typedef T& reference;
        typedef boost::lvalue_property_map_tag category;
        
        array_property_map(const IndexMap& index = IndexMap())
        : store(0), index(index)
        {}

        array_property_map(C& c,const IndexMap& index = IndexMap())
        : store(c), index(index)
        {}

        IndexMap&       get_index_map()       { return index; }
        const IndexMap& get_index_map() const { return index; }
          
    public:
        // Copy ctor absent, default semantics is OK.
        // Assignment operator absent, default semantics is OK.
        // CONSIDER: not sure that assignment to 'index' is correct.
        
        reference operator[](const key_type& v) const {
            typename boost::property_traits<IndexMap>::value_type i = get(index, v);
            //if (static_cast<unsigned>(i) >= store->size()) {
            //    store->resize(i + 1, T());
            //}
            return (*store)[i[0]][i[1]];
        }
    private:
        // Conceptually, we have a array of infinite size. For practical 
        // purposes, we start with an empty array and grow it as needed.
        // Note that we cannot store pointer to array here -- we cannot
        // store pointer to data, because if copy of property map resizes
        // the array, the pointer to data will be invalidated. 
        // I wonder if class 'pmap_ref' is simply needed.
        boost::shared_ptr<C> store;
        IndexMap index;
    };




    template <typename T>
    struct blocked_property_map
        : public boost::put_get_helper<T, boost::typed_identity_property_map<T> >
    {
        typedef T key_type;
        typedef T value_type;
        typedef T reference;
        typedef boost::readable_property_map_tag category;

        size_t block_;
        
        blocked_property_map(size_t block) : block_(block) {}

        inline value_type operator[](key_type v) const {
            size_t bi=v[0]/block_;
            size_t bj=v[1]/block_;
            size_t ni=v[0]%block_;
            size_t nj=v[1]%block_;
            return v;
        }
    };


    
    template<typename C,typename T, typename IndexMap>
    array_property_map<T, IndexMap>
    make_array_property_map(C* c, IndexMap index)
    {
        return array_property_map<C,T, IndexMap>(c,index);
    }



	template<typename T, typename TRANSFORM>
	class transform_helper {};



	/*! Transforms a key into a value, using a strategy pattern.
	 *  This doesn't do much, so why does it exist? It makes this translation of
	 *  coordinates look like a simple key-value map when it is really a 
	 *  key-value, whose value is the key for the next value.
 	 */
    template<typename TR, typename V>
    struct transform_map
        : public transform_helper<V, transform_map<TR,V> >
    {
        typedef TR                                   transform_type;
        typedef typename TR::key_type                key_type;
        typedef typename std::vector<V>              store_type;
        typedef typename store_type::value_type      value_type;
        typedef typename store_type::reference       reference;
        typedef typename store_type::const_reference const_reference;
        typedef boost::read_write_property_map_tag   category;
        
		TR tr_;
		store_type vec_;

        /*! alloc_n is the total number of elements to allocate.
         *  Unlike property_map types, this does not dynamically allocate.
         */
        transform_map(TR tr, size_t alloc_n) : tr_(tr), vec_(alloc_n) {}

        inline reference operator[](const key_type& v) {
			reference r=vec_[tr_(v)];
            return r;
        }

        inline const_reference operator[](const key_type& v) const {
            const_reference vv=vec_[tr_(v)];
            return vv;
        }
    };



	template <class TransformMap, class value_type, class K>
	inline const value_type&
	get(const transform_helper<value_type,TransformMap>& pa, const K& k)
	{
		return static_cast<const TransformMap&>(pa)[k];
	}



	/*! Turns i,j into a single index n, without funny-business.
	 */
	struct transform_ij {
		typedef boost::array<size_t,2> key_type;
		typedef size_t value_type;
		size_t w_;
		
		transform_ij(size_t w,size_t block=0) : w_(w) {}
		
		value_type operator()(const key_type& k) const {
			return k[0]*w_+k[1];
		}
	};




	/*! Width-agnostic blocking. Transposes each set of b^2 elements.
	 */
	struct transform_ij_blocked {
		typedef boost::array<size_t,2> key_type;
		typedef size_t value_type;
		size_t b_;
		size_t w_;
		
		transform_ij_blocked(size_t w, size_t b) : b_(b), w_(w) {}
		
		value_type operator()(const key_type& k) const {
			size_t n=k[0]*w_+k[1];
			size_t nd=(n%(b_*b_));
			return (n-nd) + (nd/b_) + (nd%b_)*b_;			
		}
	};


	// boustrophedon?


	/*! Blocking, the way we normally think of it.
	 */
	struct transform_ij_full_blocked {
		typedef boost::array<size_t,2> key_type;
		typedef size_t value_type;
		size_t b_;  //! Block size.
		size_t w_;  //! Width of j.
		size_t wb_; //! Width in blocks.
		
		transform_ij_full_blocked(size_t w, size_t b) : b_(b), w_(w), wb_((w_-1)/b_+1) {}
		
		value_type operator()(const key_type& k) const {
			return ((k[0]/b_)*wb_+k[1]/b_)*b_*b_ + (k[0]%b_)*b_+(k[1]%b_);
		}
	};



	struct transform_morton_ij {
	    typedef boost::array<size_t,2> key_type;
        typedef size_t value_type;

        transform_morton_ij() {}

        value_type operator()(const key_type& k) const {
            return morton_calculations::combine_xy<typename key_type::value_type,value_type>(k);
        }
	};



    struct transform_morton {
        typedef size_t key_type;
        typedef size_t value_type;

        transform_morton() {}

        value_type operator()(const key_type& k) const {
            return k;
        }
    };
}


#endif // _ARRAY_STORE_H_
