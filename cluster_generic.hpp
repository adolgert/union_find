#ifndef _CLUSTERS_GENERIC_H_
#define _CLUSTERS_GENERIC_H_

#include <functional>
#include <algorithm>
#include <set>
#include <cassert>
#include <boost/property_map/property_map.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/unordered_map.hpp>
#include <boost/array.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/utility.hpp>
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range2d.h"

#include "gather_clusters.hpp"


namespace raster_stats
{

	template<class Basis>
	class split_basis : public Basis {
		typedef typename Basis::size_type                size_type;
		typedef typename tbb::blocked_range2d<size_type> range_t;
		
		const Basis& basis_;
		range_t      range_;
	public:
		split_basis(const Basis& b, size_t granularity) : basis_(b),
            range_(0,b.size1(),granularity,0,b.size2(),granularity) {}
		split_basis(const split_basis& b)
            : basis_(b.basis_), range_(b.range_) {}
		~split_basis() {}

		split_basis( split_basis& r, tbb::split ) :
			basis_(r.basis_), range_(r.range_,tbb::split())
		{
		}

		bool empty() const { return range_.empty(); }
		bool is_divisible() const { return range_.is_divisible(); }
	};




    // The ublas matrix has the connectivity and data
    // The ublas matrix has a method that can construct a submatrix class
    // The tbb blocked_range2d can choose the size of the submatrix




    class adjacent_iterator :
        public boost::iterator_facade<adjacent_iterator,
            size_t const,boost::incrementable_traversal_tag>
    {
    public:
        //! bounds_type is (row start, row end, col start, col end)
        typedef boost::array<size_t,4> bounds_type;
        //! loc_type is (row,col)
        typedef boost::array<size_t,2> loc_type;
    private:
        const bounds_type& whole_;
        const bounds_type& bounds_;
        const loc_type& center_;
        value_type neighbor_;
        int direction_;
    public:
        adjacent_iterator(const bounds_type& whole, const bounds_type& bounds,
                          const loc_type& loc, int direction)
            : whole_(whole),bounds_(bounds), center_(loc),
              direction_(direction) {
            increment();
        }
        friend class boost::iterator_core_access;
        void increment() {
            while (direction_!=4) {
                direction_++;
                switch(direction_) {
                case 0:
                    if (center_[1]!=bounds_[3]-1) {
                        set_neighbor(center_[0],center_[1]+1);
                        return;
                    }
                    break;
                case 1:
                    if (center_[0]!=bounds_[1]-1) {
                        set_neighbor(center_[0]+1,center_[1]);
                        return;
                    }
                    break;
                case 2:
                    if (center_[1]!=bounds_[2]) {
                        set_neighbor(center_[0],center_[1]-1);
                        return;
                    }
                    break;
                case 3:
                    if (center_[0]!=bounds_[0]) {
                        set_neighbor(center_[0]-1,center_[1]);
                        return;
                    }
                    break;
                default:
                    assert(direction_<=4 && direction_>=0);
                }
            }
        }

        bool equal(adjacent_iterator const& other) const
        {
            return direction_==other.direction_;
        }

        reference dereference() const {
            return neighbor_;
        }
    private:
        void set_neighbor(value_type i, value_type j)
        {
            neighbor_=i*(whole_[3]-whole_[2])+j;
        }
    };




    /*! This iterates over a one-dimensional index into a subregion of a
     *  two-dimensional array. The index value is in the larger region,
     *  but differences and advancement act only within the subregion.
     *  You can make a non-const version of this too. See
     *  http://www.boost.org/doc/libs/1_48_0/libs/iterator/doc/iterator_facade.html
     */
    class array_iterator :
        public boost::iterator_facade<array_iterator,
                   size_t const,boost::random_access_traversal_tag>
    {
    public:
        //! bounds_type is (row start, row end, col start, col end)
        typedef boost::array<size_t,4> bounds_type;
        //! loc_type is (row,col)
        typedef boost::array<size_t,2> loc_type;
    private:
        const bounds_type whole_;
        const bounds_type bounds_;
        loc_type loc_;
        value_type index_;
    public:
        array_iterator(const bounds_type& whole,const bounds_type& bounds,
                       const loc_type& loc)
            : whole_(whole), bounds_(bounds),loc_(loc) {
            index_=(loc_[0]-whole_[0])*(whole_[3]-whole_[2])+(loc_[1]-whole_[2]);
        }
        friend class boost::iterator_core_access;

		boost::array<adjacent_iterator,2> adjacent() {
			return {{adjacent_iterator(whole_,bounds_,loc_,-1),
                        adjacent_iterator(whole_,bounds_,loc_,4)}};
		}

        void increment() {
            if (loc_[1]==bounds_[3]-1) {
                loc_[0]++;
                loc_[1]=bounds_[2];
                index_+=(whole_[3]-whole_[2])-(bounds_[3]-bounds_[2])+1;
            } else {
                loc_[1]++;
                index_++;
            }
        }
        void decrement() {
            if (loc_[1]==bounds_[2]) {
                loc_[1]=bounds_[3]-1;
                loc_[0]--;
                index_-=(whole_[3]-whole_[2])-(bounds_[3]-bounds_[2])+1;
            } else {
                loc_[1]--;
                index_--;
            }
        }

        //! The ptrdiff is within the inner region, not relative to whole.
        void advance(ptrdiff_t n) {
            // Calculate as [(i,j)-(j-j1)] + [n+(j-j1)]
            n+=(loc_[1]-bounds_[2]);
            ptrdiff_t rows=n/(bounds_[3]-bounds_[2]);
            ptrdiff_t cols=n%(bounds_[3]-bounds_[2]);
            loc_[0]+=rows;
            loc_[1]=bounds_[2]+cols;
            index_=(loc_[0]-whole_[0])*(whole_[3]-whole_[2])+
                (loc_[1]-whole_[2]);
        }

        //! distance_to is measured within the inner region, too.
        difference_type distance_to(const array_iterator& z) {
            difference_type rows=z.loc_[0]-loc_[0];
            difference_type cols=z.loc_[1]-loc_[1];
            return rows*(bounds_[3]-bounds_[2])+cols;
        }

        bool equal(array_iterator const& other) const
        {
            return loc_ == other.loc_;
        }

        reference dereference() const {
            reference val = index_;
            return val;
        }
    };



    struct ij {
        static size_t index(const boost::array<size_t,4>& bounds,
                            size_t i, size_t j)
        {
            return i*(bounds[3]-bounds[2])+j;
        }
    };



    class array_basis;

    /*! This iterates over all edges joining two regions of an array.
     *  The Vertex template argument might be an int or size_t.
     */
    class edge_iterator :
        public boost::iterator_facade<edge_iterator,
            boost::array<size_t,2> const,boost::incrementable_traversal_tag>
    {
    public:
        const boost::array<size_t,4>& whole_;
    private:
        value_type edge_;
        bool vertical_;
        boost::array<size_t,2> start_;
        size_t count_;
        size_t idx_;
    public:
        edge_iterator(const array_basis& a, const array_basis& b,bool end);
        friend class boost::iterator_core_access;

        void increment() {
            idx_++;
            if (vertical_) {
                edge_[0]=ij::index(whole_,start_[0]+idx_,start_[1]);
                edge_[1]=ij::index(whole_,start_[0]+idx_,start_[1]+1);
            } else {
                edge_[0]=ij::index(whole_,start_[0]  ,start_[1]+idx_);
                edge_[1]=ij::index(whole_,start_[0]+1,start_[1]+idx_);
            }
        }

        bool equal(edge_iterator const& b) const
        {
            return (vertical_==b.vertical_) && (start_==b.start_) &&
                (count_==b.count_) && (idx_==b.idx_);
        }

        reference dereference() const {
            return edge_;
        }

        boost::array<boost::array<size_t,2>,2> coords() {
            boost::array<boost::array<size_t,2>,2> here;
            if (vertical_) {
                here[0][0]=start_[0]+idx_;
                here[0][1]=start_[1];
                here[1][0]=start_[0]+idx_;
                here[1][1]=start_[1]+1;
            } else {
                here[0][0]=start_[0];
                here[0][1]=start_[1]+idx_;
                here[1][0]=start_[0]+1;
                here[1][1]=start_[1]+idx_;
            }
            return here;
        }
    };


    void print_bounds(const boost::array<size_t,4>& b) {
        std::cout << b[0] << ":" << b[1] << ":" << b[2] << ":" << b[3];
    }


    /*! This is a two-dimensional array with nearest-neighbors.
     *  This class brings together the 2D iterator with the
     *  splittable range concept from the TBB, so what you get
     *  is a 2D domain that can be broken into chunks.
     */
	class array_basis {
        typedef boost::array<size_t,4> bounds_type;
        typedef boost::array<size_t,2> loc_type;
		typedef size_t                 size_type;
		typedef typename tbb::blocked_range2d<size_type> range_t;
        typedef array_iterator iterator;
		
		range_t     range_;
        bounds_type whole_;
		bounds_type bounds_;
	public:
        friend class edge_iterator;

        array_basis(const bounds_type whole, size_t granularity) :
            whole_(whole),bounds_(whole),
            range_(whole[0],whole[1],granularity,whole[2],whole[3],granularity)
        {
        }
        array_basis(const array_basis& b) : whole_(b.whole_),bounds_(b.bounds_),
                                            range_(b.range_) {}
		~array_basis() {}

		array_basis( array_basis& r, tbb::split ) :
            whole_(r.whole_), range_(r.range_,tbb::split())
		{
            bounds_={{range_.rows().begin(),range_.rows().end(),
                            range_.cols().begin(),range_.cols().end() }};
			r.bounds_=	{{r.range_.rows().begin(),r.range_.rows().end(),
                          r.range_.cols().begin(),r.range_.cols().end() }};
            std::cout << "array_basis::array_basis: split ";
            std::cout << *this << "," << r << std::endl;
		}

		bool empty() const {
            bool b_empty=range_.empty();
            std::cout << "empty: "<< b_empty << " " << *this << std::endl;
            return b_empty;
        }

		bool is_divisible() const {
            bool b_divis=range_.is_divisible();
            std::cout << "divis: "<< b_divis << " " << *this << std::endl;
            return b_divis;
        }

        iterator begin() const {
            loc_type loc = {{bounds_[0],bounds_[2]}};
            return array_iterator(whole_,bounds_,loc);
        }

        iterator end() const {
            // The end isn't the upper right-hand corner but
            // a row above the upper left-hand because of what ++ does.
            loc_type loc = {{bounds_[1],bounds_[2]}};
            return array_iterator(whole_,bounds_,loc);
        }

        friend std::ostream& operator<<(std::ostream& os,const array_basis&a){
            os << a.bounds_[0] << ":" << a.bounds_[1] << ":"
               << a.bounds_[2] << ":" << a.bounds_[3] << " " <<
                a.range_.rows().begin() << ":" << a.range_.rows().end() <<
                ":" << a.range_.cols().begin() << ":" << a.range_.cols().end();
            return os;
        }
	};




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



    template<class Region,class Compare>
	class disjoint_set_cluster
	{
    public:
        //! Maps from element to count of elements in set.
        typedef boost::unordered_map<size_t,size_t>   rank_t;
        //! Maps from element to parent of element.
        typedef boost::unordered_map<size_t,size_t>   parent_t;
        typedef boost::associative_property_map<rank_t> rank_pmap_t;
        typedef boost::associative_property_map<parent_t> parent_pmap_t;
      
        typedef boost::disjoint_sets<rank_pmap_t,parent_pmap_t> dset_t;

        rank_t        rank_map_;
        parent_t      parent_map_;
        rank_pmap_t   rank_pmap_;
        parent_pmap_t parent_pmap_;
        dset_t        dset_;

        std::list<Region,std::allocator<Region>> seen_;

        Compare compare_;
    public:

        disjoint_set_cluster(Compare compare) : compare_(compare),
                rank_pmap_(rank_map_), parent_pmap_(parent_map_),
                dset_(rank_pmap_,parent_pmap_)
        {
        }

	    disjoint_set_cluster(disjoint_set_cluster& b,tbb::split)
            : compare_(b.compare_),
              rank_pmap_(rank_map_), parent_pmap_(parent_map_),
              dset_(rank_pmap_,parent_pmap_)
        {
            std::cout << "disjoint_set_cluster::split" << std::endl;
        }

        ~disjoint_set_cluster() { }

		/*! This is the reduction step where two instances of this class are joined.
         *  Called by tbb::internal::finish_reduce<Body>::execute() in
         *  parallel_reduce.h:62.
		 */
		void join(disjoint_set_cluster& b) {
            rank_map_.insert(b.rank_map_.begin(),b.rank_map_.end());
            parent_map_.insert(b.parent_map_.begin(),b.parent_map_.end());
            std::cout << "disjoint_set_cluster::join: enter" << std::endl;

            for (auto neighbor=b.seen_.begin();
                 neighbor!=b.seen_.end();
                 neighbor++) {
                 join_edges(*neighbor);
            }
		}

        void join_edges(const Region& neighbor) {
            for (auto local=seen_.begin();
                 local!=seen_.end();
                 local++) {
                std::cout << "join: a b" << std::endl;
                auto common=edge_iterator(neighbor,*local,false);
                auto common_end=edge_iterator(neighbor,*local,true);
                while (common!=common_end) {
                    std::cout << "join: edge" << std::endl;
                    auto& edge = *common;
                    if (compare_(edge[0],edge[1])) {
                        dset_.union_set(edge[0],edge[1]);
                    }
                    common++;
                }
            }
            seen_.push_back(neighbor);
        }

		/*! This acts on a subregion of the domain.
		 */
		void operator()(const Region& region) {
            std::cout << "disjoint_set_cluster::operator(): enter" << std::endl;
            auto vertex=region.begin();
            auto vertex_end=region.end();

            while (vertex!=vertex_end) {
                dset_.make_set(*vertex);
                boost::array<adjacent_iterator,2> adj = vertex.adjacent();
                while (adj[0]!=adj[1]) {
                    if (parent_map_.find(*adj[0])!=parent_map_.end()) {
                        if (compare_(*adj[0],*vertex)) {
                            dset_.union_set(*adj[0],*vertex);
                        }
                    }
                    ++adj[0];
                }
                ++vertex;
            }
            join_edges(region);
        }
	};
	
	
	template<class Landscape> class cluster_raster;
	
	template<class Landscape>
	size_t count(cluster_raster<Landscape>& cs);
	
	template<class Landscape>
	class cluster_raster {
		size_t result;
	public:
		cluster_raster() {}
		void operator()(const Landscape& raster) {
            boost::array<size_t,4> bounds =
                {{0,raster.size1(),0,raster.size2()}};
            array_basis gridlines(bounds,32);

            // Create a property map out of the ublas Matrix in
            // order to separate the basis from values defined on it.
            typedef typename Landscape::value_type value_type;
            boost::identity_property_map direct;
            boost::iterator_property_map<const value_type*,
                    boost::identity_property_map,value_type,const value_type&>
                land_use(&raster.data()[0],direct);
            
            typedef AreEqual<value_type,decltype(land_use)> AreEqual_t;
            AreEqual_t comparison(land_use);
			disjoint_set_cluster<array_basis,AreEqual_t>
                dsc(comparison);
            std::cout << "cluster_raster::operator(): Calling parallel_reduce"
                      << std::endl;
			tbb::parallel_reduce(gridlines,dsc);

            auto clusters = gather_clusters(dsc.parent_pmap_,
                                                 dsc.dset_,raster.size1(),
                                                 raster.size2());

			result = clusters->size();
		}
		friend size_t count<Landscape>(cluster_raster<Landscape>&);
	};


	
	
	/*! This is an extractor for cluster_raster to get results from it. */
	template<class Landscape>
	size_t count(cluster_raster<Landscape>& cs) {
		return cs.result;
	}
}


#endif // _CLUSTERS_GENERIC_H_
