#ifndef _GRID2D_H_
#define _GRID2D_H_

#include <boost/array.hpp>
#include "tbb/blocked_range2d.h"


namespace raster_stats {


    template<class VT,class BT>
    class four_adjacent_iterator :
        public boost::iterator_facade<four_adjacent_iterator<VT,BT>,
            VT const,boost::incrementable_traversal_tag>
    {
    public:
        //! bounds_type is (row start, row end, col start, col end)
        typedef BT bounds_type;
        typedef VT vertex_type;
    private:
        const bounds_type& whole_;
        const bounds_type& bounds_;
        const vertex_type& center_;
        vertex_type neighbor_;
        int direction_;
    public:
        four_adjacent_iterator(const bounds_type& whole,
                               const bounds_type& bounds,
                               const vertex_type& loc, int direction)
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
                    if (center_[1]!=bounds_[1][1]-1) {
                        set_neighbor(center_[0],center_[1]+1);
                        return;
                    }
                    break;
                case 1:
                    if (center_[0]!=bounds_[0][1]-1) {
                        set_neighbor(center_[0]+1,center_[1]);
                        return;
                    }
                    break;
                case 2:
                    if (center_[1]!=bounds_[1][0]) {
                        set_neighbor(center_[0],center_[1]-1);
                        return;
                    }
                    break;
                case 3:
                    if (center_[0]!=bounds_[0][0]) {
                        set_neighbor(center_[0]-1,center_[1]);
                        return;
                    }
                    break;
                default:
                    assert(direction_<=4 && direction_>=0);
                }
            }
        }

        bool equal(four_adjacent_iterator const& other) const
        {
            return direction_==other.direction_;
        }

        vertex_type const& dereference() const {
            return neighbor_;
        }
    private:
        void set_neighbor( typename vertex_type::value_type i,
                           typename vertex_type::value_type j)
        {
            neighbor_[0]=i;
            neighbor_[1]=j;
        }
    };




     /*! This walks from 0 to N.
     *  It is just a simple walk. The key is something with a ++operator.
     *  The full implementation should walk a sub-region of the array,
     *  but this walks the whole array.
     */
    template<class BASIS>
    class linear_iterator :
        public boost::iterator_facade<linear_iterator<BASIS>,
           typename BASIS::vertex_type const,
           boost::random_access_traversal_tag>
    {
    public:
        //! bounds_type is (row start, row end, col start, col end)
        typedef typename BASIS::bounds_type bounds_type;
        typedef typename BASIS::vertex_type vertex_type;
    private:
        bounds_type whole_;
        bounds_type bounds_;
        vertex_type loc_;
    public:
        linear_iterator() {}
        linear_iterator(const bounds_type& whole,const bounds_type& bounds,
                       const vertex_type& loc)
            : whole_(whole), bounds_(bounds),loc_(loc) {
        }
        friend class boost::iterator_core_access;

        linear_iterator(const linear_iterator& b)
            :  whole_(b.whole_), bounds_(b.bounds_),loc_(b.loc_) {}

        linear_iterator& operator=(const linear_iterator& b) {
            whole_=b.whole_;
            bounds_=b.bounds_;
            loc_=b.loc_;
            return *this;
        }

        void increment() {
            ++loc_;
        }
        void decrement() {
            --loc_;
        }

        //! The ptrdiff is within the inner region, not relative to whole.
        void advance(ptrdiff_t n) {
            loc_+=n;
        }

        //! distance_to is measured within the inner region, too.
        size_t distance_to(const linear_iterator& z) {
            return z-loc_;
        }

        bool equal(linear_iterator const& other) const
        {
            return loc_ == other.loc_;
        }

        vertex_type const& dereference() const {
            vertex_type const& val = loc_;
            return val;
        }
    };




    /*! This iterates over a one-dimensional index into a subregion of a
     *  two-dimensional array. The index value is in the larger region,
     *  but differences and advancement act only within the subregion.
     *  You can make a non-const version of this too. See
     *  http://www.boost.org/doc/libs/1_48_0/libs/iterator/doc/iterator_facade.html
     */
    template<class BASIS>
    class array_iterator :
        public boost::iterator_facade<array_iterator<BASIS>,
           typename BASIS::vertex_type const,
           boost::random_access_traversal_tag>
    {
    public:
        //! bounds_type is (row start, row end, col start, col end)
        typedef typename BASIS::bounds_type bounds_type;
        typedef typename BASIS::vertex_type vertex_type;
    private:
        bounds_type whole_;
        bounds_type bounds_;
        vertex_type loc_;
    public:
        array_iterator() {}
        array_iterator(const bounds_type& whole,const bounds_type& bounds,
                       const vertex_type& loc)
            : whole_(whole), bounds_(bounds),loc_(loc) {
        }
        friend class boost::iterator_core_access;

        array_iterator(const array_iterator& b)
            :  whole_(b.whole_), bounds_(b.bounds_),loc_(b.loc_) {}

        array_iterator& operator=(const array_iterator& b) {
            whole_=b.whole_;
            bounds_=b.bounds_;
            loc_=b.loc_;
            return *this;
        }

        void increment() {
            if (loc_[1]==bounds_[1][1]-1) {
                loc_[0]++;
                loc_[1]=bounds_[1][0];
            } else {
                loc_[1]++;
            }
        }
        void decrement() {
            if (loc_[1]==bounds_[1][0]) {
                loc_[1]=bounds_[1][1]-1;
                loc_[0]--;
            } else {
                loc_[1]--;
            }
        }

        //! The ptrdiff is within the inner region, not relative to whole.
        void advance(ptrdiff_t n) {
            // Calculate as [(i,j)-(j-j1)] + [n+(j-j1)]
            n+=(loc_[1]-bounds_[1][0]);
            ptrdiff_t rows=n/(bounds_[1][1]-bounds_[1][0]);
            ptrdiff_t cols=n%(bounds_[1][1]-bounds_[1][0]);
            loc_[0]+=rows;
            loc_[1]=bounds_[1][0]+cols;
        }

        //! distance_to is measured within the inner region, too.
        size_t distance_to(const array_iterator& z) {
            size_t rows=z.loc_[0]-loc_[0];
            size_t cols=z.loc_[1]-loc_[1];
            return rows*(bounds_[3]-bounds_[2])+cols;
        }

        bool equal(array_iterator const& other) const
        {
            return loc_ == other.loc_;
        }

        vertex_type const& dereference() const {
            vertex_type const& val = loc_;
            return val;
        }
    };


    template<class BASIS>
    boost::array<four_adjacent_iterator<typename BASIS::vertex_type,
                                        typename BASIS::bounds_type>,2>
    make_four_adjacent(const BASIS& basis,
                       const typename BASIS::vertex_type& loc_)
    {
        return {{
                four_adjacent_iterator<typename BASIS::vertex_type,
                                       typename BASIS::bounds_type>
                    (basis.whole_,basis.bounds_,loc_,-1),
                    four_adjacent_iterator<typename BASIS::vertex_type,
                                           typename BASIS::bounds_type>
                    (basis.whole_,basis.bounds_,loc_,4)}};
    }



    /*! This is a two-dimensional array with nearest-neighbors.
     *  This class brings together the 2D iterator with the
     *  splittable range concept from the TBB, so what you get
     *  is a 2D domain that can be broken into chunks.
     *  VT is the storage type of the vertex, size_t or unsigned.
     */
    template<class VT>
	class array_basis {
    public:
        typedef boost::array<VT,2>          vertex_type;
        typedef boost::array<vertex_type,2> bounds_type;
		typedef VT                          size_type;

		typedef typename tbb::blocked_range2d<size_type> range_t;
        typedef array_iterator<array_basis> iterator;
		
		range_t     range_;
    public:
        bounds_type whole_;
		bounds_type bounds_;
	public:

        array_basis(const bounds_type whole, size_t granularity) :
            whole_(whole),bounds_(whole),
            range_(whole[0][0],whole[0][1],granularity,
                   whole[1][0],whole[1][1],granularity)
        {
        }
        array_basis(const array_basis& b) : whole_(b.whole_),bounds_(b.bounds_),
                                            range_(b.range_) {}
		~array_basis() {}

		array_basis( array_basis& r, tbb::split ) :
            whole_(r.whole_), range_(r.range_,tbb::split())
		{
            bounds_[0][0]=range_.rows().begin();
            bounds_[0][1]=range_.rows().end();
            bounds_[1][0]=range_.cols().begin();
            bounds_[1][1]=range_.cols().end();
            r.bounds_[0][0]=r.range_.rows().begin();
            r.bounds_[0][1]=r.range_.rows().end();
            r.bounds_[1][0]=r.range_.cols().begin();
            r.bounds_[1][1]=r.range_.cols().end();
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
            vertex_type loc = {{bounds_[0],bounds_[2]}};
            return array_iterator<array_basis>(whole_,bounds_,loc);
        }

        iterator end() const {
            // The end isn't the upper right-hand corner but
            // a row above the upper left-hand because of what ++ does.
            vertex_type loc = {{bounds_[1],bounds_[2]}};
            return array_iterator<array_basis>(whole_,bounds_,loc);
        }

        friend std::ostream& operator<<(std::ostream& os,const array_basis&a){
            os << a.bounds_[0][0] << ":" << a.bounds_[0][1] << ":"
               << a.bounds_[1][0] << ":" << a.bounds_[1][1] << " " <<
                a.range_.rows().begin() << ":" << a.range_.rows().end() <<
                ":" << a.range_.cols().begin() << ":" << a.range_.cols().end();
            return os;
        }
	};



    template<class BASIS> 
    boost::array<typename BASIS::iterator,2>
    make_vertex_iterator(const BASIS& basis) {
        typename BASIS::vertex_type vertex_ll=
            {{basis.bounds_[0][0],basis.bounds_[1][0]}};
        typename BASIS::vertex_type vertex_rr=
            {{basis.bounds_[0][1],basis.bounds_[1][0]}};
        boost::array<typename BASIS::iterator,2> arr;
        arr[0]=typename BASIS::iterator(basis.whole_,basis.bounds_,vertex_ll);
        arr[1]=typename BASIS::iterator(basis.whole_,basis.bounds_,vertex_rr);
        return arr;
    }


}

#endif // _GRID2D_H_

