#ifndef _CLUSTER_G_
#define _CLUSTER_G_ 1
#include <boost/tuple.hpp>
#include <boost/array.hpp>
#include "tbb/parallel_reduce.h"



namespace raster_stats
{

class tbb_range_archetype {
public:
	tbb_range_archetype(const tbb_range_archetype& b) {}
	~tbb_range_archetype() {}
	
	tbb_range_archetype( tbb_range_archetype& r, split ) {}
	
	bool empty() const { return true; }
	bool is_divisible() const { return false; }
};



/*! This helper base class ensures the global free functions apply only to the
 *  classes we designate.
 */
template<class VertexType>
struct provides_adjacent {}

template<class Iterator,class VertexType>
  std::tuple<Iterator,Iterator> adjacent_vertices(provides_adjacent<VertexType>& r, VertexType& v)
  { return r.adjacent(v); }


template<class VertexType>
class region_archetype :
 public provides_adjacent<VertexType>,
 public tbb_range_archetype
{
      std::vector<VertexType> vertex_;
public:
	typedef VertexType vertex_type;
    typedef std::vector<VertexType>::iterator iterator;
    typedef std::tuple<iterator,iterator> iterator_pair;
	
	region_archetype() {}
	~region_archetype() {}
	
	iterator begin() { return vertex_.begin(); }
    iterator end() { return vertex_.end() };

    iterator_pair adjacent(const vertex& v) { return iterator_pair(begin(),end()); }

};



    template<class Vertex,class Property>
      struct AreEqual
      : std::binary_function(Vertex,Vertex,void)
    {
      Property& p_;
    AreEqual(Property& p) : p_(p) {}
      bool operator()(Vertex& a, Vertex& b) {
        return get(property,a)==get(property,b);
      }
    };


    template<class RankPA,class ParentPA,
      class FindCompress=find_with_full_path_compression>
      class joinable_disjoint_sets :
    public disjoint_sets<RankPA,ParentPA,FindCompress>
      {
      public:
        void join(const joinable_disjoint_sets& b) {
          rank.insert(b.rank->begin(),b.rank->end());
          rep.insert(b.rep->begin(),b.rep->end());
        }
      };



/*! This function iterates through all vertices, calling a function on every
 *  pair of vertices that have been seen.
 */
template<class Region,class NeighborAction>
  void walk_seen_neighbors(Region r, NeighborAction f)
{
  typedef typename Region::vertex_type Vertex;
  std::set<Vertex> seen;
  for( auto walk=begin(r); walk!=end(r), ++walk ) {
    tie(neighbor_begin,neighbor_end) = adjacent_vertices(r,*walk);
    for (auto neighbor=neighbor_begin; neighbor!=neighbor_end; ++neighbor) {
        if (seen.find(n)!=seen.end()) {
          f(*walk,*neighbor);
        }
      }
    seen.insert(*walk);
  }
}



template<class Region,class EdgeFunction>
void walk_boundary(Region a, Region b, EdgeFunction f)
{
    // connecting_edges should get edges starting in a
    // and finishing in b.
	tie(edge_begin,edge_end) = connecting_edges(a,b);
	for (auto edge=edge_begin; edge!=edge_end; ++edge) {
        f(start_vertex(a,edge),finish_vertex(b,edge));
    }
}



template<class DisjoinSet,class Compare,class Vertex>
class union_if_equal
: std::binary_function(Vertex,Vertex,void)
{
	DisjointSet& dset_;
	Compare& comp_;
	
	union_if_equal(DisjointSet& dset,Compare comp) : dset_(dset),comp_(comp) {}
	void operator()(Vertex a, Vertex b) {
		if (comp(a,b)) {
            dset_->union_set(a,b);
		}
	}
}



/*! Clustering algorithm using disjoint sets.
 *  When we write out this algorithm, it looks like a strategy
 *  passed to the parallel reduce algorithm.
 *
 *  The DisjointSet template argument is the type of a particular instantiation
 *  of disjoint_set, together with a join() method that will unify the maps of
 *  two disjoint sets.
 *
 *  The Region template argument is only the support, not the data on that
 *  support. Only CompareVertices needs to know the data values.
 *  The Region concept for this class includes both how to walk the data
 *  and how to split it into parts. TBB separates splitting as a strategy.
 */
template<class Region,class DisjointSet,class CompareVertices>
class disjoint_set_cluster
{
    typedef typename Region::vertex_type vertex_type;
    typedef union_if_equal<DisjointSet,CompareVertices,vertex_type> union_type;

	DisjointSet& dset_;
    union_type& union_if;
	std::set<Range> seen_;

 public:

    disjoint_set_cluster(DisjointSet& dset, CompareVertices comp_)
      : dset_(dset), union_if_equal(dset_,comp_)
      {
      }

    disjoint_set_cluster(DisjointSetCluster& b,split)
      : dset_(new DisjointSet), union_if_equal(b.union_if_equal)
      {
      }

	/*! This is the reduction step where two instances of this class are joined.
	 */
	void join(const disjoint_set_cluster& b) {
		dset_.join(b.dset_);

		for(auto region=begin(b.seen_); region!= end(b.seen_); ++region) {
          tie(ar_begin,ar_end)=adjacent_regions(*region);
          for (auto neighbor=ar_begin; neighbor!=ar_end; ++neighbor) {
            if (seen_.find(*neighbor)!=seen.end()) {
              walk_boundary(*region,*neighbor,union_if);
            }
          }
        }

        for(auto add=begin(b.seen_); add!=end(b.seen_);++add) {
          seen_.insert(*add);
        }
	}

	/*! This acts on a subregion of the domain.
	 */
	void operator()(Region region) {
		walk_seen_neighbors(region,union_if_equal(dset_));
        tie(ar_begin,ar_end) = adjacent_regions(region);
		for (auto adjacent=ar_begin; adjacent!=ar_end; ++adjacent) {
          if (seen_.find(*adjacent)!=seen_.end()) {
            walk_boundary(region,*adjacent,union_if);
          }
        }
		seen_.insert(region);
	}
	
	
	
	
};



template<class Landscape>
class cluster_raster {
	cluster_raster() {}
	void operator()(const Landscape& raster) {
		auto cs=disjoint_set_cluster(raster);
	    tbb::parallel_reduce( blocked_range2d<Landscape::size_type>(
	                           0,raster.size1(),32,
	                           0,raster.size2(),32),
	                  cs
	                  );

	    return gather_clusters(*(cs.m_parent_pmap), *(cs.m_dset), raster.size1(),
	                           raster.size2());
	}    
};


}

#endif // _CLUSTER_G_
