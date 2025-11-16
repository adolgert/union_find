/*! cluster.cpp
 *  Given a landscape as a raster image, find clusters within it.
 */

#include <set>
#include <map>
#include <algorithm>
#include <iostream>
#include <vector>
#include <iterator>
#include <memory>
#include <boost/functional.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include "raster.hpp"
#include "cluster.hpp"


using namespace std;


namespace raster_stats {

	template<typename T>
	struct store_unique {
		set<T> vals;
		void operator()(const T& val) {
			this->vals.insert(val);
		};
	};



/*!
 * This uses an array to do an order-1 sort, since we are making a set
 * of unsigned char.
 */
set<arr_type> unique_values_direct(const landscape_t& raster) {
	typedef unsigned char arr_type;

	vector<int> slots(256,0);
	
	set<arr_type> uniques;
	for (auto iter2=raster.begin2(); iter2!=raster.end2(); ++iter2) {
		slots[*iter2]=1;
	}
	for (size_t slot_idx=0; slot_idx<slots.size(); slot_idx++) {
		if (0!=slots[slot_idx]) {
			uniques.insert(slot_idx);
		}
	}
	
	return uniques;
}



/*! This iterator walks through the i,j values of a 2d array.
 */
class iterator_i_j : public iterator<input_iterator_tag,int>
{
public:
	typedef std::pair<size_t,size_t> loc_t;
private:
	size_t ib,ie,jb,je;
	loc_t ij;
public:
	iterator_i_j(size_t i_min,size_t i_max, size_t j_min,size_t j_max)
		: ib(i_min), ie(i_max), jb(j_min), je(j_max), ij(ib,jb) {}
	iterator_i_j(size_t i_min,size_t i_max, size_t j_min,size_t j_max,size_t i, size_t j)
		: ib(i_min), ie(i_max), jb(j_min), je(j_max), ij(i,j) {}
	iterator_i_j(const iterator_i_j& it) : ib(it.ib),ie(it.ie),jb(it.jb),je(it.je),ij(it.ij) {}
	iterator_i_j& operator++() {ij.second++; if (ij.second==je) {ij.second=jb; ij.first++;}}
	iterator_i_j operator++(int) {iterator_i_j tmp(*this); operator++(); return tmp;}
	bool operator==(const iterator_i_j& rhs) {return ij==rhs.ij;}
	bool operator!=(const iterator_i_j& rhs) {return ij!=rhs.ij;}
	loc_t& operator*() {return ij;}
};



    /*! Find clusters in a raster. Uses (i,j) pairs to identify each location.
     *  This version uses std::pair(i,j) for each point in the raster.
     *  It assumes the input is a multiarray.
     */
cluster_loc_t find_clusters_pair(const landscape_t& raster)
{
	typedef pair<size_t,size_t> loc_t;
	typedef map<loc_t,size_t>   rank_t;
	typedef map<loc_t,loc_t>    parent_t;

	rank_t rank_map;
	boost::associative_property_map<rank_t>   rank_pmap(rank_map);
	parent_t parent_map;
	boost::associative_property_map<parent_t> parent_pmap(parent_map);

	boost::disjoint_sets<boost::associative_property_map<rank_t>,boost::associative_property_map<parent_t> >
	 	dset(rank_pmap,parent_pmap);

	// Every node gets to be its own set.
	for (size_t iset=0; iset<raster.size1(); iset++) {
		for (size_t jset=0; jset<raster.size2(); jset++) {
			dset.make_set(loc_t(iset,jset));
		}
	}

	// Then we connect neighbors with same values.
	for (size_t irow=0; irow<raster.size1()-1; irow++) {
		for (size_t jrow=0; jrow<raster.size2(); jrow++) {
			if (raster(irow,jrow)==raster(irow+1,jrow)) {
				dset.union_set(loc_t(irow,jrow),loc_t(irow+1,jrow));
			}
		}
	}

	for (size_t icol=0; icol<raster.size1(); icol++) {
		for (size_t jcol=0; jcol<raster.size2()-1; jcol++) {
			if (raster(icol,jcol)==raster(icol,jcol+1)) {
				dset.union_set(loc_t(icol,jcol),loc_t(icol,jcol+1));
			}
		}
	}

	// Pull out the sets.
	size_t imax=raster.size1();
	size_t jmax=raster.size2();

	dset.compress_sets(iterator_i_j(0,imax,0,jmax),
		iterator_i_j(0,imax,0,jmax,imax,0));

	cluster_loc_t clusters;
	for_each(iterator_i_j(0,imax,0,jmax),
			iterator_i_j(0,imax,0,jmax,imax,0),[&](const loc_t& sp) {
				clusters[dset.find_set(sp)].push_back(sp);
				});

	return clusters;
}




    /*! Find clusters, using size_t to identify each location.
     *  
     */
cluster_t find_clusters(const landscape_t& raster)
{
    //! Maps from element to count of elements in set.
	typedef map<size_t,size_t>   rank_t;
    //! Maps from element to parent of element.
	typedef map<size_t,size_t>   parent_t;

	rank_t rank_map;
	boost::associative_property_map<rank_t>   rank_pmap(rank_map);
	parent_t parent_map;
	boost::associative_property_map<parent_t> parent_pmap(parent_map);

	boost::disjoint_sets<boost::associative_property_map<rank_t>,
                         boost::associative_property_map<parent_t> >
	 	dset(rank_pmap,parent_pmap);

	size_t icnt=raster.size1();
	size_t jcnt=raster.size2();

	// Every node gets to be its own set.
	for (size_t make_set=0; make_set<icnt*jcnt; make_set++) {
		dset.make_set(make_set);
	}

	// Then we connect neighbors with same values.
	for (size_t i=0; i<icnt-1; i++) {
		for (size_t j=0; j<jcnt; j++) {
			if (raster(i,j)==raster(i+1,j)) {
				dset.union_set(i*jcnt+j,(i+1)*jcnt+j);
			}
		}
	}

	for (size_t i=0; i<icnt; i++) {
		for (size_t j=0; j<jcnt-1; j++) {
			if (raster(i,j)==raster(i,j+1)) {
				dset.union_set(i*jcnt+j,i*jcnt+j+1);
			}
		}
	}

	// At this point, each element points to a parent.
	// The return value is a list of lists of elements.
	// The temporary map will associate a parent with the list of its children.
	cluster_t clusters; // a list of lists

	typedef map<size_t,cluster_t::iterator> ptl_t;
	ptl_t parent_to_list; // The map from parent to list of children.

	for (size_t pull=0; pull<icnt*jcnt; pull++) {
		size_t parent = dset.find_set(pull);
		ptl_t::iterator plist = parent_to_list.find(parent);
		if (plist==parent_to_list.end()) {
			cluster_t::iterator nlist = clusters.insert(clusters.end(),list<size_t>());
			parent_to_list[parent]=nlist;
			nlist->push_back(pull);
		} else {
			plist->second->push_back(pull);
		}
	}
/*	cluster_t clusters;
	for (size_t pull=0; pull<icnt*jcnt; pull++) {
		clusters[dset.find_set(pull)].push_back(pull);
	}
*/	return clusters;
}






/*! Find clusters, using size_t to identify each location.
 *  Uses two passes in total.
 */
cluster_t find_clusters_twopass(const landscape_t& raster)
{
	typedef map<size_t,size_t>   rank_t; //! Maps from element to count of elements in set.
	typedef map<size_t,size_t>   parent_t; //! Maps from element to parent of element.

	rank_t rank_map;
	boost::associative_property_map<rank_t>   rank_pmap(rank_map);
	parent_t parent_map;
	boost::associative_property_map<parent_t> parent_pmap(parent_map);

	boost::disjoint_sets<boost::associative_property_map<rank_t>,boost::associative_property_map<parent_t> >
	 	dset(rank_pmap,parent_pmap);

	const size_t icnt=raster.size1();
	const size_t jcnt=raster.size2();

	for (size_t i=0; i<icnt; i++) {
		for (size_t j=0; j<jcnt; j++) {
			if (i==0 || j==0) {
				dset.make_set(i*jcnt+j);
			}

			if (i<icnt-1) {
				dset.make_set((i+1)*jcnt+j);
				if (raster(i,j)==raster(i+1,j)) {
					dset.union_set(i*jcnt+j,(i+1)*jcnt+j);
				}
			}

			if (j<jcnt-1) {
				dset.make_set(i*jcnt+j+1);
				if (raster(i,j)==raster(i,j+1)) {
					dset.union_set(i*jcnt+j,i*jcnt+j+1);
				}
			}
		}
	}

	// At this point, each element points to a parent.
	// The return value is a list of lists of elements.
	// The temporary map will associate a parent with the list of its children.
	cluster_t clusters; // a list of lists

	typedef map<size_t,cluster_t::iterator> ptl_t;
	ptl_t parent_to_list; // The map from parent to list of children.

	for (size_t pull=0; pull<icnt*jcnt; pull++) {
		size_t parent = dset.find_set(pull);
		ptl_t::iterator plist = parent_to_list.find(parent);
		if (plist==parent_to_list.end()) {
			cluster_t::iterator nlist = clusters.insert(clusters.end(),list<size_t>());
			parent_to_list[parent]=nlist;
			nlist->push_back(pull);
		} else {
			plist->second->push_back(pull);
		}
	}
	return clusters;
}





/*! The same as find_clusters_twopass, but returning a pointer.
 */
std::shared_ptr<cluster_t> find_clusters_pointer(const landscape_t& raster)
{
	typedef std::map<size_t,size_t>   rank_t; //! Maps from element to count of elements in set.
	typedef std::map<size_t,size_t>   parent_t; //! Maps from element to parent of element.

	rank_t rank_map;
	boost::associative_property_map<rank_t>   rank_pmap(rank_map);
	parent_t parent_map;
	boost::associative_property_map<parent_t> parent_pmap(parent_map);

	boost::disjoint_sets<boost::associative_property_map<rank_t>,boost::associative_property_map<parent_t> >
	 	dset(rank_pmap,parent_pmap);

	const size_t icnt=raster.size1();
	const size_t jcnt=raster.size2();

	for (size_t i=0; i<icnt; i++) {
		for (size_t j=0; j<jcnt; j++) {
			if (i==0 || j==0) {
				dset.make_set(i*jcnt+j);
			}

			if (i<icnt-1) {
				dset.make_set((i+1)*jcnt+j);
				if (raster(i,j)==raster(i+1,j)) {
					dset.union_set(i*jcnt+j,(i+1)*jcnt+j);
				}
			}

			if (j<jcnt-1) {
				dset.make_set(i*jcnt+j+1);
				if (raster(i,j)==raster(i,j+1)) {
					dset.union_set(i*jcnt+j,i*jcnt+j+1);
				}
			}
		}
	}

	return gather_clusters(parent_pmap,dset, icnt, jcnt);
}



//! Copied from boost::disjoint_sets::detail.
template <class ParentPA, class Vertex>
Vertex
find_representative_with_full_compression(ParentPA parent, Vertex v)
{
  Vertex old = v;
  Vertex ancestor = get(parent, v);
  while (ancestor != v) {
    v = ancestor;
    ancestor = get(parent, v);
  }
  v = get(parent, old);
  while (ancestor != v) {
    put(parent, old, ancestor);
    old = v;
    v = get(parent, old);
  }
  return ancestor;
}


/*! Find clusters, using size_t to identify each location.
 *  Uses two passes in total. Compare with twopass version. This
 *  version loops through the found clusters not by (i,j) but
 *  by going straight through the associative map.
 */
cluster_t find_clusters_remap(const landscape_t& raster)
{
	typedef map<size_t,size_t>   rank_t; //! Maps from element to count of elements in set.
	typedef map<size_t,size_t>   parent_t; //! Maps from element to parent of element.

	rank_t rank_map;
	boost::associative_property_map<rank_t>   rank_pmap(rank_map);
	parent_t parent_map;
	boost::associative_property_map<parent_t> parent_pmap(parent_map);

	boost::disjoint_sets<boost::associative_property_map<rank_t>,boost::associative_property_map<parent_t> >
	 	dset(rank_pmap,parent_pmap);

	const size_t icnt=raster.size1();
	const size_t jcnt=raster.size2();

	for (size_t i=0; i<icnt; i++) {
		for (size_t j=0; j<jcnt; j++) {
			if (i==0 || j==0) {
				dset.make_set(i*jcnt+j);
			}

			if (i<icnt-1) {
				dset.make_set((i+1)*jcnt+j);
				if (raster(i,j)==raster(i+1,j)) {
					dset.union_set(i*jcnt+j,(i+1)*jcnt+j);
				}
			}

			if (j<jcnt-1) {
				dset.make_set(i*jcnt+j+1);
				if (raster(i,j)==raster(i,j+1)) {
					dset.union_set(i*jcnt+j,i*jcnt+j+1);
				}
			}
		}
	}

	// At this point, each element points to a parent.
	// The return value is a list of lists of elements.
	// The temporary map will associate a parent with the list of its children.
	cluster_t clusters; // a list of lists

	typedef map<size_t,cluster_t::iterator> ptl_t;
	ptl_t parent_to_list; // The map from parent to list of children.

	for (auto read=parent_map.begin(); read!=parent_map.end(); read++) {
		size_t parent = read->second;
		if (parent!=read->first) {
			parent = find_representative_with_full_compression(parent_pmap,read->first);
		}
		ptl_t::iterator plist = parent_to_list.find(parent);
		if (plist==parent_to_list.end()) {
			cluster_t::iterator nlist = clusters.insert(clusters.end(),list<size_t>());
			parent_to_list[parent]=nlist;
			nlist->push_back(read->first);
		} else {
			plist->second->push_back(read->first);
		}
	}
	
	return clusters;
}


/*
find_clusters()
{
    avalanche_start = next_free();
    while (avalanche_start) {
        same_type = avalanche_start;
        avalanche;
        while (same_type) {
            mark_used(same_type);
            add_to_cluster(same_type,avalanche);
            same_type = next_same_type();
        }
        add_to_avalanches(avalanche);
        avalanche_start = next_free();
    }
    
    return blah;
}
*/

/*
template<typename T>
struct next_free_t {
	T& m_used;
	size_t m_i, m_j;
	size_t m_ci, m_cj;
	next_free_t(T& used) : m_used(used), m_ci(used.size1()), m_cj(used.size2()) { }
	pair<size_t,size_t> operator()() {
		while (m_used(m_i,m_j)!=true && m_j!=m_cj) m_increment();
		return pair(m_i,m_j);
	}
	void m_increment() {
		m_i++;
		if (m_i==m_ci) {
			m_i=0;
			m_j++;
		}
	}
};


template<typename T, typename U>
struct next_in_set_t {
	T& m_mat;
	U& m_used;
	pair<size_t,size_t> m_start;
	pair<size_t,size_t> m_loc;
	T::value_type m_representative;
	
	next_in_set_t(T& mat, U& used, pair<size_t,size_t> start) : m_mat(mat), m_used(used),m_start(start) {
		m_representative = m_mat(start.first,start.second);
		m_loc = start;
	}
	pair<size_t,size_t> operator()() {
		
	}
};




int find_clusters(const landscape_t& raster) {
	typedef unsigned char arr_type;
	boost::numeric::ublas::matrix<bool> used(raster.size1(),raster.size2());
	next_free_t next_free(used);
	
	typedef loc_t pair<size_t,size_t>;
	
	const loc_t used_end = pair(raster.size1(),raster.size2());

	loc_t avalanche_start = next_free();
	while (avalanche_start != used_end) {
		arr_type representative = raster(avalanche_start.first,avalanche_start.second);
		
		loc_t in_avalanche = avalanche_start;
		while (in_avalanche)
		
		avalanche_start = next_free();
	}
	// begin1 is a column-wise iterator.
	// begin2 is a row-wise iterator.
	for (landscape_t::const_iterator1 search=raster.begin1();
		search != raster.end1();
		search++) {
			uniques.insert(*search);
		}
	for_each( raster.begin1(), raster.end1(), [&](const arr_type& x) {
		uniques.insert(x); 
		});
	cout << "found uniques " << uniques.size() << endl;
	
	for_each(uniques.begin(),uniques.end(),[](const arr_type& x) {
		cout << (int)x << endl;
		});
	return 0;
}
*/


}
