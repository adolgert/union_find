/*! cluster_tbb.cpp
 *  This file uses the TBB to do clustering in parallel.
 */

#include <set>
#include <map>
#include <algorithm>
#include <iostream>
#include <vector>
#include <iterator>
#include <memory>
#include <boost/unordered_map.hpp>
#include <boost/array.hpp>
#include <boost/functional.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range2d.h"

#include "raster.hpp"
#include "cluster.hpp"

using namespace tbb;
using namespace std;


namespace raster_stats {


// Prints arrays of four size_t to stdout.
  ostream& operator<<(ostream& s, const boost::array<size_t,4>& arr) {
    return s << arr[0] << " " << arr[1] << " " << arr[2] << " " << arr[3];
}




/*! Connects elements from a grid into sets.
 *  An object of this class is passed to parallel_for
 *  so that it can work on a smaller region.
 */
struct ConnectSets
{
    //! Maps from element to count of elements in set.
    typedef boost::unordered_map<size_t,size_t>   rank_t;
    //! Maps from element to parent of element.
    typedef boost::unordered_map<size_t,size_t>   parent_t;
    typedef boost::associative_property_map<rank_t> rank_pmap_t;
    typedef boost::associative_property_map<parent_t> parent_pmap_t;

    typedef boost::disjoint_sets<rank_pmap_t,parent_pmap_t> dset_t;

    const landscape_t& m_raster;
    std::shared_ptr<rank_t> m_rank_map;
    std::shared_ptr<parent_t> m_parent_map;
    std::shared_ptr<rank_pmap_t> m_rank_pmap;
    std::shared_ptr<parent_pmap_t> m_parent_pmap;
    std::shared_ptr<dset_t> m_dset;
    boost::array<size_t,4> m_range;

    size_t m_row_cnt;
    
    typedef boost::array<size_t,2> coord_t;
    typedef map<coord_t,size_t> edge_t;
    edge_t m_rows;
    edge_t m_cols;

    ConnectSets(const landscape_t& raster) : m_raster(raster) {
        this->create_dset();
    }

    ~ConnectSets() {
    }

    /*! Splitting constructor for TBB to create another thread. */
    ConnectSets(ConnectSets& b, split) : m_raster(b.m_raster) {
        this->create_dset();
    }

    void create_dset() {
        m_range.fill(0);
        m_rank_map = std::make_shared<rank_t>();
        m_parent_map = std::make_shared<parent_t>();
        m_rank_pmap = std::make_shared<rank_pmap_t>(*m_rank_map);
        m_parent_pmap = std::make_shared<parent_pmap_t>(*m_parent_map);
        m_dset = std::make_shared<dset_t>(*m_rank_pmap,*m_parent_pmap);

        m_row_cnt = m_raster.size2();
    }

    void join( const ConnectSets& b ) {
        //cout << "ConnectSets::join " << m_range << " to " <<
        //    b.m_range << endl;
        m_rank_map->insert(b.m_rank_map->begin(), b.m_rank_map->end());
        m_parent_map->insert(b.m_parent_map->begin(), b.m_parent_map->end());
        for_each(b.m_rows.begin(), b.m_rows.end(),
                 [&](const edge_t::value_type& val) {
                     this->add_row(val); } );
        for_each(b.m_rows.begin(), b.m_rows.end(),
                 [&](const edge_t::value_type& val) {
                     this->add_col(val); } );
    }

    /*! This adds all four edges of the region to a list of rows
     *  and columns whose regions have been added to the disjoint
     *  set. If these edges were seen before, then it is time to
     *  look for connections across the edges. This ASSUMES that
     *  each edge is the same length, so we can identify them just
     *  by their starting points.
     */
    void add_edges(const blocked_range<size_t>& rows,
                   const blocked_range<size_t>& cols) {
        coord_t lower_left  = {{ rows.begin(), cols.begin() }};
        coord_t upper_left  = {{ rows.end(),   cols.begin() }};
        coord_t lower_right = {{ rows.begin(), cols.end() }};

        // Both top and bottom rows start at the cols.begin()
        // but one is situated at rows.begin(), the other at rows.end().
        auto bottom = edge_t::value_type(lower_left,cols.end());
        auto top = edge_t::value_type(upper_left,cols.end());
        this->add_row(bottom);
        this->add_row(top);

        auto left = edge_t::value_type(lower_left,rows.end());
        auto right = edge_t::value_type(lower_right,rows.end());
        this->add_col(left);
        this->add_col(right);
    }


    void add_row(const edge_t::value_type& row_entry) {
        const auto& row = row_entry.first;
        size_t end = row_entry.second;

        auto found = m_rows.find(row);
        if (found == m_rows.end()) {
            m_rows.insert(row_entry);
        } else {
            //cout << "Ready to zip row " << row[0] << " " <<
            //    row[1] << " up to " << end << endl;
            const size_t jcnt=m_raster.size2();
            size_t i = row[0];
            for (size_t j = row[1]; j<end; j++) {
                union_if_equal(i,j,i-1,j);
            }
            m_rows.erase(found);
        }
    }



    void add_col(const edge_t::value_type& col_entry) {
        const auto& col = col_entry.first;
        size_t end = col_entry.second;

        auto found = m_cols.find(col);
        if (found == m_cols.end()) {
            m_cols.insert(col_entry);
        } else {
            //cout << "Ready to zip col " << col[0] << " " <<
            //    col[1] << " up to " << end << endl;
            const size_t jcnt=m_raster.size2();
            size_t j = col[1];
            for (size_t i = col[0]; i<end; i++) {
                union_if_equal(i,j,i,j-1);
            }
            m_cols.erase(found);
        }
    }



    void operator() (const blocked_range2d<landscape_t::size_type>& r ) {
        //cout << "ConnectSets::operator() " << m_range << " to ";
        boost::array<size_t,4> range_init = {{ r.rows().begin(),r.rows().end(),
                r.cols().begin(), r.cols().end() }};
        m_range=range_init;
        //cout << m_range << endl;
        const size_t jcnt=m_raster.size2();

        m_dset->make_set(r.rows().begin()*jcnt+r.cols().begin());

        for (size_t fr_idx=r.rows().begin()+1; fr_idx<r.rows().end();
                                                        fr_idx++) {
            m_dset->make_set(fr_idx*jcnt+r.cols().begin());
            union_if_equal(fr_idx,r.cols().begin(),fr_idx-1,r.cols().begin());
        }
        for (size_t fc_idx=r.cols().begin()+1; fc_idx<r.cols().end();fc_idx++) {
            m_dset->make_set(r.rows().begin()*jcnt+fc_idx);
            union_if_equal(r.rows().begin(),fc_idx,r.rows().begin(),fc_idx-1);
        }

        // The invariant for this loop is that all lower rows are already
        // sets and all previous column entries in the same row are sets.
        // The cursor for this loop is at the new (i,j), which is not yet
        // among the sets. It compares with previous (i,j) to see if it
        // should be unioned with them.
        for (size_t i=r.rows().begin()+1; i<r.rows().end(); i++) {
            for (size_t j=r.cols().begin()+1; j<r.cols().end(); j++) {
                m_dset->make_set(i*jcnt+j);
                
                // Only need to add this gridpoint to one set.
                // Check the row first to get most locality.
                if (!union_if_equal(i,j,i,j-1)) {
                    union_if_equal(i,j,i-1,j);
                }
            }
        }
        this->add_edges(r.rows(),r.cols());
    }

    bool union_if_equal(size_t ai, size_t aj,size_t bi, size_t bj) {
        bool added=false;
        if (m_raster(ai,aj)==m_raster(bi,bj)) {
            m_dset->union_set(ai*m_row_cnt+aj,bi*m_row_cnt+bj);
            added=true;
        }
        return added;
    }


};



/*! TBB version 0 of clustering algorithm.
 */
std::shared_ptr<cluster_t> clusters_tbb0(const landscape_t& raster)
{
    // This needs to be a reduce, so we can combine dsets at each
    // reduce step.
    auto cs=ConnectSets(raster);
    parallel_reduce( blocked_range2d<landscape_t::size_type>(
                           0,raster.size1(),32,
                           0,raster.size2(),32),
                  cs
                  );

    return gather_clusters(*(cs.m_parent_pmap), *(cs.m_dset), raster.size1(),
                           raster.size2());
}



} // namespace
