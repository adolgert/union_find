#ifndef _BLOCKED_ARRAY_H_
#define _BLOCKED_ARRAY_H_ 1

#include <vector>
#include <array>
#include <stdexcept>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/type_traits/has_trivial_constructor.hpp>
#include <boost/iterator/iterator_facade.hpp>


namespace raster_stats {


  /*! This walks through blocks as though they were in order.
   *  A is the array this iterates over.
   *  http://www.boost.org/doc/libs/1_48_0/libs/iterator/doc/iterator_facade.html
   */
  template<class T,size_t B>
  class blocked_iterator :
  public boost::iterator_facade<blocked_iterator,
    typename A::value_type,
    boost::random_access_traversal_tag>
    {
      T* data_;
      size_t bi_; //! i of block
      size_t bj_; //! j of block
      size_t ri_; //! running i within the block
      size_t rj_; //! running j withing the block
      size_t bw_; //! blocks width
      size_t rw_; //! width within the last block
    public:
      blocked_iterator()
        : bi_(0), bj_(0), ri_(0), rj_(0), bw_(0), rw_(0)
        {
          data_ = 0;
        }
      explicit blocked_iterator(T* array,size_t width,size_t i, size_t j)
        : 
      {
        data_ = array;
        bi_=i/(B*B);
        bj_=j/(B*B);
        bw_=(width-1)/B + 1;
        rw_=width/bw_;
        ri_=i%(B*B);
        rj_=j%(B*B);
      }
    private:
      friend class boost::iterator_core_access;

      void increment() {
        ++rj_;
        // The row can end mid-block for the last column.
        if (rj_==B || (rj_==rw_ && bj_==bw_-1)) {
          rj_=0;
          ++bj_;
          if (bj_==bw_) {
            bj_=0;
            ++ri_;
            if (ri_==B) {
              ri_=0;
              ++bi_;
            }
          }
        }
      }

      bool equal(blocked_iterator const& b) const
      {
        return ri_==b.ri_ && rj_==b.rj_ && bi_==b.bi_ && bj_==b.bj_ &&
          data_==b.data_;
      }

      T& dereference() const {
        return data_[((bi_*bw_+bj_)*B+ri_)*B+rj_];
      }

      void decrement() {
        if (rj_>0) {
          --rj_;
        } else {
          if (bj_>0) {
            --bj_;
            rj_=B-1;
          } else {
            rj_=rw_-1;

            if (ri_>0) {
              --ri_;
            } else {
              --bi_;
              ri_=B;
            }
          }
        }
      }

      void advance(size_t n)
      {
        // This would be the index into a regular 1-d storage.
        size_t logical_idx=((bi_*(bw_-rw_)+bj_)*B+ri_)*B+rj_;
        logical_idx+=n;
        bi_=logical_idx/(B*B);
        bj_=logical_idx/(B*B);
        bw_=(width-1)/B + 1;
        rw_=width/bw_;
        ri_=logical_idx%(B*B);
        rj_=logical_idx%(B*B);
      }
      size_t distance_to(size_t j)
      {
        size_t logical_idx=((bi_*(bw_-rw_)+bj_)*B+ri_)*B+rj_;
        return j-logical_idx;
      }
    };


/*! This implements a storage array using blocks for boost::ublas::matrix.
 *  This stores all data in NxN subarrays which are arranged in Morton
 *  z-order. It can be accessed like a regular array.
 *
 *  T is the element type.
 *  B is the length of a side of a square block.
 *  N is the number of a side of value_type in the square array.
 *
 *  The ublas::matrix needs a linear storage type, but this array needs to
 *  know matrix dimensions. The solution, for now, is for to specify 
 *  dimensions in a template specification.
 *
 *  Finding elements within the array will be computationally expensive.
 *  This is a small experiment to determine whether improvement in locality
 *  from blocked storage is more important or whether the cost of computing
 *  indices overwhelms this.
 *
 *  Implementation follows blocked_array from storage.hpp in Boost.
 */

  template<class T,size_t B, size_t N,size_t M,class ALLOC>
    class blocked_array :
  public boost::numeric::ublas::storage_array<blocked_array<T,B,N,M,ALLOC> > {
    typedef blocked_array<T,B,N,M,ALLOC> self_type;

  public:
    typedef ALLOC allocator_type;
    typedef typename ALLOC::size_type size_type;
    typedef typename ALLOC::difference_type difference_type;
    typedef T value_type;
    typedef const T &const_reference;
    typedef T &reference;
    typedef const T *const_pointer;
    typedef T *pointer;
    typedef const_pointer const_iterator;
    typedef pointer iterator;

    typedef B blocksize;

  private:
    typedef std::vector<std::array<value_type,(B*B)>> storage_type;
    static const size_t NBLOCKS=((N-1)/B+1);
    static const size_t MBLOCKS=((M-1)/B+1);

    ALLOC alloc_;
    size_type size_;
    size_type cnt_;
    pointer data_;

  public:
    explicit blocked_array( const ALLOC &a = ALLOC()) :
    alloc_(a), size_(0), cnt_(NBLOCKS*MBLOCKS*B*B) {
    }
    explicit blocked_array(size_type size, const ALLOC& a = ALLOC()) :
    alloc_(a), size_(0), cnt_(NBLOCKS*MBLOCKS*B*B) {
      if (size == 0) {
        data_ = 0;
      } else if (size_ < cnt_) {
        data_ = alloc_.allocate(cnt_);
      } else {
        throw std::runtime_error("Cannot use this array size.");
      }
    }

    blocked_array (const blocked_array &c) :
    storage_array<blocked_array<T,B,N,M,ALLOC> >(),
      alloc_(c.alloc_), size_(c.size_) {
      if (size_) {
        data_ = alloc_.allocate(cnt_);
        std::uninitialized_copy(c.begin(), c.end(), begin());
      } else {
        data_ = 0;
      }
    }

    ~blocked_array() {
      if (size_) {
        size_t cnt = NBLOCKS*MBLOCKS*B*B;
        alloc_.deallocate(data_,cnt);
      }
    }

  private:
    void resize_internal(const size_type size, const value_type init,
                         const bool preserve) {
      if (size==size_) return;
      pointer p_data = data_;
      if (size) {
        data_ = alloc.allocate(size);
      }
      if (size_) {
        alloc_.deallocate(p_data,cnt_);
      }
      if (!size_) {
        data_ = 0;
      }
      size_ = size;
    }
  public:
    void resize(size_type size) {
      resize_internal(size,value_type(), false);
    }
    void resize(size_type size, value_type init) {
      resize_internal(size,init,true);
    }
    size_type max_size() const {
      return ALLOC().max_size();
    }
    bool empty() const {
      return size_ == 0;
    }
    size_type size() const {
      return size_;
    }
    /*! This returns the array element. It is here we decide that
     *  this array stores in blocks but not in Morton order.
     */
    const_reference operator[](size_type idx) const {
      size_t i=idx%row,  j=idx/row;
      size_t bi=i/(B*B), bj=j/(B*B);
      return data_[(bi*NBLOCKS+bj)*(B*B) + (i%B)*B+(j%B)];
    }
    reference operator[](size_type idx) {
      size_t i=idx%row,  j=idx/row;
      size_t bi=i/(B*B), bj=j/(B*B);
      return data_[(bi*NBLOCKS+bj)*(B*B) + (i%B)*B+(j%B)];
    }
    blocked_array& operator=(const blocked_array& a) {
      if (this != &a) {
        resize(a.size_);
        std::copy(a.data_,a.data_+a.cnt_,data_);
      }
      return *this;
    }
    blocked_array& assign_temporary(blocked_array &a) {
      swap (a);
      return *this;
    }

    void swap(blocked_array &a) {
      if (this != &a) {
        std::swap(size_,a.size_);
        std::swap(data_,a.data_);
      }
    }
    friend void swap(blocked_array &a1, blocked_array &a2) {
      a1.swap(a2);
    }
    const_iterator begin() const {
      return data_;
    }


  };

}


#endif // _BLOCKED_ARRAY_H_
