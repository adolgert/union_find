#ifndef _ASSOCIATIVE_ARRAY_H_
#define _ASSOCIATIVE_ARRAY_H_ 1

#include <boost/ublas/matrix.hpp>

namespace raster_stats {

/*! AssociativeArray is a map concept that uses a matrix underneath.
 */
  template<typename _Key, typename _Tp, typename _Compare = std::less<_Key>,
    typename _Allocator = std::allocator<std::pair<const _Key,_Tp> > >
class AssociativeArray {
  typedef boost::ublas::matrix<_Tp> _Base;
  typedef typename _Base::const_iterator _Base_const_iterator;
  typedef typename _Base::iterator _Base_iterator;

  _Base _m;

  public:
  typedef _Key key_type;
  typedef _Tp mapped_type;
  typedef std::pair<const _Key, _Tp> value_type;
  typedef _Compare key_compare;
  typedef _Allocator allocator_type;
  //typedef typename _Base::reference reference;
  //typedef typename _Base::const_reference const_reference;

      typedef typename _Base::size_type             size_type;
      typedef typename _Base::difference_type       difference_type;
      typedef typename _Base::pointer               pointer;
      typedef typename _Base::const_pointer         const_pointer;
      typedef std::reverse_iterator<iterator>       reverse_iterator;
      typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
  

};


}

#endif // _ASSOCIATIVE_ARRAY_H_

