#ifndef _UNIQUE_VALUES_H_
#define _UNIQUE_VALUES_H_ 1

#include <functional>
#include <boost/property_map/property_map.hpp>
#include <boost/concept/assert.hpp>


namespace raster_stats {

  /*! This iterates a unary operator over an array.
   *  Because we may use Multiarrays or ublas arrays,
   *  different specializations of this function should
   *  iterate over both.
   */
  template<class Array, class Function>
    void for_array(const Array& arr, Function f)
  {
    BOOST_CONCEPT_ASSERT((boost::UnaryFunction<Function,void,
                          typename Array::value_type>) );

    for (auto iter1=arr.begin1(); iter1!=arr.end1(); ++iter1) {
      for (auto iter2=iter1.begin(); iter2!=iter1.end(); ++iter2) {
        f(*iter2);
      }
    }
  }

  template<class LandscapeType,class WritableSet>
    struct insert_unique :
  public std::unary_function<typename LandscapeType::value_type,void>
  {
    WritableSet& uniques_;
  insert_unique(WritableSet& uniques) : uniques_(uniques) {}
    void operator()(typename LandscapeType::value_type landuse) {
      uniques_.insert(landuse);
    }
  };


  template<class LandscapeType, class SetType>
    void unique_values(const LandscapeType& raster,
                                           SetType& uniques)
  {
    auto inserter = insert_unique<LandscapeType,SetType>(uniques);
    for_array(raster,inserter);
  }

}

#endif // _UNIQUE_VALUES_H_

