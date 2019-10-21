/*! raster_wrap.coo
 *  This file is the main wrapper around the code to find clusters in landscapes.
 */
#include <functional>
#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/numeric.hpp>
#include <boost/multi_array.hpp>
#ifdef USE_TCMALLOC
#include "google/heap-profiler.h"
#endif
#include "numpy/arrayobject.h"
#include "timing_harness.hpp"
#include "timing.hpp"
#include "raster.hpp"
#include "cluster.hpp"

using namespace raster_stats;
using namespace std;
using namespace boost::python;


boost::array<landscape_t::size_type, 2> numpy_dimensions(PyObject* array_py) {
	Py_intptr_t dimension_cnt = PyArray_NDIM(array_py);
	if (dimension_cnt<1) {
		throw runtime_error("Could not find number of dimensions in array.");
	}
	
	boost::array<landscape_t::size_type, 2> dimensions;
	Py_intptr_t* py_dimensions = PyArray_DIMS(array_py);
	if (0==py_dimensions) {
		throw runtime_error("Could not read each dimension of array.");
	}
	for (int cpy_idx=0; cpy_idx<dimension_cnt; cpy_idx++) {
		dimensions[cpy_idx] = py_dimensions[cpy_idx];
	}
	
	return dimensions;
}


/*! Traits type for extracting from Numpy arrays.
 *  kind: b=bool, i=signed int, u=unsigned, f=float, c=complex float, 'S'=8-bit string
 *        'U' = 32-bit unicode string, 'V'=arbitrary
 *        constexpr string type_num[] = {
 *        "bool","int8","uint8","int16","uint16","int32","uint32","int64",
 *        "uint64","float32","float64","float128","complex64","complex128",
 *        "complex256","object","string","unicode","void" };
 */
template<typename T>
struct numpy_type {
	static const char kind;
	static const char type_num;
};


template<>
struct numpy_type<unsigned char> {
	static const char kind     ='u';
	static const char type_num =2;
};




template<typename T>
const boost::multi_array<T,2> numpy_array_extract(PyObject* array)
{
	PyArray_Descr* description = PyArray_DESCR(array);
	if (description->kind != numpy_type<T>::kind) {
		throw runtime_error("Tried to extract from array of different kind.");
	}
	if (description->type_num != numpy_type<T>::type_num) {
		throw runtime_error("Tried to extract from array of different type.");
	}
	if (description->elsize != sizeof(T)) {
		throw runtime_error("Tried to extract from array with elements of different size.");
	}

	const arr_type* val = static_cast<const T*>(PyArray_DATA(array));
	typedef boost::const_multi_array_ref<T,2> clandscape_t;
	clandscape_t image(val,numpy_dimensions(array));
	return image;
}


/*! This class is a Python iterator created to present a list<list<size_t>> as a list of python lists.
 */
struct ClusterIter {
	cluster_t::iterator m_begin, m_end;
	ClusterIter(cluster_t::iterator begin, cluster_t::iterator end)
		: m_begin(begin), m_end(end) {}

	boost::python::list next() {
		boost::python::list retlist;
		if (m_begin==m_end) {
			PyErr_SetString(PyExc_StopIteration, "No more data.");
			boost::python::throw_error_already_set();
		} else {
			const cluster_t::value_type& element_list = *m_begin;
			for (auto elem_idx=element_list.begin(); elem_idx!=element_list.end(); elem_idx++) {
				retlist.append(*elem_idx);
			}			
		}
		++m_begin;
		return retlist;
	}	
};


/*! Python needs some way to read the clusters that are returned.
 *  This class is a simple wrapper around that. This class is not
 *  for daily use but for unit tests in Python.
 */
struct ClusterWrap {
	cluster_t m_clusters;
	ClusterWrap() = default; //! Wrapping by Boost.Python requires a default constructor.
	ClusterWrap(cluster_t& clusters) : m_clusters(clusters) {}
	ClusterWrap(const ClusterWrap& other) : m_clusters(other.m_clusters) {}
	size_t size() { return m_clusters.size(); }
	
	//sublist_iterator begin() { return sublist_iterator(m_clusters.begin()); }
	//sublist_iterator end() { return sublist_iterator(m_clusters.end()); }

	ClusterIter get_iterator() { return ClusterIter(m_clusters.begin(),m_clusters.end()); }
	
	boost::python::list at(size_t idx) {
		cluster_t::const_iterator iter = m_clusters.begin();
		// This is an N^2 behavior. Method is just for testing.
		while (iter!=m_clusters.end() && idx!=0) {
			++iter;
			--idx;
		}
		boost::python::list retlist;
		for (auto elem_idx=iter->begin(); elem_idx!=iter->end(); elem_idx++) {
			retlist.append(*elem_idx);
		}
		return retlist;
	}
};



struct ClusterPointerWrap {
  shared_ptr<cluster_t> m_clusters;
	ClusterPointerWrap() = default; //! Wrapping by Boost.Python requires a default constructor.
  ClusterPointerWrap(shared_ptr<cluster_t> clusters) : m_clusters(clusters) {}
	ClusterPointerWrap(const ClusterPointerWrap& other) : m_clusters(other.m_clusters) {}
	size_t size() { return m_clusters->size(); }
	
	ClusterIter get_iterator() { return ClusterIter(m_clusters->begin(),m_clusters->end()); }
	
	boost::python::list at(size_t idx) {
		cluster_t::const_iterator iter = m_clusters->begin();
		// This is an N^2 behavior. Method is just for testing.
		while (iter!=m_clusters->end() && idx!=0) {
			++iter;
			--idx;
		}
		boost::python::list retlist;
		for (auto elem_idx=iter->begin(); elem_idx!=iter->end(); elem_idx++) {
			retlist.append(*elem_idx);
		}
		return retlist;
	}
};




struct ClusterWrapPair {
	cluster_loc_t m_clusters;
	ClusterWrapPair() = default;
	ClusterWrapPair(cluster_loc_t& clusters) : m_clusters(clusters) {}
	ClusterWrapPair(const ClusterWrapPair& other) : m_clusters(other.m_clusters) {}
	size_t size() { return m_clusters.size(); }
	
	boost::python::list at(size_t idx) {
		cluster_loc_t::const_iterator iter = m_clusters.begin();
		// This is an N^2 behavior. Method is just for testing.
		while (iter!=m_clusters.end() && idx!=0) {
			++iter;
			--idx;
		}
		boost::python::list retlist;
		const cluster_loc_t::mapped_type& element_list = iter->second;
		for (auto elem_idx=element_list.begin(); elem_idx!=element_list.end(); elem_idx++) {
			retlist.append(boost::python::make_tuple(elem_idx->first, elem_idx->second));
		}
		return retlist;
	}
};



ClusterWrapPair* find_clusters_pair_wrap(object raster_object) {

	const landscape_t raster = numpy_array_extract<arr_type>(raster_object.ptr());
	cluster_loc_t clusters = find_clusters_pair(raster);
	ClusterWrapPair* clusters_wrap = new ClusterWrapPair(clusters);
	
	return clusters_wrap;
}


ClusterWrap* find_clusters_wrap(object raster_object) {

	const landscape_t raster = numpy_array_extract<arr_type>(raster_object.ptr());
#ifdef USE_TCMALLOC
	HeapProfilerStart("find_clusters");
#endif
	cluster_t clusters = find_clusters(raster);
#ifdef USE_TCMALLOC
	HeapProfilerDump("profile");
	HeapProfilerStop();
#endif
	ClusterWrap* clusters_wrap = new ClusterWrap(clusters);

	return clusters_wrap;
}


ClusterWrap* find_clusters_twopass_wrap(object raster_object) {

	const landscape_t raster = numpy_array_extract<arr_type>(raster_object.ptr());
#ifdef USE_TCMALLOC
	HeapProfilerStart("find_clusters_twopass");
#endif
	cluster_t clusters = find_clusters_twopass(raster);
#ifdef USE_TCMALLOC
	HeapProfilerDump("profile");
	HeapProfilerStop();
#endif
	ClusterWrap* clusters_wrap = new ClusterWrap(clusters);

	return clusters_wrap;
}

long long find_clusters_pair_time_wrap(size_t n, object raster_object) {
	const landscape_t raster = numpy_array_extract<arr_type>(raster_object.ptr());
	return timeit([&raster](){ find_clusters_pair(raster); }, n).count();
}


long long find_clusters_time_wrap(size_t n, object raster_object) {
	const landscape_t raster = numpy_array_extract<arr_type>(raster_object.ptr());
	return timeit([&raster](){ find_clusters(raster); }, n).count();
}


long long find_clusters_twopass_time_wrap(size_t n, object raster_object) {
	const landscape_t raster = numpy_array_extract<arr_type>(raster_object.ptr());
	return timeit([&raster](){ find_clusters_twopass(raster); }, n).count();
}

ClusterPointerWrap* find_clusters_pointer_wrap(object raster_object) {

	const landscape_t raster = numpy_array_extract<arr_type>(raster_object.ptr());
#ifdef USE_TCMALLOC
	HeapProfilerStart("profile");
#endif
	shared_ptr<cluster_t> clusters = find_clusters_pointer(raster);
#ifdef USE_TCMALLOC
	HeapProfilerDump("profile");
	HeapProfilerStop();
#endif
	ClusterPointerWrap* clusters_wrap = new ClusterPointerWrap(clusters);

	return clusters_wrap;
}

long long find_clusters_pointer_time_wrap(size_t n, object raster_object) {
	const landscape_t raster = numpy_array_extract<arr_type>(raster_object.ptr());
	return timeit([&raster](){ find_clusters_pointer(raster); }, n).count();
}


ClusterWrap* find_clusters_remap_wrap(object raster_object) {

	const landscape_t raster = numpy_array_extract<arr_type>(raster_object.ptr());
#ifdef USE_TCMALLOC
	HeapProfilerStart("profile");
#endif
	cluster_t clusters = find_clusters_remap(raster);
#ifdef USE_TCMALLOC
	HeapProfilerDump("profile");
	HeapProfilerStop();
#endif
	ClusterWrap* clusters_wrap = new ClusterWrap(clusters);

	return clusters_wrap;
}

long long find_clusters_remap_time_wrap(size_t n, object raster_object) {
	const landscape_t raster = numpy_array_extract<arr_type>(raster_object.ptr());
	return timeit([&raster](){ find_clusters_remap(raster); }, n).count();
}

boost::python::list get_list(void) {
	boost::python::list retlist;
	retlist.append(boost::python::make_tuple(3,4));
	retlist.append(boost::python::make_tuple(7,6));
	retlist.append(24);
	return retlist;
}


boost::python::list unique_values_wrap(object raster_object) {
	const landscape_t raster = numpy_array_extract<arr_type>(raster_object.ptr());
	set<arr_type> vals = unique_values_direct(raster);
	boost::python::list retlist;
	for (auto idx=vals.begin(); idx!=vals.end(); idx++) {
		retlist.append(*idx);
	}
	return retlist;
}


inline object pass_through(object const& o) { return o; }


void dump_heap()
{
#ifdef USE_TCMALLOC
  HeapProfilerDump("find_clusters_remap");
#endif
}


// Expose classes and methods to Python
BOOST_PYTHON_MODULE(raster_stats) {
	boost::python::numeric::array::set_module_and_type("numpy","ndarray");
	class_<ClusterWrap>("Clusters")
			.def("__len__", &ClusterWrap::size)
			.def("__iter__", &ClusterWrap::get_iterator)
			.def("at", &ClusterWrap::at)
			;

	class_<ClusterPointerWrap>("Clusters")
			.def("__len__", &ClusterPointerWrap::size)
			.def("__iter__", &ClusterPointerWrap::get_iterator)
			.def("at", &ClusterPointerWrap::at)
			;

	class_<ClusterIter>("ClusterIter", no_init)
	        .def("next", &ClusterIter::next)
	        .def("__iter__", pass_through)
	      ;

	def( "dump_heap", dump_heap ) ;
	def( "unique_values", unique_values_wrap) ;

	def( "find_clusters_fourpass", find_clusters_wrap, return_value_policy<manage_new_object>() ) ;
	def( "find_clusters_fourpass_time", find_clusters_time_wrap ) ;
	def( "find_clusters", find_clusters_twopass_wrap, return_value_policy<manage_new_object>() ) ;
	def( "find_clusters_time", find_clusters_twopass_time_wrap ) ;
	def( "find_clusters_pointer", find_clusters_pointer_wrap, return_value_policy<manage_new_object>() ) ;
	def( "find_clusters_pointer_time", find_clusters_pointer_time_wrap ) ;
	def( "find_clusters_remap", find_clusters_remap_wrap, return_value_policy<manage_new_object>() ) ;
	def( "find_clusters_remap_time", find_clusters_remap_time_wrap ) ;

	def("get_list", get_list) ;
}
