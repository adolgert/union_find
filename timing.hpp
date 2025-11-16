#ifndef _TIMING_H_
#define _TIMING_H_ 1
/*! timing.cpp
 *  This file does timing on other functions.
 */
#include <time.h>
#include <sys/resource.h>
#include <memory>
#include <boost/chrono.hpp>



/*! Do standard clock timing.
 */
template<typename F>
double time_clock(F f, size_t run_cnt) {
	clock_t start = clock();
	for (size_t i=0; i<run_cnt; i++) {
		f();
	}
	clock_t finish = clock();
	double secs_per_iteration = ((finish-start)/((double) CLOCKS_PER_SEC)*run_cnt);
	return secs_per_iteration;
}


typedef boost::chrono::nanoseconds ns_t;

template<typename F>
ns_t timeit(F f, size_t run_cnt) {
	auto start = boost::chrono::high_resolution_clock::now();	
	
	for (size_t i=0; i<run_cnt; i++) {
		f();
	}
	
	ns_t sec = boost::chrono::high_resolution_clock::now() - start;
	
	return sec/run_cnt;
}



/*! Get rusage.
 */
template<typename F>
std::shared_ptr<struct rusage> usage(F f, size_t run_cnt) {
	struct rusage begin_usage;
	int res_start = getrusage(RUSAGE_SELF,&begin_usage);
	if (0 != res_start) {
		if (EFAULT == errno) {
			throw std::runtime_error("rusage reuturned efault, so the address is invalid.");
		} else if (EINVAL == errno) {
			throw std::runtime_error("rusage returne EINVAL, so the who parameter is incorrect.");
		} else {
			throw std::runtime_error("rusage returned a completely unknown error");
		}
	}
	
	for (size_t i=0; i<run_cnt; i++) {
		f();
	}
	
	std::shared_ptr<struct rusage> end_usage(new struct rusage);
	int res_finish = getrusage(RUSAGE_SELF,end_usage.get());
	if (0 != res_finish) {
		if (EFAULT == errno) {
			throw std::runtime_error("rusage reuturned efault, so the address is invalid.");
		} else if (EINVAL == errno) {
			throw std::runtime_error("rusage returne EINVAL, so the who parameter is incorrect.");
		} else {
			throw std::runtime_error("rusage returned a completely unknown error");
		}
	}
	return end_usage;
}


#endif // _TIMING_H_
