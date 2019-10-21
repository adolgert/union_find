#ifndef _TIMING_HARNESS_H_
#define _TIMING_HARNESS_H_ 1

#include <string>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace raster_stats {


/*! This runs functions by name. The point is to load C++ functions
 *  into this harness and then paw through them in Python, running
 *  what we want.
 */

    class timing_harness {
    public:
        virtual long long time(size_t n) const = 0;
        virtual const std::string& name() const = 0;
    };



    template<class SUBJECT>
    class timing_harness_test : public timing_harness {
        SUBJECT subject_;
        std::string name_;

    public:
        timing_harness_test(SUBJECT subject, std::string name) :
            subject_(subject), name_(name) {}
        
        virtual long long time(size_t n) const
        {
            return timeit(subject_,n).count();
        }
	
        virtual const std::string& name() const { return name_; };
    };


    template<class SUBJECT>
    boost::shared_ptr<timing_harness_test<SUBJECT>>
    make_timing(SUBJECT subject, std::string name) {
        return boost::make_shared<timing_harness_test<SUBJECT>>(subject,name);
    }

}

#endif
