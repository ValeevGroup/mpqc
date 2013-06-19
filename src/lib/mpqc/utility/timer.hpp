#ifndef MPQC_UTIL_TIMER_HPP
#define MPQC_UTIL_TIMER_HPP

#include <iostream>
#include <boost/date_time/posix_time/posix_time_types.hpp>

namespace mpqc {

    struct timer {
	typedef boost::posix_time::time_duration value_type;
	boost::posix_time::ptime start_;
	typedef boost::posix_time::microsec_clock microsec_clock;
	timer() {
	    this->reset();
	}
	void reset() {
	    start_ = microsec_clock::universal_time();
	}
	boost::posix_time::time_duration duration() const {
	    return microsec_clock::universal_time() - start_;
	}
	operator value_type() const {
	    return duration();
	}
	operator double() const {
	    return duration().total_microseconds()/1e6;
	}
	long total_seconds() const {
	    return duration().total_seconds();
	}
    };

    inline std::ostream& operator << (std::ostream &os, const timer &t) {
	//os << boost::posix_time::to_simple_string(t.duration());
	os << double(t);
	return os;// << boost::posix_time::to_simple_string(t.duration());
    }

}

#endif // MPQC_UTIL_TIMER_HPP
