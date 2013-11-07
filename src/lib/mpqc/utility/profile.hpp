#ifndef MPQC_PROFILE_HPP
#define MPQC_PROFILE_HPP

// #include "boost/utility/timer.hpp"

#include "boost/date_time/posix_time/posix_time_types.hpp"

#include <string>
#include <sstream>
#include <map>
#include <typeinfo>

#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/typeof/typeof.hpp>

/**
   @defgroup UtilityProfiler Profile
   @ingroup CoreUtility

   Profile module is a set of macros to measure time spent in a particular scope on the fly.
   To use it, first compile the code with MPQC_PROFILE_ENABLE defined.
   Register the thread before you start profiling by putting somewhere in the beginning
   @code
   MPQC_PROFILE_REGISTER_THREAD;
   @endcode
   Insert MPQC_PROFILE_LINE in the scope of the code you want profiled, eg
   @code
   {
       MPQC_PROFILE_LINE;
       do_lots_of_work_function();
   }
   @endcode
   Dump (and reset/clear) profiling data with MPQC_PROFILE_DUMP(cout), eg:
   @code
   MPQC_PROFILE_DUMP(std::cout);
   @endcode
   The output will be grouped by thread id and will contain file:line keys followed by total
   time spent in the scope and the number of time the scope was entered.
   
   @warning Dump/reset/clear profiler only after the events have finished, otherwise
   there will be difficult to track segfaults.  Better yet, one needs to fix profiler
   to generate error in those cases.

*/

#ifndef DOXYGEN

#define MPQC_PROFILE__(tag) BOOST_PP_CAT(boost_profile__,  tag)

#ifdef MPQC_PROFILE_ENABLE
#define MPQC_PROFILE_FUNCTION(...)					\
    mpqc::profiler::event MPQC_PROFILE__(__LINE__)		\
    (mpqc::global_profiler()					\
     [mpqc::detail::profiler_event				\
      (__PRETTY_FUNCTION__)(__VA_ARGS__)])

#define MPQC_PROFILE_TYPE(T, ...)					\
    mpqc::profiler::event MPQC_PROFILE__(__LINE__)		\
    (mpqc::global_profiler()					\
     [mpqc::detail::profiler_event				\
      (typeid(T).name())(__VA_ARGS__)])
	
#define MPQC_PROFILE(...)					\
    mpqc::profiler::event MPQC_PROFILE__(__LINE__)	\
    (mpqc::global_profiler()				\
     [mpqc::detail::profiler_event(__VA_ARGS__)])

#define MPQC_PROFILE_LINE							\
    mpqc::profiler::event MPQC_PROFILE__(__LINE__)			\
    (mpqc::global_profiler()						\
     [mpqc::detail::profiler_event(std::string(__FILE__) + ":" +	\
					     (BOOST_PP_STRINGIZE(__LINE__)))])

#define MPQC_PROFILE_REGISTER_THREAD			\
    mpqc::global_profiler().register_thread()

#else
#define MPQC_PROFILE_FUNCTION(...)
#define MPQC_PROFILE_TYPE(T, ...)
#define MPQC_PROFILE_LINE
#define MPQC_PROFILE(...)
#define MPQC_PROFILE_REGISTER_THREAD
#endif

#define MPQC_PROFILE_DUMP(stream)			\
    mpqc::global_profiler().dump((stream))

#define MPQC_PROFILE_RESET			\
    mpqc::global_profiler().clear()

namespace mpqc {
namespace detail {

    struct profiler_event {
	profiler_event(const std::string &key) : data_(key) {}
	operator const std::string&() const { return data_; }
	profiler_event& operator()(const std::string &key) {
	    data_ += (":" + key);
	    return *this;
	}
	profiler_event& operator()() { return *this; }
    private:
	std::string data_;
    };

    template<class>
    struct profiler {

	typedef std::string event_key;

	struct event_data {
	    event_data(): size_(0), value_(0) {}
	    event_data(const event_data &e)
		: size_(e.size_), value_(e.value_) {} 
	    event_data& operator+=(double t) {
		++size_;
		value_ += t;
		return *this;
	    }
	    event_data& operator++() { return (*this += 1); }
	    template<class O>
	    O& to_stream(O &ostream) const {
		ostream << value_ << "/" << size_;
		return ostream;
	    }
	private:
	    size_t size_;
	    double value_;
	};

	struct event {
	    event(event_data *data)
		: data_(data),
		  start_(microsec_clock::universal_time()) {}
	    ~event() {
		boost::posix_time::time_duration t = 
		    (microsec_clock::universal_time() - start_);
		// std::cout << timer_ << std::endl;
		if (data_) {
		    *data_ += t.total_microseconds()/1.0e6;
		}
	    }
	private:
	    typedef boost::posix_time::microsec_clock microsec_clock;
	    event_data *data_;
	    boost::posix_time::ptime start_;
	    // utility::timer timer_;
	};

	void register_thread() {
	    boost::lock_guard<boost::mutex> lock(mutex_);
	    //std::cout << "register thread " << thread();
	    events_[thread()];
	}

	event_data* operator[](const event_key &key) {
	    //std::cout << key << std::endl;
	    std::string thread = this->thread();
	    boost::lock_guard<boost::mutex> lock(mutex_);
	    if (events_.find(thread) == events_.end()) return NULL;
	    //std::cout << thread << ": " << key << std::endl;
	    return &events_[thread][key];
	}

	void clear() {
	    boost::lock_guard<boost::mutex> lock(mutex_);
	    BOOST_AUTO(it, events_.begin());
	    while (it != events_.end()) {
		it++->second.clear();
	    }
	}

	template<class S>
	S& to_stream(S &ostream) const {
	    boost::lock_guard<boost::mutex> lock(mutex_);
	    typedef std::map<event_key, event_data> events_t;
	    typename std::map<std::string, events_t>::const_iterator
		e = events_.begin();
	    //std::cout << events_.size() << std::endl;
	    while (e != events_.end()) {
		ostream << "Thread " << e->first << ":" << std::endl;
		typename events_t::const_iterator it = e->second.begin();
		while (it != e->second.end()) {
		    ostream << "    " << it->first << ": ";
		    it->second.to_stream(ostream);
		    ostream << std::endl;
		    ++it;
		}
		++e;
	    }
	    return ostream;
	}

	template<class S>
	void dump(S &stream) {
	    this->to_stream(stream);
	    clear();
	}

    private:
	mutable boost::mutex mutex_;
	std::map<
            std::string, std::map<event_key, event_data>
	> events_;

	static std::string thread() {
	    std::stringstream ss;
	    ss << boost::this_thread::get_id();
	    return ss.str();
	}
    };

    // template<typename T>
    // profiler<T> profiler<T>::global;

    template<typename T>
    inline std::ostream& operator<<(std::ostream &ostream, const profiler<T> &p) {
	return p.to_stream(ostream);
    }
    
} // namespace detail
} // namespace mpqc


namespace mpqc {

    typedef detail::profiler<void> profiler;
    inline profiler& global_profiler() {
	static profiler p;
	return p;
    }

}

#endif // DOXYGEN

#endif // MPQC_PROFILE_HPP
