#ifndef CC_UTILITY_HPP
#define CC_UTILITY_HPP

#include <stdexcept>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/multi_array/base.hpp>

namespace cchem {
namespace cc {

    template<typename T>
    struct Buffer {
	struct vector {
	    operator T*() { return data(); }
	    T* data() { return &data_[0]; }
	    void resize(size_t size) { data_.reserve(size); }
	    size_t size() const { return data_.capacity(); }
	    size_t bytes() const { return size()*sizeof(T); }
	private:
	    std::vector<T> data_;
	};
	vector& operator[](int i) {
	    //assert(data_.count(i) > 0);
	    return data_[i];
	}
	void clear() { data_.clear(); }
    private:
	std::map<int, vector> data_;
    };

    // template<typename T>
    // struct Buffer {
    // 	std::vector<T> data_;
    // 	operator T*() { return data(); }
    // 	T* data() { return &data_[0]; }
    // 	void reserve(size_t size) { data_.reserve(size); }
    // };

    template<class C>
    struct Symbol {

	struct range : boost::numeric::ublas::range {
	    typedef boost::numeric::ublas::range base;
	    range(const size_t start, const size_t finish)
		: base(start, finish) {}
	    template<typename R>
	    explicit range(const R &r, typename boost::enable_if<
			   boost::is_fundamental<R> >::type* = 0)
		: base(r, r+1) {}
	    template<typename R>
	    explicit range(const R &r, typename boost::disable_if<
			   boost::is_fundamental<R> >::type* = 0)
		: base(r.start(), r.start()+r.size()) {}
	    size_t finish() const { return (this->start()+this->size()); }
	    typedef boost::multi_array_types::extent_range extent_range;
	    operator extent_range() const {
		return extent_range(this->start(), this->finish());
	    }
	};

	struct no_such_key : std::runtime_error {
	    no_such_key(const std::string &s) : std::runtime_error(s) {}
	};

	typedef C& reference;
	typedef const C& const_reference;
	typedef typename C::pointer pointer;
	reference operator[](const std::string &key) {
	    if (!contains(key)) throw no_such_key(key);
	    return *data_.find(key)->second;
	}
	// const_reference operator[](const std::string &key) const {
	//     
	//     return data_.find(key)->second;
	// }

	bool contains(const std::string &key) const {
	    return (data_.count(key) > 0);
	}

	void erase(std::string key) {
	    data_.erase(key);
	}

	void clear() {
	    data_.clear();
	}

	std::vector<std::string> keys() const {
	    std::vector<std::string> keys;
	    BOOST_AUTO(it, data_.begin());
	    while (it != data_.end()) {
		keys.push_back(it->first);
		++it;
	    }
	    return keys;
	}

	template<class E>
	reference set(std::string key, pointer data, const E &extents) {
	    erase(key);
	    data_.insert(key, new C(data, extents));
	    return this->operator[](key);
	}

	reference set(std::string key, pointer data) {
	    return this->set(key, data, boost::extents[0][0][0][0]);
	}

	template<class R1, class R2, class R3, class R4>
	reference set(std::string key, pointer data, 
		      const R1 &r1, const R2 &r2, const R3 &r3, const R4 &r4) {
	    size_t r[] = { range(r1).size(), range(r2).size(),
			   range(r3).size(), range(r4).size() };
	    return this->set(key, data, boost::extents[r[3]][r[2]][r[1]][r[0]]);
	}

	template<class E>
	reference reset(std::string key, const E &extents) {
	    pointer data = this->operator[](key).data();
	    return this->set(key, data, extents);
	}

	template<class R1, class R2, class R3, class R4>
	reference reset(std::string key,
			const R1 &r1, const R2 &r2, const R3 &r3, const R4 &r4) {
	    pointer data = this->operator[](key).data();
	    return this->set(key, data, r1, r2, r3, r4);
	}

	template<class R1, class R2, class R3, class R4, class A>
	reference load(const std::string &key, pointer data, 
		       const R1 &r1, const R2 &r2, const R3 &r3, const R4 &r4,
		       const A &a) {
	    reference ref = set(key, data, r1, r2, r3, r4);
	    range r[] = { range(r1), range(r2), range(r3), range(r4) };
	    size_t start[4], stop[4];
	    //std::cout << data << " ";
	    for (size_t i = 0; i < 4; ++i) {
		start[i] = r[i].start();
		stop[i] = start[i] + r[i].size();	    
		//std::cout <<  stop[i] - start[i] << " ";
	    }
	    //std::cout << std::endl;
	    a.get(ref.data(), start, stop);
	    return ref;
	}

	template<class R1, class R2, class R3, class R4, class A>
	reference load(const std::string &key, const A &a, 
		       const R1 &r1, const R2 &r2, const R3 &r3, const R4 &r4) {
	    BOOST_AUTO(&s, this->operator[](key));
	    return this->load(key, s.data(), r1, r2, r3, r4, a);
	}

	template<class R1, class R2, class R3, class R4, class A>
	void store(const std::string &key,
		   const R1 &r1, const R2 &r2, const R3 &r3, const R4 &r4,
		   A &a) {
	    range r[] = { range(r1), range(r2), range(r3), range(r4) };
	    size_t start[4], stop[4];
	    for (size_t i = 0; i < 4; ++i) {
		start[i] = r[i].start();
		stop[i] = start[i] + r[i].size();
	    }
	    a.put(this->operator[](key).data(), start, stop);	    
	}
	void rename(std::string name, std::string rename) {
	    this->operator[](name);
	    BOOST_AUTO(ptr, data_.release(data_.find(name)));
	    data_.insert(rename, ptr.release());
	}
	void swap(std::string a, std::string b) {
	    this->operator[](a);
	    this->operator[](b);
	    BOOST_AUTO(pa, data_.release(data_.find(a)));
	    BOOST_AUTO(pb, data_.release(data_.find(b)));
	    data_.insert(a, pb.release());
	    data_.insert(b, pa.release());
	}
    private:
	boost::ptr_map<std::string, C> data_;
    };

} // namespace cc
} // namespace cchem

#endif /* CC_UTILITY_HPP */
