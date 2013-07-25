#ifndef MPQC_FILE_HPP
#define MPQC_FILE_HPP

#include <string>
#include <cassert>
#include <fstream>

#include <stdlib.h>
#include <hdf5.h>

#include <memory>

#include <boost/array.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_set.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/type_traits.hpp>
#include <boost/assert.hpp>

#include <boost/typeof/typeof.hpp>
#include <boost/fusion/include/make_map.hpp>
#include <boost/fusion/include/at_key.hpp>

#include "mpqc/range.hpp"
#include "mpqc/range/operator.hpp"

#include "mpqc/utility/foreach.hpp"
#include "mpqc/utility/timer.hpp"

#include <util/misc/exenv.h>

/**
 * File stores hierarchical data in an
 * <a href="http://www.hdfgroup.org/HDF5/">HDF5</a> backend.
 * Technically, File is an std::set of File::Handle objects,
 * each of which is a handle to a HDF5 (POSIX) file.
 */

namespace mpqc {
  namespace file {

#define MPQC_FILE_VERIFY(expr) BOOST_VERIFY((expr) >= 0)

    template<typename T>
    static hid_t h5_predtype() {
      using boost::fusion::make_map;
      BOOST_AUTO(pt,
                 (make_map<int, double>(H5T_NATIVE_INT, H5T_NATIVE_DOUBLE)));
      typedef typename boost::remove_const<T>::type U;
      return boost::fusion::at_key < U > (pt);
    }

    struct Object {
        Object() :
            id_(0) {
        }
        Object(const Object &o) {
          *this = o;
        }
        template<class F>
        Object(const Object &parent, hid_t id, F close, bool increment) {
          assert(id);
          parent_.reset(new Object(parent));
          update(id, close, increment);
        }
        ~Object() {
          if (id_) {
            assert(close_);
            close_(id_);
          }
        }
        void operator=(const Object &o) {
          Object *parent = o.parent_.get();
          parent_.reset(parent ? new Object(*parent) : NULL);
          update(o.id_, o.close_, true);
        }
        hid_t id() const {
          return id_;
        }
        const Object& parent() const {
          return *parent_;
        }
        hid_t file() const {
          return H5Iget_file_id(id_);
        }
        static std::string filename(hid_t id) {
          std::vector<char> str(H5Fget_name(id, NULL, 0) + 1);
          MPQC_FILE_VERIFY(H5Fget_name(id, &str[0], str.size()));
          return std::string(&str[0]);
        }
        operator bool() const {
          return (this->id_ != 0);
        }
      protected:
        std::auto_ptr<Object> parent_;
        hid_t id_;
        void (*close_)(hid_t);
        template<class F>
        void update(hid_t id, F close, bool increment) {
          if (id && increment)
            H5Iinc_ref(id);
          id_ = id;
          close_ = close;
        }
    };

    struct Attribute: Object {
//         template<typename T>
//         static Attribute create(Object parent, const std::string &name) {
//             hid_t space = H5Screate(H5S_SCALAR);
//             MPQC_FILE_VERIFY(space);
//             hid_t id =
// #if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
//             H5Acreate(parent.id(), name.c_str(), h5_predtype<T>(), space, H5P_DEFAULT);
// #else
//             H5Acreate1(parent.id(), name.c_str(), h5_predtype<T>(), space, H5P_DEFAULT);
// #endif
//             MPQC_FILE_VERIFY(id);
//             MPQC_FILE_VERIFY(H5Sclose(space));
//             return Attribute(Object(parent, id, &Attribute::close));
//         }
//         static void close(hid_t *id) {
//             H5Aclose(*id);
//         }
//         template<typename T>
//         operator T() const {
//             T value;
//             MPQC_FILE_VERIFY(H5Aread(this->id(), h5_predtype<T>(), &value));
//             return value;
//         }
//         template<typename T>
//         void operator=(const T &value) {
//             MPQC_FILE_VERIFY(H5Awrite(this->id(), h5_predtype<T>(), &value));
//         }
//     private:
//         Attribute(Object o) : Object(o) {}
    };

    template<class Container>
    struct static_container {
        static Container data;
    };
    
    template<class Container>
    Container static_container<Container>::data;

  } // namespace file
} // namespace mpqc

namespace mpqc {

  struct File: file::Object {

      typedef file::Object Object;
      typedef file::Attribute Attribute;
      struct Group;

      template<typename T>
      struct Dataset;
      template<typename T>
      struct Dataspace;

      struct Driver {
          struct Core;
          Driver() :
              fapl_(H5P_DEFAULT) {
          }
          hid_t fapl() const {
            return fapl_;
          }
        private:
          hid_t fapl_;
      };

      // struct Driver::Core : Driver {
      //     explicit Core(size_t increment = 64<<20,
      // 		  hbool_t backing_store = false) {
      // 	printf("HDF5 Core increment = %lu\n", increment);
      // 	Driver::fapl_ = H5Pcreate(H5P_FILE_ACCESS);
      // 	H5Pset_fapl_core(Driver::fapl_, increment, backing_store);
      //     }
      // };

      File() {
      }

      /// creates a File
      /// @param[in] name file name
      explicit File(const std::string &name, const Driver &driver = Driver()) {
        initialize(name, driver);
      }

      static void close(hid_t id) {
        // H5Fclose(id);
        // if (!H5Iget_ref(id)) files_::data.erase(id);
      }

      Group group(const std::string &name = "/"); //  {
      //     return Group::create(*this, name);
      // }

    private:

      typedef file::static_container<std::set<hid_t> > files_;

      void initialize(const std::string &name, const Driver &driver) {
        Object o = File::open(name, driver);
        if (!o)
          o = File::create(name, driver);
        Object::operator=(o);
      }

      explicit File(hid_t id, bool increment) :
          Object(Object(), id, &File::close, increment) {
        assert(id > 0);
      }

      static std::string realpath(const std::string &name) {
        return name;
        // char *str = ::realpath(name.c_str(), NULL);
        // std::cout << "realpath: " << str << std::endl;
        // std::string path(str ? str : "");
        // free(str);
        // return path;
      }

      static File open(const std::string &name, const Driver &driver) {
        std::string path = realpath(name);
        if (path.empty())
          return File();
        foreach (auto id, files_::data) {
          if (path == realpath(Object::filename(id)))
          return File(id, true);
        }
        hid_t fapl = driver.fapl();
        hid_t id = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
        MPQC_FILE_VERIFY(id);
        files_::data.insert(id);
        return File(id, false);
        //return File(H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), false);
      }

      static File create(const std::string &name, const Driver &driver) {
        hid_t fapl = driver.fapl();
        hid_t id = H5Fcreate(name.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, fapl);
        MPQC_FILE_VERIFY(id);
        files_::data.insert(id);
        return File(id, false);
      }

  };

  template<typename T_>
  struct File::Dataspace: range_operator_base<Dataspace<T_>*, Dataspace<T_>,
      Dataspace<const T_> > {

      typedef range_operator_base<Dataspace*, Dataspace<T_>, Dataspace<const T_> > range_base;
      using range_base::operator();

      size_t size() const {
        return size_;
      }

      typedef typename boost::remove_const<T_>::type T;
      Dataspace operator[](size_t i) {
        std::vector<range> r = range_;
        r[ndims_ - 1] = range(i, i + 1);
        return Dataspace(parent_, base_, r, ndims_ - 1);
      }

      Dataspace<T> operator()(const std::vector<range> &r) {
        return Dataspace<T>(parent_, base_, extend(r), ndims_);
      }

      void write(const T *buffer) {
        MPQC_PROFILE_LINE;
        timer t;
        apply(&H5Dwrite, this->parent_.id(), rebase(range_), (T*) buffer);
        // printf("File::write size=%lu bytes, %f mb/s\n",
        //        (sizeof(T)*size_), (sizeof(T)*size_)/(t*1e6));
      }

      void read(T *buffer) {
        MPQC_PROFILE_LINE;
        timer t;
        apply(&H5Dread, this->parent_.id(), rebase(range_), buffer);
        // printf("File::read size=%lu bytes, rate=%f mb/s\n",
        //        (sizeof(T)*size_), (sizeof(T)*size_)/(t*1e6));
      }

      const std::vector<range>& extents() const {
        return range_;
      }

    private:
      Object parent_;
      size_t ndims_, size_;
      std::vector<size_t> base_;
      std::vector<range> range_;
      friend class Dataset<T> ;

      Dataspace(const Object &parent, const std::vector<size_t> &base,
                const std::vector<range> &r, size_t ndims) :
          parent_(parent), range_base(this), base_(base), range_(r), ndims_(
              ndims) {
        assert(ndims <= r.size());
        size_ = (r.size() ? 1 : 0);
        for (int i = 0; i < r.size(); ++i) {
          size_ *= r[i].size();
        }
      }

      std::vector<range> extend(const std::vector<range> &r) const {
        //std::cout << r.size() << " " << ndims_ << std::endl;
        assert(r.size() == ndims_);
        std::vector<range> x = r;
        for (size_t i = ndims_; i < range_.size(); ++i) {
          x.push_back(range_[i]);
        }
        return x;
      }

      std::vector<range> rebase(const std::vector<range> &r) const {
        assert(r.size() == base_.size());
        std::vector<range> v;
        for (int i = 0; i < base_.size(); ++i) {
          int begin = *r[i].begin() - base_[i];
          v.push_back(range(begin, begin + r[i].size()));
        }
        return v;
      }

      template<class F>
      static void apply(F f, hid_t dset, const std::vector<range> &r,
                        T *buffer) {
        hid_t fspace = H5Dget_space(dset);
        size_t size = select(fspace, r);
        hsize_t mdims[] = { size };
        hid_t mspace = H5Screate_simple(1, mdims, NULL);
        MPQC_FILE_VERIFY(mspace);
        MPQC_FILE_VERIFY(H5Sselect_all(mspace));
        hid_t type = file::h5_predtype<T>();
        timer t;
        MPQC_FILE_VERIFY(f(dset, type, mspace, fspace, H5P_DEFAULT, buffer));
        //std::cout << "hdf5 op rate: " << (size*sizeof(T)/1e6)/double(t) << " mb/s" << std::endl;
        MPQC_FILE_VERIFY(H5Sclose(mspace));
        MPQC_FILE_VERIFY(H5Sclose(fspace));
      }

      static size_t select(hid_t space, const std::vector<range> &r) {
        size_t N = H5Sget_simple_extent_ndims(space);
        MPQC_FILE_VERIFY(N);
        //printf("id = %i, N = %lu\n", space, N);
        hsize_t fstart[N];
        hsize_t fstride[N]; // Stride of hyperslab
        hsize_t fcount[N]; // Block count
        hsize_t fblock[N]; // Block sizes
        size_t size = 1;
        //printf("select [ ");
        for (size_t i = 0, j = N - 1; i < N; ++i, --j) {
          fstart[i] = *r[j].begin();
          fcount[i] = r[j].size();
          //printf("%i:%i,", fstart[i], fstart[i]+fcount[i]);
          fstride[i] = 1;
          fblock[i] = 1;
          size *= fcount[i];
        }
        //printf(" ], size = %i\n", size);
        MPQC_FILE_VERIFY(
            H5Sselect_hyperslab (space, H5S_SELECT_SET, fstart, fstride, fcount, fblock));
        return size;
      }

  };

  template<typename T>
  struct File::Dataset: File::Object, range_operator_base<Dataset<T>*,
      Dataspace<T>, Dataspace<const T> > {
      typedef range_operator_base<Dataset<T>*, Dataspace<T>, Dataspace<const T> > range_base;
      using range_base::operator();

      Dataset() :
          range_base(this) {
      }

      template<typename Extents>
      Dataset(const Object &parent, const std::string &name,
              const Extents &extents, const std::vector<size_t> &chunk =
                  std::vector<size_t>()) :
          Object(Dataset::create(parent, name, extents, chunk)), range_base(
              this) {
        assert(id() > 0);
foreach      (auto e, extents) {
        range r = extent(e);
        base_.push_back(*r.begin());
        dims_.push_back(r.size());
      }
    }

    size_t rank() const {
      return dims_.size();
    }

    Dataspace<T> operator[](size_t index) {
      std::vector<range> r;
      for (int i = 0; i < this->rank(); ++i) {
        r.push_back(range(base_[i], base_[i] + dims_[i]));
      }
      return Dataspace<T>(*this, base_, r, r.size())[index];
    }

    Dataspace<T> operator()(const std::vector<range> &r) {
      //std::cout << this->extents_.size() << " " << r.size() << std::endl;
      assert(this->rank() == r.size());
      return Dataspace<T>(*this, base_, r, r.size());
    }

    Dataspace<const T> operator()(const std::vector<range> &r) const {
      assert(this->rank() == r.size());
      return Dataspace<const T>(*this, base_, r, r.size());
    }

#define MPQC_FILE_DATASET_OPERATOR(Z, N, TEXT)                                  \
        template<BOOST_PP_ENUM_PARAMS(N, class R)>                              \
        Dataspace<T> operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const R, &r)) {  \
            return this->operator()(rangify(BOOST_PP_ENUM_PARAMS(N, r)));  \
        }

    MPQC_FILE_DATASET_OPERATOR(nil, 3, nil)

    template<typename Extents>
    static Object create(const Object &parent,
        const std::string &name,
        const Extents &extents,
        const std::vector<size_t> &chunk) {

      std::vector<hsize_t> dims;
      foreach (auto e, extents) {
        dims.push_back(extent(e).size());
      }
      std::reverse(dims.begin(), dims.end());

      hid_t fspace = H5Screate_simple(dims.size(), &dims[0], NULL);
      hid_t type = file::h5_predtype<T>();
      hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
      if (chunk.size()) {
        assert(chunk.size() == dims.size());
        std::vector<hsize_t> block(chunk.rbegin(), chunk.rend());
        H5Pset_chunk(dcpl, block.size(), &block[0]);
      }
      //printf("parent.id() = %i, name=%s\n", parent.id(), name.c_str());
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
      hid_t id = H5Dcreate
#else
      hid_t id = H5Dcreate1
#endif
      (parent.id(), name.c_str(), type, fspace, dcpl);
      MPQC_FILE_VERIFY(id > 0);
      H5Pclose(dcpl);
      return Object(parent, id, &Dataset::close, false);
    }

    // static Dataset open(Object parent, const std::string &name) {
    //     return Dataset(H5Dopen1(parent.id(), name.c_str()));
    // }

    static void close(hid_t id) {
      // printf("H5Dclose(%i), ref=%i\n", id, H5Iget_ref(id));
      // H5Dclose(id);
    }

    private:
    std::vector<size_t> base_, dims_;

    static range extent(const range &r) {
      return r;
    }

    template <typename E>
    static range extent(const E &e) {
      return range(0,e);
    }

  };

  struct File::Group: File::Object {
      static Group create(Object parent, const std::string &name) {
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
        hid_t id = H5Gopen(parent.id(), name.c_str());
#else
        hid_t id = H5Gopen1(parent.id(), name.c_str());
#endif
        return Group(Object(parent, id, &Group::close, false));
      }
      static void close(hid_t id) {
        //printf("H5Gclose(%i), ref=%i\n", id, H5Iget_ref(id));
        //H5Gclose(id);
      }
      template<typename T, typename Dims>
      Dataset<T> dataset(const std::string &name, const Dims &dims) {
        return Dataset<T>(*this, name, dims);
      }
    private:
      Group(Object o) :
          Object(o) {
      }
  };

  inline File::Group File::group(const std::string &name) {
    return Group::create(*this, name);
  }

  template<typename T, class A>
  void operator>>(File::Dataspace<T> ds, A &a) {
    ds.read(a.data());
  }

  template<typename T, class A>
  void operator<<(File::Dataspace<T> ds, const A &a) {
    ds.write(a.data());
  }

} // namespace mpqc

#endif // MPQC_FILE_HPP
