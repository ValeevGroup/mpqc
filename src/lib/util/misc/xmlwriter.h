//
// xmlwriter.h
//
// Copyright (C) 2013 David Hollman
//
// Author: David Hollman
//
// This file is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// This file is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with the MPQC; see the file COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#ifndef _util_misc_xmlwriter_h
#define _util_misc_xmlwriter_h

#define BOOST_PARAMETER_MAX_ARITY 15

#include <algorithm>
#include <stack>
#include <tuple>

// boost archive iterators for base64 conversion
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/insert_linebreaks.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/archive/iterators/ostream_iterator.hpp>

// boost::property_tree
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

// boost type traits and mpl
#include <boost/type_traits.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/and.hpp>

// boost parameter library
#include <boost/parameter/name.hpp>
#include <boost/parameter/preprocessor.hpp>
#include <boost/version.hpp>

// Eigen
#include <Eigen/Dense>

#include <util/misc/runnable.h>
#include <util/misc/exenv.h>
#include <util/misc/consumableresources.h>
#include <util/keyval/keyval.h>
#include <util/misc/scexception.h>
#include <util/misc/xml.h>

#ifndef NO_USE_BOOST_ENDIAN
#  ifdef __has_include
#    if __has_include(<boost/detail/endian.hpp>)
#     include <boost/detail/endian.hpp>
#    else
#       include <boost/endian.hpp>
#    endif
#  else
#     include <boost/endian.hpp>
#  endif
#  if defined(BOOST_LITTLE_ENDIAN)
#    define IS_BIG_ENDIAN false
#  elif defined(BOOST_BIG_ENDIAN)
#    define IS_BIG_ENDIAN true
#  else
#    include <arpa/inet.h>
#    define IS_BIG_ENDIAN htonl(47) == 47
#  endif
#else
#  include <arpa/inet.h>
#  define IS_BIG_ENDIAN htonl(47) == 47
#endif


namespace sc {

  using boost::property_tree::ptree;
#if BOOST_VERSION >= 105600
  typedef boost::property_tree::xml_writer_settings<std::string> xml_writer_settings;
#else
  typedef boost::property_tree::xml_writer_settings<char> xml_writer_settings;
#endif

  ////////////////////////////////////////////////////////////////////////////////

  namespace parameter {
    // Parameters for boost::parameter parameters are always wrapped in namespace
    //   sc::parameter in MPQC.

    // Parameters for generic function that quickly writes an object to a given
    //   XML file; extremely useful for debugging and matching old code.  These
    //   are basically just the keyval arguments to XMLWriter
    BOOST_PARAMETER_NAME(object)
    BOOST_PARAMETER_NAME(outfile)
    BOOST_PARAMETER_NAME(out)
    BOOST_PARAMETER_NAME(compress_data)
    BOOST_PARAMETER_NAME(pretty_print)
    BOOST_PARAMETER_NAME(indent_spaces)
    BOOST_PARAMETER_NAME(use_tabs)
    BOOST_PARAMETER_NAME(human_readable)
    BOOST_PARAMETER_NAME(fold_in_class_names)
    BOOST_PARAMETER_NAME(root_name)
    BOOST_PARAMETER_NAME(name)
    BOOST_PARAMETER_NAME(attributes)
  }

  ////////////////////////////////////////////////////////////////////////////////

  class SCVector3;
  class Grid;
  class Units;
  class XMLWritable;
  class SCVector;

  ////////////////////////////////////////////////////////////////////////////////

  namespace {
    template<typename T, bool do_it>
    struct destroy_data {
        void operator()(T* data) const { };
    };

    template<typename T>
    struct destroy_data<T, true> {
      void operator()(T* data) const {
        deallocate(data);
      }
    };

  }

  template<typename T, bool deallocate_when_destroyed=not std::is_const<T>::value>
  class XMLDataStream {
  private:
    T* data_;
    unsigned long n_;
    bool human_readable_;
    bool pretty_print_;

  public:
    XMLDataStream(
        T* data,
        unsigned long n,
        bool human_readable=false,
        bool pretty_print=false
    ) :
      data_(data),
      n_(n),
      human_readable_(human_readable),
      pretty_print_(pretty_print)
    { }

    ~XMLDataStream() {
      destroy_data<T, deallocate_when_destroyed>()(data_);
    }

    unsigned long n() const { return n_; }
    T* data() const { return data_; }
    bool human_readable() const { return human_readable_; }
    bool pretty_print() const { return pretty_print_; }

  private:

  };


  ////////////////////////////////////////////////////////////////////////////////

  class XMLWriter;

  template <typename T>
  typename boost::enable_if<
    boost::is_base_of<XMLWritable, typename boost::decay<T>::type>,
    ptree&
  >::type
  write_xml(
      const Ref<T>& obj,
      ptree& parent,
      const XMLWriter& writer
  );

  template <typename T>
  ptree&
  write_xml(
      const Ref<T>& obj,
      typename boost::disable_if_c<
        boost::is_base_of<XMLWritable, typename boost::decay<T>::type>::value
        or not boost::is_base_of<RefCount, typename boost::decay<T>::type>::value,
        ptree&
      >::type const& parent,
      const XMLWriter& writer
  );

  ptree& write_xml(XMLWritable&, ptree&, const XMLWriter&);
  ptree& write_xml(const SCVector3&, ptree&, const XMLWriter&);
  ptree& write_xml(const SCVector&, ptree&, const XMLWriter&);
  ptree& write_xml(const Grid&, ptree&, const XMLWriter&);
  ptree& write_xml(const Units&, ptree&, const XMLWriter&);
  ptree& write_xml(const Eigen::MatrixXd&, ptree&, const XMLWriter&);
  ptree& write_xml(const Eigen::VectorXd&, ptree&, const XMLWriter&);
  ptree& write_xml(const std::vector<double>&, ptree&, const XMLWriter&);

  template<typename Derived>
  ptree& write_xml(const Eigen::MatrixBase<Derived>&, ptree&, const XMLWriter&);

  ////////////////////////////////////////////////////////////////////////////////

  // TODO gzip compression, turned on and off by the writer
  // TODO write as data is received?
  // TODO perfect forwarding of objects; make sure copies aren't happening
  class XMLWriter : public Runnable {

    protected:

      std::ostream* out_;

      ptree* current_root_;
      std::stack<ptree*> pt_stack_;

      bool delete_out_;
      bool delete_pt_ = false;
      bool compress_data_;
      bool pretty_print_;
      bool human_readable_;
      bool fold_in_class_name_;
      bool writing_done_ = false;

      int pretty_print_spaces_;
      char pretty_print_space_char_;
      std::vector< Ref<XMLWritable> > data_;
      xml_writer_settings write_settings_;

      void init();

      void init_filename(const std::string& filename);

      std::string filename_;

    public:

      XMLWriter() = delete;

      XMLWriter(const XMLWriter&) = delete;

      XMLWriter(const Ref<KeyVal>& keyval);

      XMLWriter(std::ostream& out=ExEnv::out0());

      XMLWriter(ptree* pt, std::ostream& out=ExEnv::out0());

      XMLWriter(const std::string& filename);

      virtual ~XMLWriter();

      bool fold_in_class_name() const { return fold_in_class_name_; }
      const std::string& filename() const { return filename_; }
      bool writing_done() const { return writing_done_; }
      bool has_active_context() const { return pt_stack_.size() > 1; }

      template<typename T>
      void
      put_binary_data(
          ptree& pt,
          T* data,
          long ndata,
          bool preserve_data_pointer = false
      ) const
      {
        assert(!compress_data_); // Not implemented
        if(!human_readable_){
          pt.put("<xmlattr>.ndata", ndata);
          pt.put("<xmlattr>.datum_size", sizeof(T));
          pt.put("<xmlattr>.big_endian", IS_BIG_ENDIAN);
        }
        else{
          pt.put("<xmlattr>.ndata", ndata);
          pt.put("<xmlattr>.human_readable", true);
        }

        if(preserve_data_pointer) {
          XMLDataStream<const T> xds(const_cast<const T*>(data), ndata,
              /* human_readable = */ human_readable_,
              /* pretty_print = */ pretty_print_
          );
          pt.put_value(xds);
        }
        else {
          XMLDataStream<T> xds(data, ndata,
              /* human_readable = */ human_readable_,
              /* pretty_print = */ pretty_print_
          );
          pt.put_value(xds);
        }
      }

      template <typename T>
      ptree& insert_child(
          ptree& parent,
          T&& obj
      ) const
      {
        return this->write_to_xml(obj, parent);
      }

      template <typename T>
      ptree& insert_child(
          ptree& parent,
          T&& obj,
          std::string wrapper_name
      ) const
      {
        // Wrap the object in a ptree named wrapper_name
        ptree& child_tree = parent.add_child(wrapper_name, ptree());
        this->write_to_xml(obj, child_tree);
        return child_tree;
      }

      template <typename T, typename MapType>
      ptree& insert_child(
          ptree& parent,
          T&& obj,
          std::string wrapper_name,
          const MapType& attrs
      ) const
      {
        // Wrap the object in a ptree named wrapper_name
        ptree& child_tree = parent.add_child(wrapper_name, ptree());
        this->write_to_xml(obj, child_tree);
        for(auto it : attrs) {
          child_tree.put("<xmlattr>." + it.first, it.second);
        }
        return child_tree;
      }

      template <typename T, typename... Args>
      ptree& insert_child_default(
          T&& obj,
          Args... args
      ) const
      {
        assert(not writing_done_);
        return insert_child(*current_root_, obj, args...);
      }

      void add_data(const Ref<XMLWritable>& datum) { data_.push_back(datum); }

      void do_write();

      void run();

      template <typename T>
      ptree& write_to_xml(T&& obj, ptree& parent) const
      {
        // Call the non-intrusive interface
        return write_xml(std::forward<T>(obj), parent, *this);
      }

      template<typename T>
      inline typename boost::enable_if<boost::is_base_of<RefCount, T>, ptree&>::type
      write_to_xml(const Ref<T>& obj, ptree& parent) const {
        return write_to_xml_impl(obj, parent, boost::is_base_of<XMLWritable, T>());
      }

      template<typename T>
      ptree& write_to_xml(T&& obj) const {
        assert(not writing_done_);
        return write_to_xml(std::forward<T>(obj), *current_root_);
      }

      //----------------------------------------------------------------------------//
      // Transparent context for writing objects

      void begin_writing_context(const std::string& root_name);

      void end_writing_context();

      static Ref<XMLWriter> current_writer;
      static std::stack<Ref<XMLWriter>> writer_stack;
      static std::string current_context_name;
      static std::stack<std::string> context_name_stack;

    private:

      template<typename T>
      inline ptree& write_to_xml_impl(const Ref<T>& obj, ptree& parent, const boost::true_type is_writable) const {
          Ref<XMLWritable> obj_write = obj;
          return write_to_xml(obj_write, parent);
      }

      template<typename T>
      inline ptree& write_to_xml_impl(const Ref<T>& obj, ptree& parent, const boost::false_type is_writable) const {
          return write_to_xml(*obj, parent);
      }


  };

  ////////////////////////////////////////////////////////////////////////////////

  template <typename T>
  typename boost::enable_if<
    boost::is_base_of<XMLWritable, typename boost::decay<T>::type>,
    ptree&
  >::type
  write_xml(
      const Ref<T>& obj,
      ptree& parent,
      const XMLWriter& writer
  )
  {
    return obj->write_xml(parent, writer);
  }

  template <typename T>
  ptree&
  write_xml(
      const Ref<T>& obj,
      typename boost::disable_if_c<
        boost::is_base_of<XMLWritable, typename boost::decay<T>::type>::value
        or not boost::is_base_of<RefCount, typename boost::decay<T>::type>::value,
        ptree&
      >::type const& parent,
      const XMLWriter& writer
  )
  {
    return write_xml(*obj, parent, writer);
  }

  template<typename Derived>
  ptree& write_xml(const Eigen::MatrixBase<Derived>& obj, ptree& parent, const XMLWriter& writer)
  {
    typedef Eigen::MatrixBase<Derived> MatrixType;
    typedef typename Eigen::internal::traits<Derived>::Scalar Scalar;

    ptree& child = parent.add_child("EigenDerived", ptree());
    const int ninner = obj.innerSize();
    const int nouter = obj.outerSize();
    child.put("<xmlattr>.ninner", ninner);
    child.put("<xmlattr>.nouter", nouter);
    child.put("<xmlattr>.row_major", int(MatrixType::IsRowMajor));
    child.put("<xmlattr>.is_vector", int(MatrixType::IsVectorAtCompileTime));
    child.put("<xmlattr>.signed",
        int(std::is_signed<Scalar>::value)
    );

    // Just iterate over everything so we don't have to think about strides and such
    Scalar* data = allocate<Scalar>(nouter*ninner);
    for(int i = 0; i < nouter; ++i){
      for(int j = 0; j < ninner; ++j){
        if(MatrixType::IsRowMajor) {
          data[i*ninner + j] = obj(i, j);
        }
        else {
          data[i*ninner + j] = obj(j, i);
        }
      }
    }

    // Note: the XMLDataStream created by put_binary_data now owns
    //   the pointer 'data'
    writer.put_binary_data(child.add_child("data", ptree()), data, nouter*ninner);
    return child;
  }

  namespace detail {

    template<typename Derived, unsigned int ViewMode>
    struct TriangleWriter { };

    template<typename Derived>
    struct TriangleWriter<Derived, Eigen::Lower> {
      ptree& operator()(
          const Eigen::TriangularView<Derived, Eigen::Lower>& obj,
          ptree& pt,
          const XMLWriter& writer
      ) const
      {
        const int nrows = obj.rows();
        const int ncols = obj.cols();
        if(nrows == ncols) {
          ptree& my_tree = pt.add_child("EigenDerived", ptree());
          my_tree.put("<xmlattr>.n", nrows);
          my_tree.put("<xmlattr>.lower_triangle", true);
          ptree& data_tree = my_tree.add_child("data", ptree());
          long ndata = nrows * (nrows+1) / 2;
          double* data = allocate<double>(ndata);
          double* data_spot = data;
          for(int irow = 0; irow < nrows; ++irow) {
            for(int icol = 0; icol <= irow; ++icol) {
              (*data_spot) = obj.coeff(irow, icol);
              ++data_spot;
            }
          }
          // The XMLDataStream object created by this function call
          //   owns the pointer data after this.
          writer.put_binary_data<double>(data_tree, data, ndata);
          return my_tree;
        }
        else {
          typedef Derived MatrixType;
          // Just write the unraveled version

          ptree& child = pt.add_child("EigenDerived", ptree());
          child.put("<xmlattr>.ninner", ncols);
          child.put("<xmlattr>.nouter", nrows);
          child.put("<xmlattr>.row_major", true);
          child.put("<xmlattr>.is_vector", int(MatrixType::IsVectorAtCompileTime));
          ptree& data_tree = child.add_child("data", ptree());
          double* data = allocate<double>(nrows*ncols);
          double* data_spot = data;
          for(int irow = 0; irow < nrows; ++irow) {
            for(int icol = 0; icol < ncols; ++icol) {
              if(icol <= irow) {
                (*data_spot) = obj.coeff(irow, icol);
              }
              else if(icol < nrows) {
                (*data_spot) = obj.coeff(icol, irow);
              }
              else {
                (*data_spot) = 0.0;
              }
              ++data_spot;
            }
          }
          // The XMLDataStream object created by this function call
          //   owns the pointer data after this.
          writer.put_binary_data<double>(data_tree, data, nrows*ncols);
          return child;
        }
      }

    };

  }

  template<typename Derived, unsigned int ViewMode>
  ptree& write_xml(const Eigen::TriangularView<Derived, ViewMode>& obj, ptree& parent, const XMLWriter& writer) {
    return detail::TriangleWriter<Derived, ViewMode>()(obj, parent, writer);
  }


  template<template<typename...> class Container, template<typename...> class TupleType, typename... NumTypes>
  typename boost::enable_if_c<
    boost::is_convertible<
      Container<TupleType<NumTypes...>>,
      std::vector<TupleType<NumTypes...>>
    >::value
    and
    boost::mpl::and_<
      boost::mpl::or_<
        boost::is_integral<NumTypes>,
        boost::is_floating_point<NumTypes>
      >...
    >::value,
    ptree&
  >::type
  write_xml(
      const Container<TupleType<NumTypes...>>& obj,
      ptree& parent,
      const XMLWriter& writer
  )
  {
    typedef Container<TupleType<NumTypes...>> input_type;
    typedef std::vector<TupleType<NumTypes...>> vect_type;

    // Cheesy but effective way to do automatic type conversion without
    //   unnecessary copies
    const auto& the_function = [](
      const vect_type& obj,
      ptree& parent,
      const XMLWriter& writer
    ) -> ptree&
    {
      ptree* child_ptr;
      if(writer.fold_in_class_name()) {
        parent.put("<xmlattr>.type", "std::vector<double>");
        child_ptr = &(parent);
      }
      else{
        ptree& tmp = parent.add_child("data", ptree());
        child_ptr = &tmp;
      }
      ptree& child = *child_ptr;
      writer.put_binary_data(child, obj.data(), obj.size(), true);
      child.put("<xmlattr>.nperdatum", sizeof...(NumTypes));
      return child;
    };

    return the_function(obj, parent, writer);
  }

  ////////////////////////////////////////////////////////////////////////////////



  using boost::is_convertible;
  namespace mpl = boost::mpl;

  #define XMLWRITER_FILENAME_NOT_GIVEN "__FILENAME_NOT_GIVEN__"

  // optional parameters mirroring the keyval contructor parameters for XMLWriter
  #define XMLWRITER_KV_OPTIONAL_PARAMS \
    (outfile, *(is_convertible<mpl::_, std::string>), XMLWRITER_FILENAME_NOT_GIVEN) \
    (compress_data, (bool), false) \
    (pretty_print, (bool), true) \
    (indent_spaces, (int), 2) \
    (use_tabs, (bool), false) \
    (human_readable, (bool), false) \
    (fold_in_class_names, (bool), true)

  #define XMLWRITER_CREATE_FROM_PARAMS(VARNAME) \
    std::string __fname = outfile; \
    Ref<AssignedKeyVal> akv = new AssignedKeyVal; \
    akv->assign("filename", \
        __fname == XMLWRITER_FILENAME_NOT_GIVEN ? "-" : __fname \
    ); \
    akv->assignboolean("compress_data", (bool)compress_data); \
    akv->assignboolean("pretty_print", (bool)pretty_print); \
    akv->assignboolean("indent_spaces", (bool)indent_spaces); \
    akv->assignboolean("use_tabs", (bool)use_tabs); \
    akv->assignboolean("human_readable", (bool)human_readable); \
    akv->assignboolean("fold_in_class_names", (bool)fold_in_class_names); \
    Ref<XMLWriter> VARNAME(new XMLWriter(akv));

  BOOST_PARAMETER_FUNCTION(
    (void),                         // Return type
    write_to_xml_file,              // Name of the function
    sc::parameter::tag,             // Namespace of the tag types

    (required                       // required parameters
      (object, *)                   // an object of any type.  Will not compile if the
                                    //   object can't be written to XML by some means
    )
    (optional
       XMLWRITER_KV_OPTIONAL_PARAMS
       (root_name,                  // root_name is the name of the root node on the xml output
           *(is_convertible<mpl::_, std::string>), std::string("mpqc")
       )
    )
  )
  {

    XMLWRITER_CREATE_FROM_PARAMS(writer)

    writer->begin_writing_context(std::string(root_name));
    writer->write_to_xml(object);
    writer->end_writing_context();
  }


  // Transparent context for writing
  BOOST_PARAMETER_FUNCTION(
    (void),
    begin_xml_context,
    parameter::tag,
    (required
      (name, *(is_convertible<mpl::_, std::string>))
    )
    (optional
      XMLWRITER_KV_OPTIONAL_PARAMS
    )
  )
  {
    std::string fname = outfile;
    std::string tag_name = name;
    bool fname_not_given = fname == XMLWRITER_FILENAME_NOT_GIVEN;
    if(!XMLWriter::current_writer.null() and not XMLWriter::current_writer->writing_done()){
      if(fname_not_given || XMLWriter::current_writer->filename() == fname){
        XMLWriter::current_writer->begin_writing_context(tag_name);
      }
      else {
        // Save the current writer, switch to new writer
        XMLWriter::writer_stack.push(XMLWriter::current_writer);

        XMLWRITER_CREATE_FROM_PARAMS(writer)

        XMLWriter::current_writer = writer;
        writer->begin_writing_context(tag_name);
      }
    }
    else {
      // TODO fail if there is an unclosed inner context that doesn't make sense
      XMLWRITER_CREATE_FROM_PARAMS(writer)

      XMLWriter::current_writer = writer;
      writer->begin_writing_context(tag_name);
    }

    XMLWriter::context_name_stack.push(XMLWriter::current_context_name);
    XMLWriter::current_context_name = tag_name;
  }

  // Transparent context for writing
  BOOST_PARAMETER_FUNCTION(
    (void),
    end_xml_context,
    parameter::tag,
    (optional
      (name, *(is_convertible<mpl::_, std::string>), XMLWRITER_FILENAME_NOT_GIVEN)
    )
  )
  {
    std::string tag_name = name;
    if(
        (tag_name == XMLWRITER_FILENAME_NOT_GIVEN or tag_name == XMLWriter::current_context_name)
        and !XMLWriter::current_writer.null()
        and not XMLWriter::context_name_stack.empty()
    ) {
      XMLWriter::current_writer->end_writing_context();
      if(not XMLWriter::current_writer->has_active_context()){
        if(XMLWriter::writer_stack.size() > 0) {
          XMLWriter::current_writer = XMLWriter::writer_stack.top();
          XMLWriter::writer_stack.pop();
        }
        else{
          XMLWriter::current_writer = 0;
        }
      }
    }
    else {
      throw ProgrammingError("Mismatched transparent XML contexts", __FILE__, __LINE__);
    }

    XMLWriter::current_context_name = XMLWriter::context_name_stack.top();
    XMLWriter::context_name_stack.pop();
  }

  template <typename T, typename MapType>
  inline typename boost::enable_if<is_convertible<MapType, bool>, void>::type
  _write_as_xml_impl(T&& object, const std::string& tag_name, const MapType& attrs)
  {
     XMLWriter::current_writer->insert_child_default(object, tag_name);
  }

  template <typename T, typename MapType>
  inline typename boost::enable_if<boost::mpl::not_<is_convertible<MapType, bool>>, void>::type
  _write_as_xml_impl(T&& object, const std::string& tag_name, const MapType& attrs)
  {
     XMLWriter::current_writer->insert_child_default(object, tag_name, attrs);
  }

  BOOST_PARAMETER_FUNCTION(
    (void),
    write_as_xml,
    parameter::tag,
    (required
      (name, *)
      (object, *)
    )
    (optional
      (attributes, *, false)
    )
  )
  {
    if(XMLWriter::current_writer.null()) {
      throw ProgrammingError("No current XML context", __FILE__, __LINE__);
    }
    else{
      std::string tag_name = name;
      _write_as_xml_impl(object, tag_name, attributes);
    }
  }

  template<typename Value=std::string>
  using attrs = std::map<std::string, Value>;

  ////////////////////////////////////////////////////////////////////////////////

  namespace detail {
    std::string
    get_human_readable_data(
        double* data,
        long ndata,
        int nperline=8,
        int precision=8,
        int width=15
    );

    std::string
    get_human_readable_data(
        int* data,
        long ndata,
        int nperline=8,
        int precision=8, // ignored
        int width=15
    );

    template<typename T>
    std::string
    get_human_readable_data(
        T* data,
        long ndata,
        int nperline=8,
        int precision=8,
        int width=15
    )
    {
      throw FeatureNotImplemented("get_human_readable_data", __FILE__, __LINE__);
      return "  "; // unreachable
    }
  }

  ////////////////////////////////////////////////////////////////////////////////

  template<typename T>
  struct XMLDataStreamTranslator {
    typedef std::string internal_type;
    typedef XMLDataStream<T> external_type;

    boost::optional<external_type> get_value(const internal_type& str)
    {
      if (!str.empty()) {
        // Require the string to be aligned to the sizeof(T)
        assert(str.length()*sizeof(char) % sizeof(T) == 0);

        // Transfer the data to it's new home
        unsigned long ndata = str.length() * sizeof(char) / sizeof(T);
        T* rv_data = new T[ndata];
        ::memcpy(rv_data, &(str.c_str()[0]), ndata * sizeof(T));
        return boost::optional<external_type>(external_type(rv_data, ndata));
      }
      else {
        return boost::optional<external_type>(boost::none);
      }

    }

    boost::optional<internal_type> put_value(const external_type& xds)
    {
      if(xds.human_readable()){
        return boost::optional<internal_type>(
            detail::get_human_readable_data(
                xds.data(),
                xds.n()
            )
        );

      }
      else{
        using namespace boost::archive::iterators;

        typedef
          base64_from_binary<
            transform_width<const char*, 6, 8>
          >
          base64_text;
        typedef insert_linebreaks<base64_text, 72> base64_text_linebreaks;

        std::stringstream sstr;
        if(xds.pretty_print()){
          std::copy(
              base64_text_linebreaks(xds.data()),
              base64_text_linebreaks(xds.data() + xds.n()),
              boost::archive::iterators::ostream_iterator<char>(sstr)
          );
          if(xds.n() * 8 * sizeof(T) % 6 == 2)
            return boost::optional<internal_type>(internal_type("\n" + sstr.str() + "==\n"));
          else if(xds.n() * 8 * sizeof(T) % 6 == 4)
            return boost::optional<internal_type>(internal_type("\n" + sstr.str() + "=\n"));
          else
            return boost::optional<internal_type>(internal_type("\n" + sstr.str() + "\n"));
        }
        else{
          std::copy(
              base64_text(xds.data()),
              base64_text(xds.data() + xds.n()),
              boost::archive::iterators::ostream_iterator<char>(sstr)
          );
          if(xds.n() * 8 * sizeof(T) % 6 == 2)
            return boost::optional<internal_type>(internal_type(sstr.str() + "=="));
          else if(xds.n() * 8 * sizeof(T) % 6 == 4)
            return boost::optional<internal_type>(internal_type(sstr.str() + "="));
          else
            return boost::optional<internal_type>(internal_type(sstr.str()));
        }
      }

    }
  };

} // end namespace sc

////////////////////////////////////////////////////////////////////////////////

namespace boost {
  namespace property_tree {

    template<typename Ch, typename Traits, typename Alloc, typename T, bool val>
    struct translator_between<std::basic_string<Ch, Traits, Alloc>, sc::XMLDataStream<T, val> >
    {
      typedef sc::XMLDataStreamTranslator<T> type;
    };

    //template<typename Ch, typename Traits, typename Alloc, bool val>
    //struct translator_between<std::basic_string<Ch, Traits, Alloc>, sc::XMLDataStream<const double, val> >
    //{
    //  typedef sc::XMLDataStreamTranslator<const double> type;
    //};



  }
}



#endif /* _util_misc_xmlwriter_h */
