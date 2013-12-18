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

// boost parameter library
#include <boost/parameter/name.hpp>
#include <boost/parameter/preprocessor.hpp>

// Eigen
#include <Eigen/Dense>

#include <util/misc/runnable.h>
#include <util/misc/exenv.h>
#include <util/misc/consumableresources.h>
#include <util/keyval/keyval.h>
#include <util/misc/scexception.h>
#include <util/misc/xml.h>

#ifndef NO_USE_BOOST_ENDIAN
#  include <boost/detail/endian.hpp>
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


using boost::property_tree::ptree;

namespace sc {

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
    BOOST_PARAMETER_NAME(attrs)
  }

  ////////////////////////////////////////////////////////////////////////////////

  class SCVector3;
  class Grid;
  class Units;
  class XMLWritable;
  class SCVector;

  ////////////////////////////////////////////////////////////////////////////////

  template<typename T>
  class XMLDataStream {
  private:
    T* data_;
    unsigned long n_;
    bool deallocate_when_destroyed_;
    bool human_readable_;
    bool pretty_print_;

  public:
    XMLDataStream(
        T* data,
        unsigned long n,
        bool human_readable=false,
        bool pretty_print=false,
        bool deallocate_when_destroyed=true
    ) :
      data_(data),
      n_(n),
      deallocate_when_destroyed_(deallocate_when_destroyed),
      human_readable_(human_readable),
      pretty_print_(pretty_print)
    {

    }

    ~XMLDataStream(){
      if(deallocate_when_destroyed_){
        deallocate(data_);
      }
    }

    unsigned long n() const { return n_; }
    T* data() const { return data_; }
    bool human_readable() const { return human_readable_; }
    bool pretty_print() const { return pretty_print_; }


  };

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
      boost::property_tree::xml_writer_settings<char> write_settings_;

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
          XMLDataStream<T> xds(data, ndata,
              /* human_readable = */ false,
              /* pretty_print = */ pretty_print_,
              /* deallocate_when_destroyed = */ !preserve_data_pointer
          );
          pt.put_value(xds);
        }
        else{
          pt.put("<xmlattr>.ndata", ndata);
          pt.put("<xmlattr>.human_readable", true);
          XMLDataStream<T> xds(data, ndata,
              /* human_readable = */ true,
              /* pretty_print = */ pretty_print_,
              /* deallocate_when_destroyed = */ !preserve_data_pointer
          );
          pt.put_value(xds);
        }
      }

      template <typename T>
      ptree& insert_child(
          ptree& parent,
          T obj
      ) const
      {
        return this->write_to_xml(obj, parent);
      }

      template <typename T>
      ptree& insert_child(
          ptree& parent,
          T obj,
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
          T obj,
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
          T obj,
          Args... args
      ) const
      {
        assert(not writing_done_);
        return insert_child(*current_root_, obj, args...);
      }


      void add_data(const Ref<XMLWritable>& datum) { data_.push_back(datum); }

      void do_write();

      void run();

      ptree& write_to_xml(const Ref<XMLWritable>& obj, ptree& parent) const;
      ptree& write_to_xml(XMLWritable& obj, ptree& parent) const;

      //----------------------------------------------------------------------------//
      // Non-intrusive interface for some types

      ptree& write_to_xml(const Eigen::VectorXd& obj, ptree& parent) const;
      ptree& write_to_xml(const Eigen::MatrixXd& obj, ptree& parent) const;
      ptree& write_to_xml(const SCVector3& obj, ptree& parent) const;
      ptree& write_to_xml(const SCVector& obj, ptree& parent) const;
      ptree& write_to_xml(const Grid& obj, ptree& parent) const;
      ptree& write_to_xml(const Units& obj, ptree& parent) const;

      template<typename T>
      inline typename boost::enable_if<boost::is_base_of<RefCount, T>, ptree&>::type
      write_to_xml(const Ref<T>& obj, ptree& parent) const {
        return write_to_xml_impl(obj, parent, boost::is_base_of<XMLWritable, T>());
      }

      /*template<typename T>
      inline typename boost::enable_if<boost::mpl::not_<boost::is_base_of<RefCount, T>>, ptree&>::type
      write_to_xml(T obj, ptree& parent) const {
        return write_to_xml(obj, parent);
      }*/

      template<typename T>
      inline ptree&
      write_to_xml(T obj) const {
        assert(not writing_done_);
        return write_to_xml(obj, *current_root_);
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
    if(XMLWriter::current_writer.nonnull() and not XMLWriter::current_writer->writing_done()){
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
        and XMLWriter::current_writer.nonnull()
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
  _write_as_xml_impl(T object, const std::string& tag_name, const MapType& attrs)
  {
     XMLWriter::current_writer->insert_child_default(object, tag_name);
  }

  template <typename T, typename MapType>
  inline typename boost::enable_if<boost::mpl::not_<is_convertible<MapType, bool>>, void>::type
  _write_as_xml_impl(T object, const std::string& tag_name, const MapType& attrs)
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
      (attrs, *, false)
    )
  )
  {
    if(XMLWriter::current_writer.null()) {
      throw ProgrammingError("No current XML context", __FILE__, __LINE__);
    }
    else{
      std::string tag_name = name;
      _write_as_xml_impl(object, tag_name, attrs);
    }
  }

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
      assert(false); // not implemented
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

    template<typename Ch, typename Traits, typename Alloc>
    struct translator_between<std::basic_string<Ch, Traits, Alloc>, sc::XMLDataStream<double> >
    {
      typedef sc::XMLDataStreamTranslator<double> type;
    };

  }
}



#endif /* _util_misc_xmlwriter_h */
