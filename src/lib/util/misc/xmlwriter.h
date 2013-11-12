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

#ifndef XMLWRITER_H_
#define XMLWRITER_H_

#include <algorithm>
#include <algorithm>
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/insert_linebreaks.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/archive/iterators/ostream_iterator.hpp>
#include <Eigen/Dense>

#include <util/misc/runnable.h>
#include <util/misc/exenv.h>
#include <util/misc/consumableresources.h>
#include <Eigen/Dense>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/type_traits.hpp>

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

  class SCVector3;
  class Grid;
  class Units;
  class XMLWritable;
  class SCVector;

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
  class XMLWriter : public Runnable {

    protected:

      std::ostream* out_;
      ptree pt_;

      bool delete_out_;
      bool compress_data_;
      bool pretty_print_;
      bool human_readable_;
      bool fold_in_class_name_;

      int pretty_print_spaces_;
      char pretty_print_space_char_;
      std::vector< Ref<XMLWritable> > data_;
      boost::property_tree::xml_writer_settings<char> write_settings_;

      void init();

      void init_filename(const std::string& filename);


    public:

      XMLWriter(const Ref<KeyVal>& keyval);

      XMLWriter(std::ostream& out=ExEnv::out0());

      XMLWriter(ptree& pt, std::ostream& out=ExEnv::out0());

      XMLWriter(const std::string& filename);

      virtual ~XMLWriter();

      bool fold_in_class_name() const { return fold_in_class_name_; }

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
          T obj,
          std::string wrapper_name=""
      ) const
      {
        //----------------------------------------//
        ptree* my_tree;
        if(not wrapper_name.empty()){
          // Wrap the object in a ptree named wrapper_name
          ptree& child_tree = parent.add_child(wrapper_name, ptree());
          my_tree = &child_tree;
        }
        else{
          my_tree = &parent;
        }
        //----------------------------------------//
        if(wrapper_name.empty()){
          return this->write_to_xml(obj, *my_tree);
        }
        else{
          this->write_to_xml(obj, *my_tree);
          return *my_tree;
        }
        //----------------------------------------//
      }

      void add_data(const Ref<XMLWritable>& datum) { data_.push_back(datum); }

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
      ptree& write_to_xml(const Ref<T>& obj, ptree& parent) const {
        return write_to_xml_impl(obj, parent, boost::is_base_of<XMLWritable, T>());
      }

    private:

      template<typename T>
      ptree& write_to_xml_impl(const Ref<T>& obj, ptree& parent, const boost::true_type is_writable) const {
          Ref<XMLWritable> obj_write = obj;
          return write_to_xml(obj_write, parent);
      }

      template<typename T>
      ptree& write_to_xml_impl(const Ref<T>& obj, ptree& parent, const boost::false_type is_writable) const {
          return write_to_xml(*obj, parent);
      }



  };

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

namespace boost {
  namespace property_tree {

    template<typename Ch, typename Traits, typename Alloc>
    struct translator_between<std::basic_string<Ch, Traits, Alloc>, sc::XMLDataStream<double> >
    {
      typedef sc::XMLDataStreamTranslator<double> type;
    };

  }
}



#endif /* XMLWRITER_H_ */
