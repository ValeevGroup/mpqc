//
// xml.h
//
// Copyright (C) 2013 MPQC authors
//
// Author: David Hollman <david.s.hollman@gmail.com
// Maintainer: DSH
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef XML_H_
#define XML_H_

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <util/misc/exenv.h>
#include <util/misc/consumableresources.h>
#include <util/misc/runnable.h>
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

//TODO add option to gzip data streams

using boost::property_tree::ptree;

namespace sc {

  // Forward Declarations
  class XMLWritable;
  class XMLReadable;

  ////////////////////////////////////////////////////////////////////////////////

  template<typename T>
  class XMLDataStream {
  private:
    T* data_;
    unsigned long n_;
    bool deallocate_when_destroyed_;
    bool human_readable_;

  public:
    XMLDataStream(
        T* data,
        unsigned long n,
        bool human_readable=false,
        bool deallocate_when_destroyed=true
    ) :
      data_(data),
      n_(n),
      deallocate_when_destroyed_(deallocate_when_destroyed),
      human_readable_(human_readable)
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

      template<typename T>
      void
      put_binary_data(
          ptree& pt,
          T* data,
          long ndata
      ) const
      {
        assert(!compress_data_); // Not implemented
        if(!human_readable_){
          pt.put("<xmlattr>.ndata", ndata);
          pt.put("<xmlattr>.datum_size", sizeof(T));
          pt.put("<xmlattr>.big_endian", IS_BIG_ENDIAN);
          XMLDataStream<T> xds(data, ndata);
          pt.put_value(xds);
        }
        else{
          pt.put("<xmlattr>.ndata", ndata);
          pt.put("<xmlattr>.human_readable", true);
          XMLDataStream<T> xds(data, ndata,
              /* human_readable = */ true
          );
          pt.put_value(xds);
        }
      }

      ptree&
      add_writable_child(
          ptree& parent,
          const std::string& name,
          const Ref<XMLWritable>& child
      ) const;


      void add_data(const Ref<XMLWritable>& datum) { data_.push_back(datum); }

      void run();


  };

  ////////////////////////////////////////////////////////////////////////////////

  class XMLReader {

  };

  ////////////////////////////////////////////////////////////////////////////////

  class XMLWritable : virtual public RefCount {
    public:
      virtual void write_xml(ptree& pt, const XMLWriter& writer) = 0;
      virtual ~XMLWritable() {}
  };

  class DescribedXMLWritable : public XMLWritable, virtual public DescribedClass {
    protected:
      ptree* my_ptree_ = 0;
      virtual ptree& get_my_ptree(ptree& parent, std::string name = "");
  };

  ////////////////////////////////////////////////////////////////////////////////

  class XMLReadable {
    public:
      virtual void read_xml(ptree& pt, const XMLReader& reader) = 0;
      virtual ~XMLReadable() {}
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
        return boost::optional<internal_type>(
            internal_type(
                (char*)xds.data(),
                xds.n()*sizeof(T)/sizeof(char)
            )
        );
      }

    }
  };

  ////////////////////////////////////////////////////////////////////////////////

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


#endif /* XML_H_ */
