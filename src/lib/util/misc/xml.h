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

//TODO add option to gzip data streams

using boost::property_tree::ptree;

namespace sc {

  template<typename T>
  inline void data_to_base64_chars(T* indata, char* outdata, int ndata, bool align=true, char fill_value=0){
    ::memcpy(outdata, indata, sizeof(T)/sizeof(char)*ndata);
    int nfill = sizeof(T)*ndata % 8*sizeof(char);
    if(align and nfill != 0){
      ::memset(outdata + (sizeof(T)*ndata / 8) * 8 * sizeof(char), fill_value, nfill);
    }
    else{
      assert(false); // not implemented
    }
  }

  class XMLWritable {
    public:
      virtual void write_xml(ptree& pt) = 0;
      virtual ~XMLWritable() {}
  };

  class XMLReadable {
    public:
      virtual void read_xml(ptree& pt) = 0;
      virtual ~XMLReadable() {}
  };

  class XMLWriter {

    protected:

      std::ostream* out_;
      ptree pt_;

      void init();

    public:

      XMLWriter(std::ostream& out=ExEnv::out0());

      XMLWriter(ptree& pt, std::ostream& out=ExEnv::out0());

      XMLWriter(std::string filename);

  };

  template<typename T>
  class XMLDataStream {
  private:
    T* data_;
    unsigned long n_;
    bool deallocate_when_destroyed_;

  public:
    XMLDataStream(
        T* data,
        unsigned long n,
        bool deallocate_when_destroyed=true
    ) : data_(data), n_(n), deallocate_when_destroyed_(deallocate_when_destroyed)
    {
    }

    ~XMLDataStream(){
      if(deallocate_when_destroyed_){
        deallocate(data_);
      }
    }

    unsigned long n() const { return n_; }
    T* data() const { return data_; }

  };

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
      return boost::optional<internal_type>(
          internal_type(
              (char*)xds.data(),
              xds.n()*sizeof(T)/sizeof(char)
          )
      );

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


#endif /* XML_H_ */
