//
// xml.h
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

#ifndef XML_H_
#define XML_H_
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <util/misc/exenv.h>
#include <util/misc/consumableresources.h>
//TODO add option to gzip data streams

namespace sc {

  class XMLWriter;

  ////////////////////////////////////////////////////////////////////////////////

  class XMLReader {

  };

  ////////////////////////////////////////////////////////////////////////////////

  class XMLWritable : virtual public RefCount {
    public:
      virtual boost::property_tree::ptree& write_xml(boost::property_tree::ptree& pt, const XMLWriter& writer) = 0;
      virtual ~XMLWritable() {}
  };

  class DescribedXMLWritable : public XMLWritable, virtual public DescribedClass {
    protected:
      boost::property_tree::ptree* my_ptree_ = 0;
      virtual boost::property_tree::ptree& get_my_ptree(boost::property_tree::ptree& parent, std::string name = "");
  };

  ////////////////////////////////////////////////////////////////////////////////

  class XMLReadable {
    public:
      virtual void read_xml(boost::property_tree::ptree& pt, const XMLReader& reader) = 0;
      virtual ~XMLReadable() {}
  };

  ////////////////////////////////////////////////////////////////////////////////

} // end namespace sc


#endif /* XML_H_ */
