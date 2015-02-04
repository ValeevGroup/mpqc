//
// xml.cc
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



#include <util/misc/xml.h>
#include <util/misc/xmlwriter.h>

using namespace std;
using namespace sc;
using boost::property_tree::ptree;
using boost::property_tree::xml_writer_settings;


////////////////////////////////////////////////////////////////////////////////
// write_human_readable_data overloads

string
detail::get_human_readable_data(
    double* data,
    long ndata,
    int nperline,
    int precision,
    int width
)
{
  std::stringstream sstr;
  // Newline at the beginning
  if (ndata > 0) sstr << "\n";
  for(long i = 0; i < ndata; ++i){
    sstr << setfill(' ') << setw(width) << setprecision(precision) << data[i];
    if((i + 1) % nperline == 0 && i + 1 != ndata){
      sstr << std::endl;
    }
  }
  // Newline at the end
  if (ndata > 0) sstr << "\n";
  return sstr.str();
}

string
detail::get_human_readable_data(
    int* data,
    long ndata,
    int nperline,
    int precision, // ignored
    int width
)
{
  std::stringstream sstr;
  // Newline at the beginning
  if (ndata > 0) sstr << "\n";
  for(long i = 0; i < ndata; ++i){
    sstr << setfill(' ') << setw(width) << data[i];
    if((i + 1) % nperline == 0 && i + 1 != ndata){
      sstr << std::endl;
    }
  }
  // Newline at the end
  if (ndata > 0) sstr << "\n";
  return sstr.str();
}

////////////////////////////////////////////////////////////////////////////////

static ClassDesc DescribedXMLWritable_cd(
  typeid(DescribedXMLWritable), "DescribedXMLWritable", 1, "public XMLWritable, virtual public DescribedClass",
  0, 0, 0);

ptree&
DescribedXMLWritable::get_my_ptree(ptree& parent, std::string name)
{
  if(my_ptree_ == 0){

    ptree& new_tree = parent.add_child(
        name.empty() ? this->class_name() : name, ptree()
    );
    my_ptree_ = &new_tree;
  }
  return *my_ptree_;

}

////////////////////////////////////////////////////////////////////////////////

//boost::optional<std::string>
//XMLDataStreamTranslator<double>::put_value(const XMLDataStream<double> xds)
//{
//
//}



