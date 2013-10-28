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



#include <util/misc/xml.h>

using namespace std;
using namespace sc;
using boost::property_tree::ptree;
using boost::property_tree::xml_writer_settings;

static ClassDesc XMLWriter_cd(
  typeid(XMLWriter), "XMLWriter", 1, "public Runnable",
  0, create<XMLWriter>, 0);

XMLWriter::XMLWriter(const Ref<KeyVal>& keyval) :
    out_(0),
    pt_(ptree()),
    delete_out_(false),
    compress_data_(false),
    pretty_print_(false),
    pretty_print_spaces_(2),
    pretty_print_space_char_(' '),
    data_()
{
  compress_data_ = keyval->booleanvalue("compress_data", KeyValValueboolean(false));
  pretty_print_ = keyval->booleanvalue("pretty_print", KeyValValueboolean(true));
  pretty_print_spaces_ = keyval->intvalue("indent_spaces", KeyValValueint(2));
  pretty_print_space_char_ = keyval->booleanvalue("use_tabs", KeyValValueboolean(false)) ? '\t' : ' ';
  human_readable_ = keyval->booleanvalue("human_readable", KeyValValueboolean(false));
  string filename = keyval->stringvalue("filename", KeyValValuestring("-"));
  if(keyval->exists("data")){
    int key_count = keyval->count("data");
    for(int idat = 0; idat < key_count; ++idat){
      Ref<XMLWritable> to_write;
      to_write << keyval->describedclassvalue("data", idat);
      data_.push_back(to_write);
    }
  }
  assert(false);

  init_filename(filename);
  init();
}


XMLWriter::XMLWriter(ostream& out) :
    out_(&out),
    pt_(ptree()),
    delete_out_(false),
    compress_data_(false),
    pretty_print_(false),
    pretty_print_spaces_(2),
    pretty_print_space_char_(' '),
    data_()
{
  init();
}


XMLWriter::XMLWriter(ptree& pt, ostream& out) :
    out_(&out),
    pt_(pt),
    delete_out_(false),
    compress_data_(false),
    pretty_print_(false),
    pretty_print_spaces_(2),
    pretty_print_space_char_(' '),
    data_()
{
  init();
}


XMLWriter::XMLWriter(const string& filename) :
    out_(0),
    pt_(ptree()),
    delete_out_(false),
    compress_data_(false),
    pretty_print_(false),
    pretty_print_spaces_(2),
    pretty_print_space_char_(' '),
    data_()
{
  init_filename(filename);
  init();
}

void
XMLWriter::init_filename(const string& filename){
  // TODO delay opening of file pointer until the data is actually written
  if (filename == "-") {
      out_ = &(ExEnv::out0());
  }
  else {
      out_ = new std::ofstream(filename.c_str());
      delete_out_ = true;
  }
}

void
XMLWriter::init(){
  if(pretty_print_){
    write_settings_ = xml_writer_settings<char>(pretty_print_space_char_, pretty_print_spaces_);
  }
  else{
    write_settings_ = xml_writer_settings<char>();
  }
}


XMLWriter::~XMLWriter()
{
  if(delete_out_){
    delete out_;
  }
}

ptree&
XMLWriter::add_writable_child(
    ptree& parent,
    const std::string& name,
    const Ref<XMLWritable>& child
) const
{
  ptree& child_tree = parent.add_child(name, ptree());
  child->write_xml(child_tree, *this);
  return child_tree;
}

void
XMLWriter::run()
{
  std::vector< Ref<XMLWritable> >::iterator it;
  for(it = data_.begin(); it != data_.end(); ++it){
    (*it)->write_xml(pt_, *this);
  }
  boost::property_tree::write_xml(*out_, pt_, write_settings_);
}

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


