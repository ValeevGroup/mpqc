//
// xmlwriter.cc
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

#include <util/misc/xmlwriter.h>
#include <util/misc/units.h>
#include <util/misc/xml.h>
#include <math/mmisc/grid.h>
#include <math/scmat/vector3.h>
#include <math/scmat/abstract.h>
#include <util/misc/consumableresources.h>

using namespace sc;
using namespace std;
using boost::property_tree::ptree;
using boost::property_tree::xml_writer_settings;

static ClassDesc XMLWriter_cd(
  typeid(XMLWriter), "XMLWriter", 1, "public Runnable",
  0, create<XMLWriter>, 0);

Ref<XMLWriter> XMLWriter::current_writer = 0;
std::stack<Ref<XMLWriter>> XMLWriter::writer_stack = {};
std::string XMLWriter::current_context_name = "";
std::stack<std::string> XMLWriter::context_name_stack = {};

XMLWriter::XMLWriter(const Ref<KeyVal>& keyval) :
    out_(0),
    current_root_(new ptree()),
    delete_out_(false),
    compress_data_(false),
    pretty_print_(false),
    data_(),
    delete_pt_(true)
{
  compress_data_ = keyval->booleanvalue("compress_data", KeyValValueboolean(false));
  pretty_print_ = keyval->booleanvalue("pretty_print", KeyValValueboolean(true));
  pretty_print_spaces_ = keyval->intvalue("indent_spaces", KeyValValueint(2));
  pretty_print_space_char_ = keyval->booleanvalue("use_tabs", KeyValValueboolean(false)) ? '\t' : ' ';
  human_readable_ = keyval->booleanvalue("human_readable", KeyValValueboolean(false));
  // Whenever reasonable, don't put class names on a seperate level; instead,
  //   use the type attribute.  This makes "manual" parsing easier by reducing
  //   the number of levels.  Might deprecate a false value at some point
  fold_in_class_name_ = keyval->booleanvalue("fold_in_class_names", KeyValValueboolean(true));

  string filename = keyval->stringvalue("filename", KeyValValuestring("-"));
  if(keyval->exists("data")){
    int key_count = keyval->count("data");
    for(int idat = 0; idat < key_count; ++idat){
      Ref<XMLWritable> to_write;
      to_write << keyval->describedclassvalue("data", idat);
      data_.push_back(to_write);
    }
  }

  init_filename(filename);
  init();
}

XMLWriter::XMLWriter(ostream& out) :
    out_(&out),
    current_root_(new ptree()),
    delete_out_(false),
    compress_data_(false),
    pretty_print_(false),
    pretty_print_spaces_(2),
    pretty_print_space_char_(' '),
    data_(),
    delete_pt_(true)
{
  init();
}

XMLWriter::XMLWriter(ptree* pt, ostream& out) :
    out_(&out),
    current_root_(pt),
    delete_out_(false),
    compress_data_(false),
    pretty_print_(false),
    pretty_print_spaces_(2),
    pretty_print_space_char_(' '),
    filename_(XMLWRITER_FILENAME_NOT_GIVEN),
    data_(),
    delete_pt_(false)
{
  init();
}

XMLWriter::XMLWriter(const string& filename) :
    out_(0),
    current_root_(new ptree()),
    delete_out_(false),
    compress_data_(false),
    pretty_print_(false),
    pretty_print_spaces_(2),
    pretty_print_space_char_(' '),
    data_(),
    delete_pt_(true)
{
  init_filename(filename);
  init();
}

void
XMLWriter::init_filename(const string& filename){
  filename_ = filename;
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

  pt_stack_.push(current_root_);
}

XMLWriter::~XMLWriter()
{
  if(pt_stack_.size() != 1) {
    throw ProgrammingError(
        "mismatched transparent contexts in XMLWriter: destructor reached with missing ends",
        __FILE__,
        __LINE__
    );
  }

  if(not writing_done_)
    do_write();

  if(delete_out_){
    out_->flush();
    delete out_;
  }

  if(delete_pt_)
    delete current_root_;
}

void
XMLWriter::run()
{
  begin_writing_context("mpqc");
  std::vector< Ref<XMLWritable> >::iterator it;
  for(it = data_.begin(); it != data_.end(); ++it){
    ptree& tmp = *current_root_;
    (*it)->write_xml(tmp, *this);
  }
  end_writing_context();
  do_write();
}

void
XMLWriter::do_write()
{
  if(writing_done_) {
    throw ProgrammingError(
        "already wrote xml to output",
        __FILE__,
        __LINE__
    );
  }

  Ref<MessageGrp> msg = MessageGrp::get_default_messagegrp();
  if(msg->me() != 0) return;
  boost::property_tree::write_xml(*out_, *current_root_, write_settings_);
  writing_done_ = true;

}

void
XMLWriter::begin_writing_context(const std::string& root_name)
{
  assert(not writing_done_);
  ptree& tmp_root = current_root_->add_child(root_name, ptree());
  current_root_ = &tmp_root;
  pt_stack_.push(current_root_);
}

void
XMLWriter::end_writing_context()
{
  assert(not writing_done_);
  if(pt_stack_.size() == 1) {
    throw ProgrammingError("mismatched transparent contexts in XMLWriter: too many ends", __FILE__, __LINE__);
  }
  pt_stack_.pop();
  current_root_ = pt_stack_.top();
}

////////////////////////////////////////////////////////////////////////////////
// Non-intrusive interface

ptree&
XMLWriter::write_to_xml(const Ref<XMLWritable>& obj, ptree& parent) const
{
  return obj->write_xml(parent, *this);
}

ptree&
XMLWriter::write_to_xml(XMLWritable& obj, ptree& parent) const
{
  return obj.write_xml(parent, *this);
}

ptree&
XMLWriter::write_to_xml(const Eigen::VectorXd& obj, ptree& parent) const
{
  ptree& child = parent.add_child("EigenVectorXd", ptree());
  const int n = obj.innerSize();
  child.put("<xmlattr>.n", n);
  // For now, just iterate over everything, since the data() method doesn't
  //   seem to work like I expect it to.
  double* data = allocate<double>(n);
  for(int i = 0; i < n; ++i){
    data[i] = obj(i);
  }
  // Note: the XMLDataStream created by put_binary_data now owns
  //   the pointer 'data'
  this->put_binary_data(child.add_child("data", ptree()), data, n);
  return child;
}

ptree&
XMLWriter::write_to_xml(const Eigen::MatrixXd& obj, ptree& parent) const
{
  ptree& child = parent.add_child("EigenMatrixXd", ptree());
  const int nrow = obj.rows();
  const int ncol = obj.cols();
  child.put("<xmlattr>.nrow", nrow);
  child.put("<xmlattr>.ncol", ncol);
  child.put("<xmlattr>.row_major", true);
  // For now, just iterate over everything, since the data() method doesn't
  //   seem to work like I expect it to.
  double* data = allocate<double>(nrow*ncol);
  int spot = 0;
  for(int row = 0; row < nrow; ++row){
    for(int col = 0; col < ncol; ++col, ++spot){
      data[spot] = obj(row, col);
    }
  }
  // Note: the XMLDataStream created by put_binary_data now owns
  //   the pointer 'data'
  this->put_binary_data(child.add_child("data", ptree()), data, nrow*ncol);
  return child;
}

ptree&
XMLWriter::write_to_xml(const SCVector3& obj, ptree& parent) const
{
  ptree& child = parent.add_child("SCVector3", ptree());
  // For now, just iterate over everything, since the data() method doesn't
  //   seem to work like I expect it to.
  double* data = allocate<double>(3);
  ::memcpy(data, obj.data(), 3*sizeof(double));
  this->put_binary_data(child.add_child("data", ptree()), data, 3);
  return child;
}

ptree&
XMLWriter::write_to_xml(const SCVector& obj, ptree& parent) const
{
  ptree& my_tree = parent.add_child("SCVector", ptree());
  my_tree.put("<xmlattr>.n", obj.n());
  ptree& data_tree = my_tree.add_child("data", ptree());
  double* data = allocate<double>(obj.n());
  obj.convert(data);
  // The XMLDataStream object created by this function call
  //   owns the pointer data after this.
  this->put_binary_data<double>(data_tree, data, obj.n());
  return my_tree;
}

ptree&
XMLWriter::write_to_xml(const Units& obj, ptree& parent) const
{
  ptree& child = parent.add_child("Units", ptree());
  child.put_value(obj.string_rep());
  return child;
}

ptree&
XMLWriter::write_to_xml(const Grid& obj, ptree& parent) const
{
  ptree& my_tree = parent.add_child("Grid", ptree());
  this->insert_child(my_tree, obj.unit);
  my_tree.put("num_x", obj.numx);
  my_tree.put("num_y", obj.numy);
  my_tree.put("num_z", obj.numz);
  this->insert_child(my_tree, obj.origin, "origin");
  this->insert_child(my_tree, obj.axisx, "x_axis");
  this->insert_child(my_tree, obj.axisy, "y_axis");
  this->insert_child(my_tree, obj.axisz, "z_axis");
  return my_tree;
}


