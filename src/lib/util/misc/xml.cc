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

XMLWriter::XMLWriter(ostream& out) :
    out_(&out), pt_(ptree())
{

}

XMLWriter::XMLWriter(ptree& pt, ostream& out) :
    out_(&out), pt_(pt)
{

}

XMLWriter::XMLWriter(string filename) :
    pt_(ptree())
{
  if (filename == "-") {
      out_ = &(ExEnv::out0());
  }
  else {
      out_ = new std::ofstream(filename.c_str());
      // TODO Delete this
  }
}

