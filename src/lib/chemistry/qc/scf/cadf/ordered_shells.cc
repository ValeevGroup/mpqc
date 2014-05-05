//
// ordered_shells.cc
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: May 5, 2014
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

#include "ordered_shells.h"

using namespace sc;

////////////////////////////////////////////////////////////////////////////////

std::ostream&
sc::operator <<(std::ostream& out, const OrderedShellList& list) {
  out << "\nOrderedShellList: {";
  for(auto ish : list) {
    out << "\n  " << ish << ", with value " << ish.value;
  }
  out << "\n  aux_value = " << list.get_aux_value();
  if(list.aux_vector_initialized_) {
    out << "\n  aux_vector = { ";
    for(int ival = 0; ival < list.aux_vector_.rows(); ++ival) {
      out << list.aux_vector_[ival] << " ";
    }
    out << "}" << std::endl;
  }
  out << "\n}" << std::endl;
  return out;
}

