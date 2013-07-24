//
// dercent.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
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

#include <chemistry/qc/basis/dercent.h>

using namespace sc;

DerivCenters::DerivCenters()
{
  centers_.reserve(4); // at most have 4-center integrals
  clear();
}

void
DerivCenters::clear()
{
  centers_.resize(0);
  omitted_center_ = std::make_pair(-1,-1);
}

void
DerivCenters::add_center(int center, int atom)
{
  centers_.push_back(std::make_pair(center, atom));
}

void
DerivCenters::add_omitted(int center, int atom)
{
  omitted_center_ = std::make_pair(center, atom);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
