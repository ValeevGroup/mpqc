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
  clear();
}

void
DerivCenters::clear()
{
  ncenter_ = 0;
  omitted_center_ = -1;
  omitted_atom_ = -1;
}

void
DerivCenters::add_center(int center, int atom)
{
  center_[ncenter_] = center;
  atom_[ncenter_] = atom;
  ncenter_++;
}

void
DerivCenters::add_omitted(int center, int atom)
{
  omitted_center_ = center;
  omitted_atom_ = atom;
}

void
DerivCenters::add_center(int center, const Ref<GaussianBasisSet> &b, int shell)
{
  add_center(center, b->shell_to_center(shell));
}

void
DerivCenters::add_omitted(int center, const Ref<GaussianBasisSet> &b, int shell)
{
  add_omitted(center, b->shell_to_center(shell));
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
