//
// appearance.cc
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

#include <util/render/appearance.h>

using namespace std;

DescribedClass_REF_def(Appearance);
#define CLASSNAME Appearance
#define HAVE_KEYVAL_CTOR
#define PARENTS public DescribedClass
#include <util/class/classi.h>
void *
Appearance::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

Appearance::Appearance()
{
  level_.set(1);
}

Appearance::Appearance(const RefKeyVal& keyval)
{
  int level = keyval->intvalue("level");
  if (keyval->error() == KeyVal::OK) level_.set(level);
}

Appearance::~Appearance()
{
}

void
Appearance::print(ostream& os) const
{
  os << "Appearance:" << endl;
  os << "  level is ";
  if (level_.is_set()) {
      os << "set to " << level_.value();
      if (level_.overrides()) {
          os << " and overrides";
        }
    }
  else
      os << "not set";
  os << endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
