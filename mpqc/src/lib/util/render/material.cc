//
// material.cc
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

#include <util/render/material.h>

using namespace std;
using namespace sc;

static ClassDesc Material_cd(
  typeid(Material),"Material",1,"public DescribedClass",
  0, create<Material>, 0);

Material::Material()
{
  diffuse_.set(Color(0.5, 0.5, 0.5));
  ambient_.set(Color(0.5, 0.5, 0.5));
}

Material::Material(const Ref<KeyVal>& keyval)
{
  if (keyval->exists("diffuse")) {
      Color c(new PrefixKeyVal(keyval, "diffuse"));
      diffuse_.set(c);
    }
  if (keyval->exists("ambient")) {
      Color c(new PrefixKeyVal(keyval, "ambient"));
      ambient_.set(c);
    }
}

Material::~Material()
{
}

void
Material::print(ostream& os) const
{
  os << "Material" << endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
