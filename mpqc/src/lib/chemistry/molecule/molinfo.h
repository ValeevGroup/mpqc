//
// molinfo.h
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

#ifndef _chemistry_molecule_molinfo_h
#define _chemistry_molecule_molinfo_h

#include <util/class/class.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/chemelem.h>

class MolInfo: public DescribedClass {
#   define CLASSNAME MolInfo
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  protected:
    RefKeyVal keyval;
  public:
    MolInfo(const RefKeyVal&);
    virtual ~MolInfo();
    double doublevalue(const char * sym, const char * property);
    double doublevalue(const char * sym, const char * property, int i);
    KeyVal::KeyValError error();
};
DescribedClass_REF_dec(MolInfo);

#define ATOMINFO_MAXZ 100
class AtomInfo: public MolInfo {
#   define CLASSNAME AtomInfo
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    // certain values are cached here for fast access
    double radius_vals[ATOMINFO_MAXZ];
    int have_radius[ATOMINFO_MAXZ];
    double rgb_vals[ATOMINFO_MAXZ][3];
    int have_rgb[ATOMINFO_MAXZ];
    int get_zindex(const ChemicalElement&);

    // these items are cached for quick lookup
    double radius_scale_factor_;
    double radius_to_bohr_;

    void preload_values();
  public:
    AtomInfo(const RefKeyVal&);
    ~AtomInfo();
    double doublevalue(const ChemicalElement& atom,
                       const char * property);
    double doublevalue(const ChemicalElement& atom,
                       const char * property, int i);

    // routines for common properties
    double radius(const ChemicalElement&);
    double rgb(const ChemicalElement&, int color);
    double red(const ChemicalElement&);
    double green(const ChemicalElement&);
    double blue(const ChemicalElement&);
};
DescribedClass_REF_dec(AtomInfo);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
